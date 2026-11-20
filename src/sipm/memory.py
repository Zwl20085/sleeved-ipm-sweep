"""Memory protection for sweeps (psutil, cross-platform).

* :class:`MemoryGate` - before each point, wait while free physical RAM is below a
  threshold; give up after a timeout so the sweep can stop cleanly (and be resumed).
* :class:`Watchdog` - while a point runs, sample the resident memory (RSS) of the solver
  process tree that this package started, record the peak, and terminate that tree if it
  exceeds a cap.

Both take injectable readers/sleepers so they can be unit tested without real memory
pressure.
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from typing import Callable, Iterable

GB = 1e9


def free_ram_gb() -> float:
    import psutil

    return psutil.virtual_memory().available / GB


def total_ram_gb() -> float:
    import psutil

    return psutil.virtual_memory().total / GB


def tree_rss_gb(pids: Iterable[int]) -> float:
    """Summed RSS of the given processes and all their children (each counted once)."""
    import psutil

    seen, total = set(), 0
    for pid in pids:
        try:
            root = psutil.Process(pid)
            procs = [root] + root.children(recursive=True)
        except psutil.Error:
            continue
        for p in procs:
            if p.pid in seen:
                continue
            seen.add(p.pid)
            try:
                total += p.memory_info().rss
            except psutil.Error:
                pass
    return total / GB


def kill_tree(pid: int) -> None:
    """Terminate a process and its children (only ever called on processes we started)."""
    import psutil

    try:
        root = psutil.Process(pid)
        procs = root.children(recursive=True) + [root]
    except psutil.Error:
        return
    for p in procs:
        try:
            p.kill()
        except psutil.Error:
            pass
    psutil.wait_procs(procs, timeout=10)


@dataclass
class MemoryGate:
    """Wait for free RAM before launching a point.

    The threshold is ``max(min_free_gb, min_free_fraction * total)``.
    """

    min_free_gb: float = 4.0
    min_free_fraction: float = 0.0
    timeout_s: float = 1800.0
    poll_s: float = 10.0
    reader: Callable[[], float] = field(default=free_ram_gb)
    total: Callable[[], float] = field(default=total_ram_gb)
    sleep: Callable[[float], None] = field(default=time.sleep)
    clock: Callable[[], float] = field(default=time.monotonic)

    @property
    def threshold_gb(self) -> float:
        frac = self.min_free_fraction * self.total() if self.min_free_fraction > 0 else 0.0
        return max(self.min_free_gb, frac)

    def wait(self, log: Callable[[str], None] | None = None) -> tuple[bool, float]:
        """Block until free RAM >= threshold. Returns (ok, free_gb); ok is False on timeout."""
        thr = self.threshold_gb
        t0 = self.clock()
        free = self.reader()
        warned = False
        while free < thr:
            if self.clock() - t0 >= self.timeout_s:
                return False, free
            if log and not warned:
                log(f"free RAM {free:.1f} GB < {thr:.1f} GB: waiting (timeout {self.timeout_s:.0f} s)")
                warned = True
            self.sleep(self.poll_s)
            free = self.reader()
        return True, free


class Watchdog:
    """Sample the RSS of the solver tree during a point; kill it above ``cap_gb``.

    ``pids`` returns the root PIDs this package started for the running point (it may be
    empty before the solver has been launched). ``cap_gb`` None = monitor only.
    """

    def __init__(self, pids: Callable[[], list[int]], cap_gb: float | None = None, poll_s: float = 1.0,
                 rss: Callable[[list[int]], float] = tree_rss_gb, kill: Callable[[int], None] = kill_tree):
        self.pids, self.cap_gb, self.poll_s = pids, cap_gb, poll_s
        self.rss, self.kill = rss, kill
        self.peak_gb = 0.0
        self.tripped = False
        self.tripped_at_gb = None
        self._stop = threading.Event()
        self._thread = None

    def sample(self) -> None:
        pids = [p for p in (self.pids() or []) if p]
        if not pids:
            return
        r = self.rss(pids)
        self.peak_gb = max(self.peak_gb, r)
        if self.cap_gb is not None and r > self.cap_gb and not self.tripped:
            self.tripped, self.tripped_at_gb = True, r
            for p in pids:
                self.kill(p)

    def _run(self):
        while not self._stop.is_set():
            try:
                self.sample()
            except Exception:
                pass
            self._stop.wait(self.poll_s)

    def __enter__(self):
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5 * self.poll_s + 1)
        return False
