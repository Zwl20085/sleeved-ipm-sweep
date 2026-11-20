import time

import pytest

from sipm.config import Config
from sipm.memory import MemoryGate, Watchdog, free_ram_gb, total_ram_gb, tree_rss_gb
from sipm.sweep import Sweep, grid


class FakeClock:
    def __init__(self):
        self.t = 0.0
        self.sleeps = []

    def clock(self):
        return self.t

    def sleep(self, s):
        self.sleeps.append(s)
        self.t += s


def gate_with(readings, timeout=60.0, poll=10.0, min_free=4.0):
    it = iter(readings)
    last = [None]

    def reader():
        try:
            last[0] = next(it)
        except StopIteration:
            pass
        return last[0]

    fc = FakeClock()
    g = MemoryGate(min_free_gb=min_free, timeout_s=timeout, poll_s=poll, reader=reader, total=lambda: 32.0,
                   sleep=fc.sleep, clock=fc.clock)
    return g, fc


class Solver:
    """Fake solver with one fake PID per point; sleeps so a watchdog can sample it."""

    def __init__(self, delay=0.0):
        self.calls = 0
        self.delay = delay
        self.current = None

    def pids(self):
        return [self.current] if self.current else []

    def __call__(self, p, wd):
        self.calls += 1
        self.current = 1000 + self.calls
        try:
            time.sleep(self.delay)
            return dict(svm_rotor_max=100.0 + p.interference)
        finally:
            self.current = None


def test_real_readers_work():
    assert 0 < free_ram_gb() <= total_ram_gb()
    import os

    assert tree_rss_gb([os.getpid()]) > 0


def test_gate_waits_then_proceeds():
    g, fc = gate_with([1.0, 2.0, 6.0])
    ok, free = g.wait()
    assert ok and free == 6.0
    assert fc.sleeps == [10.0, 10.0]


def test_gate_threshold_fraction():
    g, _ = gate_with([5.0], min_free=4.0)
    g.min_free_fraction = 0.25  # 25 % of 32 GB = 8 GB > 4 GB
    assert g.threshold_gb == 8.0
    assert g.wait()[0] is False


def test_gate_timeout_stops_sweep_cleanly_and_resume(tmp_path):
    cfg = Config(workdir=tmp_path)
    sw = Sweep("g", "mech", grid(sleeve_thickness=[1.0, 2.0], interference=[0.0, 0.1, 0.2]),
               topology="sleeved_bridged", cfg=cfg)
    g, fc = gate_with([10.0, 10.0, 1.0], timeout=30.0)  # two points pass, then RAM stays low
    s1 = Solver()
    df = sw.run(solver=s1, gate=g, progress=None)
    assert s1.calls == 2
    assert sw.stopped_reason and "resume" in sw.stopped_reason
    assert (df.status == "ok").sum() == 2 and (df.status == "pending").sum() == 4
    assert (df.status == "failed").sum() == 0
    assert fc.t >= 30.0
    # resume with enough RAM: only the four pending points are solved
    g2, _ = gate_with([10.0])
    s2 = Solver()
    df = sw.run(solver=s2, gate=g2, progress=None)
    assert s2.calls == 4 and (df.status == "ok").all() and sw.stopped_reason is None


def test_watchdog_cap_kills_and_records(tmp_path):
    cfg = Config(workdir=tmp_path)
    sw = Sweep("w", "mech", grid(interference=[0.0, 0.1]), topology="sleeved_bridged", cfg=cfg)
    killed = []
    rss = lambda pids: 20.0 if 1001 in pids else 1.5  # noqa: E731  first point's solver is too big
    make = lambda pids: Watchdog(pids, cap_gb=8.0, poll_s=0.02, rss=rss, kill=killed.append)  # noqa: E731
    g, _ = gate_with([10.0])
    df = sw.run(solver=Solver(delay=0.3), gate=g, watchdog=make, progress=None).sort_values("index")
    assert killed == [1001]
    st = df.set_index("interference")
    assert st.loc[0.0, "status"] == "failed" and st.loc[0.0, "message"].startswith("memory cap")
    assert st.loc[0.1, "status"] == "ok" and st.loc[0.1, "solver_rss_peak_gb"] == pytest.approx(1.5)


def test_watchdog_monitor_only_and_idle():
    calls = []
    wd = Watchdog(lambda: [], cap_gb=None, rss=lambda p: calls.append(p) or 99.0, kill=lambda p: None)
    wd.sample()
    assert calls == [] and wd.peak_gb == 0.0
    wd = Watchdog(lambda: [7], cap_gb=None, rss=lambda p: 99.0, kill=lambda p: pytest.fail("must not kill"))
    wd.sample()
    assert wd.peak_gb == 99.0 and not wd.tripped


def test_config_memory_validation():
    with pytest.raises(ValueError):
        Config(max_parallel=2)
    with pytest.raises(ValueError):
        Config(max_solver_ram_gb=0)
    assert Config().min_free_ram_gb == 4.0 and Config().aedt_restart_every == 10
