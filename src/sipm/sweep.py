"""Sweep runner: design table in, results table out.

* every point is identified by a hash of its fully resolved parameters and the
  physics; its result (or failure) is written atomically to
  ``<workdir>/<name>/cache/<key>.json`` as soon as it finishes;
* re-running the same sweep skips cached points, so an interrupted sweep resumes
  where it stopped (failed points are kept as failures unless ``retry_failed``);
* invalid designs (parameter or geometry errors) and solver failures are recorded
  with their message and never stop the sweep;
* EM points that share a Maxwell model (same parameters apart from ``gamma`` and
  ``current_rms``) are run back to back in one project.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import time
import traceback
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from .config import Config
from .geometry import GeometryError, build_rotor, build_stator
from .memory import MemoryGate, Watchdog
from .params import Params, coerce, default_params, field_names

CACHE_SCHEMA = 1


def grid(**axes: Iterable) -> list[dict[str, Any]]:
    """Full-factorial design table: ``grid(sleeve_thickness=[1, 2], gamma=[120, 130])``."""
    names = list(axes)
    return [dict(zip(names, vals)) for vals in itertools.product(*(list(axes[n]) for n in names))]


def point_key(p: Params, physics: str) -> str:
    s = json.dumps(dict(p=p.to_dict(), physics=physics, schema=CACHE_SCHEMA), sort_keys=True)
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def _atomic_write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, default=_json_default)
    os.replace(tmp, path)


def _json_default(o):
    import numpy as np

    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def resolve(row: dict, base: dict | None = None, topology: str | None = None) -> Params:
    """Params for one table row. ``topology`` may come from the row, ``base`` or the argument."""
    merged = dict(base or {})
    merged.update({k: v for k, v in row.items() if not (isinstance(v, float) and v != v)})
    topo = merged.pop("topology", None) or topology or "sleeved_bridged"
    unknown = set(merged) - set(field_names())
    if unknown:
        raise KeyError(f"unknown parameter column(s): {sorted(unknown)}")
    merged = {k: coerce(k, v) for k, v in merged.items()}
    return default_params(topo, **merged)


class Sweep:
    """A named sweep in a working directory.

    Parameters
    ----------
    name : folder name under ``cfg.workdir``
    physics : 'em' (torque at gamma, current), 'mech' (stress) or 'acloss'
    table : list of dicts or a DataFrame; columns are parameter names (+ optional 'topology')
    base : overrides applied to every row before the row's own values
    """

    def __init__(self, name: str, physics: str, table, base: dict | None = None,
                 topology: str | None = None, cfg: Config | None = None):
        if physics not in ("em", "mech", "acloss"):
            raise ValueError("physics must be 'em', 'mech' or 'acloss'")
        self.name, self.physics = name, physics
        self.cfg = cfg or Config()
        if isinstance(table, pd.DataFrame):
            table = table.to_dict("records")
        self.rows = [dict(r) for r in table]
        self.base = dict(base or {})
        self.topology = topology
        self.dir = Path(self.cfg.workdir) / name
        self.cache = self.dir / "cache"

    # ------------------------------------------------------------------------------------
    def _cache_path(self, key: str) -> Path:
        return self.cache / f"{key}.json"

    def load_cached(self, key: str) -> dict | None:
        f = self._cache_path(key)
        if not f.exists():
            return None
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return None  # a torn write cannot happen with os.replace, but be tolerant

    def plan(self):
        """Resolve all rows: list of (index, row, Params|None, key|None, error|None)."""
        out = []
        for i, row in enumerate(self.rows):
            try:
                p = resolve(row, self.base, self.topology)
                build_rotor(p)
                if self.physics != "mech":
                    build_stator(p)
                out.append((i, row, p, point_key(p, self.physics), None))
            except (ValueError, KeyError, GeometryError) as e:
                out.append((i, row, None, None, f"{type(e).__name__}: {e}"))
        return out

    def run(self, solver: Callable | None = None, stop_after: int | None = None, retry_failed: bool = False,
            progress: Callable[[str], None] | None = print, gate: MemoryGate | None = None,
            watchdog: Callable[..., Watchdog] | None = None) -> pd.DataFrame:
        """Solve all uncached points and return the results table (also saved as ``results.csv``).

        ``stop_after`` stops after that many newly solved points (interrupt/resume demo).
        Before every point the memory gate waits for ``cfg.min_free_ram_gb`` of free RAM; if
        that does not happen within ``cfg.ram_wait_timeout_s`` the sweep stops cleanly
        (``self.stopped_reason`` is set, cached points are kept, re-run to resume). While a
        point runs, a watchdog records the peak RSS of the solver process tree and, if
        ``cfg.max_solver_ram_gb`` is set, kills that tree above it and records the point as
        failed with reason "memory cap". ``gate``/``watchdog`` override both (tests)."""
        log = progress or (lambda s: None)
        self.stopped_reason = None
        cf = self.cfg
        gate = gate or MemoryGate(min_free_gb=cf.min_free_ram_gb, min_free_fraction=cf.min_free_ram_fraction,
                                  timeout_s=cf.ram_wait_timeout_s, poll_s=cf.ram_poll_s)
        cap = cf.max_solver_ram_gb
        make_wd = watchdog or (lambda pids: Watchdog(pids, cap_gb=cap, poll_s=1.0))
        self.dir.mkdir(parents=True, exist_ok=True)
        plan = self.plan()
        todo = []
        for i, row, p, key, err in plan:
            if err is not None:
                continue
            c = self.load_cached(key)
            if c is None or (retry_failed and c.get("status") != "ok"):
                todo.append((i, row, p, key))
        if self.physics in ("em", "acloss"):
            from .maxwell import em_model_key

            todo.sort(key=lambda t: (em_model_key(t[2]), t[0]))
        log(f"[{self.name}] {len(plan)} points, {len(todo)} to solve, "
            f"{sum(1 for t in plan if t[4])} invalid")
        solved = 0
        ctx = solver or _default_solver(self.physics, self.cfg)
        try:
            for i, row, p, key in todo:
                if stop_after is not None and solved >= stop_after:
                    log(f"[{self.name}] stopping after {solved} new points (stop_after)")
                    break
                ok, free = gate.wait(log)
                if not ok:
                    self.stopped_reason = (
                        f"free RAM stayed below {gate.threshold_gb:.1f} GB for {gate.timeout_s:.0f} s "
                        f"(last reading {free:.1f} GB); sweep stopped with {solved} new points cached - "
                        "run the same sweep again to resume")
                    log(f"[{self.name}] {self.stopped_reason}")
                    break
                t0 = time.time()
                rec = dict(key=key, index=i, physics=self.physics, inputs=row, params=p.to_dict(),
                           started=time.strftime("%Y-%m-%d %H:%M:%S"), free_ram_gb_at_start=free)
                wd = make_wd(getattr(ctx, "pids", lambda: []))
                with wd:
                    try:
                        res = ctx(p, self.dir / "points" / key)
                        rec.update(status="ok", message="", outputs={k: v for k, v in res.items()
                                                                     if not k.startswith("_")},
                                   extra={k[1:]: v for k, v in res.items() if k.startswith("_") and k != "_raw"})
                    except KeyboardInterrupt:
                        raise
                    except Exception as e:
                        rec.update(status="failed", message=f"{type(e).__name__}: {e}", outputs={},
                                   traceback=traceback.format_exc(limit=6))
                if wd.tripped:
                    rec.update(status="failed", outputs={},
                               message=f"memory cap: solver tree reached {wd.tripped_at_gb:.2f} GB RSS "
                                       f"(> max_solver_ram_gb = {wd.cap_gb} GB) and was terminated")
                if rec["status"] != "ok":
                    _reset_solver(ctx)
                rec["solver_rss_peak_gb"] = wd.peak_gb
                if rec["status"] == "ok":
                    rec["outputs"]["solver_rss_peak_gb"] = wd.peak_gb
                rec["wall_s"] = time.time() - t0
                _atomic_write_json(self._cache_path(key), rec)
                solved += 1
                log(f"[{self.name}] point {i}: {rec['status']} in {rec['wall_s']:.1f} s "
                    + (_brief(rec["outputs"]) if rec["status"] == "ok" else rec["message"][:300]))
        finally:
            _close_solver(ctx)
        df = self.results()
        df.to_csv(self.dir / "results.csv", index=False)
        return df

    def results(self) -> pd.DataFrame:
        """Table with one row per design-table row: inputs, outputs, status, message."""
        recs = []
        for i, row, p, key, err in self.plan():
            r = {"index": i, **row}
            if p is not None:
                r.setdefault("topology", p.topology)
            if err is not None:
                r.update(status="invalid", message=err)
            else:
                c = self.load_cached(key)
                if c is None:
                    r.update(status="pending", message="")
                else:
                    r.update(c.get("outputs", {}))
                    r.update(status=c["status"], message=c.get("message", ""), wall_s=c.get("wall_s"))
                r["key"] = key
                for k in ("sleeve_thickness", "interference", "gamma", "pm1_bridge_thickness",
                          "pm2_bridge_thickness", "pm1_rib_half_width", "pm2_web_width", "lip_depth",
                          "current_rms", "mech_rpm", "em_rpm"):
                    r.setdefault(k, getattr(p, k))
            recs.append(r)
        return pd.DataFrame.from_records(recs)

    def waveforms(self) -> dict[int, dict]:
        """Per-row waveform data saved by EM points (index -> {quantity: list})."""
        out = {}
        for i, _row, _p, key, _err in self.plan():
            if key:
                c = self.load_cached(key)
                if c and c.get("status") == "ok" and "wave" in c.get("extra", {}):
                    out[i] = c["extra"]["wave"]
        return out


def _brief(o: dict) -> str:
    keys = [k for k in ("torque_avg", "svm_rotor_max", "sh_sleeve_max", "kac") if k in o]
    return " ".join(f"{k}={o[k]:.4g}" for k in keys)


# ----------------------------------------------------------------------------------------
class _MechSolver:
    def __init__(self, cfg):
        self.cfg = cfg
        self._pid = None

    def pids(self):
        return [self._pid] if self._pid else []

    def _started(self, pid):
        self._pid = pid

    def __call__(self, p, wd):
        from .mapdl import solve_mech

        try:
            return solve_mech(p, self.cfg, wd, on_start=self._started)
        finally:
            self._pid = None


class _EmSolver:
    """One AEDT session, released and restarted every ``cfg.aedt_restart_every`` points."""

    def __init__(self, cfg, physics):
        self.cfg, self.physics, self.session = cfg, physics, None
        self.count = 0

    def pids(self):
        return [self.session.pid] if self.session is not None and self.session.pid else []

    def __call__(self, p, wd):
        from .maxwell import MaxwellSession

        n = self.cfg.aedt_restart_every
        if self.session is not None and n and self.count >= n:
            self.close()  # release AEDT so its memory cannot accumulate over a long sweep
        if self.session is None:
            self.session = MaxwellSession(self.cfg)
            self.session.open()
            self.count = 0
        self.count += 1
        if self.physics == "em":
            return self.session.solve(p)
        return self.session.ac_loss(p, p.em_rpm)

    def reset(self):
        self.close()

    def close(self):
        if self.session is not None:
            self.session.close()
            self.session = None


def _default_solver(physics, cfg):
    return _MechSolver(cfg) if physics == "mech" else _EmSolver(cfg, physics)


def _reset_solver(s):
    if hasattr(s, "reset"):
        try:
            s.reset()
        except Exception:
            pass


def _close_solver(s):
    if hasattr(s, "close"):
        try:
            s.close()
        except Exception:
            pass


def mtpa(df: pd.DataFrame, design_cols: list[str] | None = None) -> pd.DataFrame:
    """Peak average torque over gamma for each design (rows with status 'ok')."""
    ok = df[df["status"] == "ok"].copy()
    if design_cols is None:
        from .post import GEOMETRY_SWEEP_COLUMNS

        design_cols = [c for c in ("topology", "sleeve_thickness") + GEOMETRY_SWEEP_COLUMNS + ("current_rms",)
                       if c in ok.columns]
    idx = ok.groupby(design_cols, dropna=False)["torque_avg"].idxmax()
    return ok.loc[idx].reset_index(drop=True)
