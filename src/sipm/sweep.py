"""Sweep runner: design table in, results table out.

* every point is identified by a hash of its fully resolved parameters and the
  physics; its result (or failure) is written atomically to
  ``<workdir>/<name>/cache/<key>.json`` as soon as it finishes;
* re-running the same sweep skips cached points, so an interrupted sweep resumes
  where it stopped (failed points are kept as failures unless ``retry_failed``); a cached
  point that a later solver rule changed (:func:`point_rules`) is solved again (logged per rule
  before solving, with the number of points that run the lift-off search); a point with the
  lift-off search switched off reuses the record of the same point at the default tolerance if
  that was solved without the search (:func:`search_off_fallback_key`);
* invalid designs (parameter or geometry errors) and solver failures are recorded
  with their message and never stop the sweep;
* EM points that share a Maxwell model (same parameters apart from ``gamma`` and
  ``current_rms``) are run back to back in one project.
"""

from __future__ import annotations

import dataclasses
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
# private solver outputs kept out of the point cache (raw solver keys, per-node arrays)
NOT_CACHED = {"_raw", "_sleeve_nodes", "_contact_elems", "_cases", "_liftoff_trials"}


def grid(**axes: Iterable) -> list[dict[str, Any]]:
    """Full-factorial design table: ``grid(sleeve_thickness=[1, 2], gamma=[120, 130])``."""
    names = list(axes)
    return [dict(zip(names, vals)) for vals in itertools.product(*(list(axes[n]) for n in names))]


# Parameters added after the cache format was fixed (thermal load case, composite sleeve and
# its strengths, load-case set, layered sleeve and winding prestress, sleeve conductivity, lift-off
# search tolerance). Each is left out of the
# point key while it sits at its default, so a point that does not use the new features hashes
# exactly as before and caches written by earlier
# versions stay valid. Such old cached mech points lack the new output columns (sleeve failure
# indices, sr_sleeve_max, tau_sleeve_max); re-run with a fresh sweep name to get them. Results that
# changed without a parameter (a solver rule) are handled by point_rules below, not by the key.
# A field added in future belongs here too.
KEY_OMITTED_AT_DEFAULT = (
    "lam_cte", "pm_cte", "sleeve_cte_hoop", "sleeve_cte_radial", "sleeve_cte_axial",
    "mech_temperature", "mech_ref_temperature",
    "sleeve_fibre", "sleeve_matrix", "sleeve_vf", "sleeve_micromechanics", "sleeve_wind_angle",
    "sleeve_Xt", "sleeve_Xc", "sleeve_Yt", "sleeve_Yc", "sleeve_S",
    "mech_cases", "overspeed_factor", "mech_temperature_cold",
    "sleeve_layers", "sleeve_layer_prestress",
    "sleeve_sigma_fibre", "sleeve_sigma_transverse",
    "liftoff_search_tol",
)
_FIELD_DEFAULTS = {f.name: f.default for f in dataclasses.fields(Params)}


def point_key(p: Params, physics: str) -> str:
    d = {k: v for k, v in p.to_dict().items()
         if not (k in KEY_OMITTED_AT_DEFAULT and v == _FIELD_DEFAULTS[k])}
    s = json.dumps(dict(p=d, physics=physics, schema=CACHE_SCHEMA), sort_keys=True)
    return hashlib.sha1(s.encode()).hexdigest()[:16]


# Solver rules that changed the result of some points after their key format was fixed, without a
# parameter to put in the key. A point a rule applies to is solved under it and its cache record lists
# the rule (``rules``); a cached record of such a point that does not list it was written by earlier
# code and is solved again (same key, same file). Points no rule applies to are untouched: their key,
# their record and its reuse are as before. Rule names are never reused.
RULE_BORE_ARC = "bore-arc-mesh"        # single-layer sleeve thinner than 3 x mech_mesh (docs/methods.md, 3)
RULE_LIFTOFF_SEARCH = "liftoff-search"   # lift-off found by iteration (sipm.loadcases.search_applies)


def point_rules(p: Params, physics: str) -> tuple[str, ...]:
    """The rules of :data:`RULE_BORE_ARC` / :data:`RULE_LIFTOFF_SEARCH` that apply to the point
    (``()`` for most points: every EM point, every mech point with a sleeve at least 3 x
    ``mech_mesh`` thick or a layered sleeve, and no lift-off search). From the parameters alone:
    no geometry is built, so it is cheap and cannot raise for a point whose model would fail (that
    point fails in its own run, as before)."""
    if physics != "mech":
        return ()
    from .loadcases import search_applies
    from .mapdl import bore_arc_rule

    rules = []
    if bore_arc_rule(p):
        rules.append(RULE_BORE_ARC)
    if p.mech_cases and search_applies(p):
        rules.append(RULE_LIFTOFF_SEARCH)
    return tuple(rules)


def is_current(record: dict | None, rules: Iterable[str]) -> bool:
    """A cache record is current if it exists and lists every rule that applies to its point."""
    return record is not None and set(rules) <= set(record.get("rules", ()))


def search_off_fallback_key(p: Params, physics: str) -> str | None:
    """Second cache key of a mech point with ``liftoff_search_tol`` = 0 (search switched off), or None.

    ``liftoff_search_tol`` is in the key except at its default, so a point with the search switched
    off hashes apart from the same point at the default, which hashes as before the search existed.
    The solver reads the tolerance only to decide and run the search; a record under the default key
    that does not list :data:`RULE_LIFTOFF_SEARCH` was therefore solved without the search (code
    before it, or a point it does not apply to) and is the result of this point too: :meth:`Sweep.run`
    and :meth:`Sweep.results` use it when the point's own key has no record and the record is
    current for the point's own rules (:func:`is_current`). A record that lists the rule holds
    searched lift-off values and is not used."""
    if physics != "mech" or p.liftoff_search_tol != 0:
        return None
    return point_key(p.replace(liftoff_search_tol=_FIELD_DEFAULTS["liftoff_search_tol"]), physics)


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

    def lookup(self, p: Params, key: str, rules: Iterable[str]) -> tuple[dict | None, str | None]:
        """The cache record of a point and the key it was found under: its own ``key``, or for a
        point with the lift-off search switched off and no record of its own, the record of the
        same point at the default tolerance if that was solved without the search
        (:func:`search_off_fallback_key`). (None, None) if neither."""
        c = self.load_cached(key)
        if c is not None:
            return c, key
        fk = search_off_fallback_key(p, self.physics)
        if fk is not None:
            f = self.load_cached(fk)
            if f is not None and RULE_LIFTOFF_SEARCH not in f.get("rules", ()) and is_current(f, rules):
                return f, fk
        return None, None

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
        todo, stale, borrowed = [], {}, 0
        for i, row, p, key, err in plan:
            if err is not None:
                continue
            rules = point_rules(p, self.physics)
            c, found = self.lookup(p, key, rules)
            if c is not None and not is_current(c, rules):
                for r in set(rules) - set(c.get("rules", ())):
                    stale[r] = stale.get(r, 0) + 1
            if not is_current(c, rules) or (retry_failed and c.get("status") != "ok"):
                todo.append((i, row, p, key, rules))
            elif found != key:
                borrowed += 1
        if self.physics in ("em", "acloss"):
            from .maxwell import em_model_key

            todo.sort(key=lambda t: (em_model_key(t[2]), t[0]))
        log(f"[{self.name}] {len(plan)} points, {len(todo)} to solve, "
            f"{sum(1 for t in plan if t[4])} invalid")
        for line in _cost_notes(todo, stale, borrowed):
            log(f"[{self.name}] {line}")
        solved = 0
        ctx = solver or _default_solver(self.physics, self.cfg)
        try:
            for i, row, p, key, rules in todo:
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
                if rules:   # written only where a rule applies: every other record is as before
                    rec["rules"] = list(rules)
                wd = make_wd(getattr(ctx, "pids", lambda: []))
                with wd:
                    try:
                        res = ctx(p, self.dir / "points" / key)
                        rec.update(status="ok", message="", outputs={k: v for k, v in res.items()
                                                                     if not k.startswith("_")},
                                   extra={k[1:]: v for k, v in res.items()
                                          if k.startswith("_") and k not in NOT_CACHED})
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
                rules = point_rules(p, self.physics)
                c, _found = self.lookup(p, key, rules)
                if c is None:
                    r.update(status="pending", message="")
                elif not is_current(c, rules):
                    missing = sorted(set(rules) - set(c.get("rules", ())))
                    r.update(status="pending", message=f"cached result predates solver rule(s) {', '.join(missing)}; "
                                                       "run the sweep again to re-solve it")
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


def _cost_notes(todo: list, stale: dict[str, int], borrowed: int) -> list[str]:
    """Log lines on the solver work that is not plain: points that run the lift-off search, cached
    points solved again because of a solver rule, records reused from the default tolerance."""
    from .loadcases import SEARCH_MAX_SOLVES

    lines = []
    n = sum(1 for t in todo if RULE_LIFTOFF_SEARCH in t[4])
    if n:
        lines.append(f"{n} of the {len(todo)} points to solve run the lift-off search (rule {RULE_LIFTOFF_SEARCH}): "
                     f"up to {SEARCH_MAX_SOLVES} extra MAPDL runs per temperature each, "
                     f"{2 * SEARCH_MAX_SOLVES} per point; liftoff_search_tol = 0 switches it off (lift-off NaN)")
    for rule in sorted(stale):
        hint = ("; with liftoff_search_tol = 0 a record that lacks only this rule is reused"
                if rule == RULE_LIFTOFF_SEARCH else "")
        lines.append(f"{stale[rule]} cached point(s) solved again: the record predates solver rule {rule}{hint}")
    if borrowed:
        lines.append(f"{borrowed} point(s) with liftoff_search_tol = 0 reuse the record of the same point at the "
                     "default tolerance, solved without the search")
    return lines


def _brief(o: dict) -> str:
    keys = [k for k in ("torque_avg", "sleeve_loss_w", "svm_rotor_max", "sh_sleeve_max", "kac") if k in o]
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
