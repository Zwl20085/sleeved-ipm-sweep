"""Mechanical load-case set: result columns, worst-case envelope and lift-off speed.

Pure functions on the result dicts of :func:`sipm.mapdl.solve_mech` (one dict per load case);
the solver side is :func:`sipm.mapdl.solve_mech_cases`. See docs/methods.md section 3.1.

Columns written for a point with ``mech_cases`` set:

* ``<output>__<case>`` for every output of a single-case point and every requested case;
* the un-suffixed ``<output>`` columns hold the case named in ``mech_case``: ``rated_hot``
  (= ``mech_rpm`` at ``mech_temperature``, what they hold without ``mech_cases``) when it is
  requested, otherwise the first case of the list;
* ``<output>__worst`` and ``<output>__worst_case`` (the case that governs) for the outputs of
  :data:`ENVELOPE`; ties go to the earlier case of :data:`sipm.params.MECH_CASES`, so the
  envelope does not depend on the order of the list;
* ``liftoff_rpm_cold`` / ``liftoff_rpm_hot`` with ``liftoff_basis_cold`` / ``_hot`` (the two
  cases used, or why there is no value): see :func:`liftoff_at_temperature`; for a rotor whose
  contact is not linear in speed (:func:`search_applies`) the value comes from
  :func:`liftoff_search`, run by the solver with extra solutions at trial speeds.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence

import numpy as np

from .params import CASE_SEPARATOR, MECH_CASES, Params

# outputs with a worst-case column and the direction that is worse
ENVELOPE = (
    ("svm_lam_max", "max"),
    ("sh_sleeve_max", "max"),
    ("sr_sleeve_max", "max"),
    ("sleeve_fi_max_stress", "max"),
    ("sleeve_fi_tsai_wu", "max"),
    ("sleeve_fi_hashin_fibre", "max"),
    ("sleeve_fi_hashin_matrix", "max"),
    ("sleeve_sf", "min"),
    ("sleeve_open_fraction", "max"),
    ("contact_p_min", "min"),
    ("tau_sleeve_interface_max", "max"),   # layered sleeves only (absent outputs are skipped)
)
WORST = "worst"
CLOSED_STATUS = 2.0      # CONT,STAT of CONTA172: 0/1 open, 2 sliding, 3 sticking
CENTROID_TOL = 1e-6      # mm, element centroids of two runs of the same mesh must agree
# Largest deviation of a middle closed case's element pressure from the two-point line, relative
# to that element's standstill pressure p0_e, accepted as linear. Measured (docs/validation.md
# section 11): 2.8e-5 at a closed middle speed and 4e-4 at 0.95 x lift-off in the bridge-retained
# rotor; 6.6 % in the bridgeless rotor, where the linearity fails. 1 % lies 25x above the former
# and about 7x below the latter.
LINEARITY_TOL = 1e-2
NOT_LINEAR_TOPOLOGIES = {
    "sleeved_bridgeless": "the sleeve contact of sleeved_bridgeless is not linear in speed (loose lamination "
                          "pieces, partly open edges at the OD openings; docs/validation.md section 11)",
}


def case_column(output: str, case: str) -> str:
    """Column of ``output`` in load case ``case``: ``<output>__<case>``."""
    return f"{output}{CASE_SEPARATOR}{case}"


def worst_column(output: str) -> str:
    return f"{output}{CASE_SEPARATOR}{WORST}"


def worst_case_column(output: str) -> str:
    return f"{output}{CASE_SEPARATOR}{WORST}_case"


def base_case(cases: Sequence[str]) -> str:
    """The case whose values fill the un-suffixed columns: rated_hot if requested, else the first."""
    if not cases:
        raise ValueError("no load cases")
    return "rated_hot" if "rated_hot" in cases else cases[0]


def canonical(cases: Sequence[str]) -> list[str]:
    """``cases`` in the order of :data:`MECH_CASES`."""
    return [c for c in MECH_CASES if c in cases]


def _number(v) -> float | None:
    if isinstance(v, bool) or not isinstance(v, (int, float, np.integer, np.floating)):
        return None
    v = float(v)
    return None if math.isnan(v) else v


def envelope(results: Mapping[str, Mapping]) -> dict:
    """Worst value over the cases of each :data:`ENVELOPE` output and the case that governs.

    ``results`` maps case name -> result dict. Cases are visited in canonical order and a later
    case replaces the current worst only if it is strictly worse, so ties go to the earlier
    canonical case. Missing or NaN values are skipped; an output with no value is left out."""
    out = {}
    for name, how in ENVELOPE:
        best, where = None, None
        for case in canonical(list(results)):
            v = _number(results[case].get(name))
            if v is None:
                continue
            if best is None or (v > best if how == "max" else v < best):
                best, where = v, case
        if best is not None:
            out[worst_column(name)] = best
            out[worst_case_column(name)] = where
    return out


# ----------------------------------------------------------------------------- lift-off
def fully_closed(elems: Mapping[str, np.ndarray]) -> bool:
    """Every sleeve contact element closed (status sliding or sticking) with positive pressure."""
    return bool(len(elems["pres"]) and np.all(elems["pres"] > 0) and np.all(elems["status"] >= CLOSED_STATUS))


def liftoff_from_elements(s1: float, e1: Mapping[str, np.ndarray], s2: float, e2: Mapping[str, np.ndarray]
                          ) -> tuple[float, int]:
    """Speed at which the first contact element's pressure reaches zero, from two solutions.

    ``s1`` != ``s2`` are the speeds (any unit, e.g. rpm; the result is in the same unit) of two
    solutions at the same temperature, ``e1``, ``e2`` their per-element contact arrays
    (:func:`sipm.mapdl.parse_contact_elements`; same mesh, matched by element number and
    checked by centroid). Each element's pressure is taken as linear in the squared speed,
    p_e(s) = p0_e - c_e s^2, through the two points (exact while the contact stays closed and
    frictionless and the magnets are bonded: the model is then linear), and the result is
    min over the elements with c_e > 0 of sqrt(p0_e / c_e): element by element, because the
    element with the lowest pressure can change with speed. Returns (speed, element number);
    (inf, -1) if no element's pressure falls with speed, (0, element) if one is already <= 0
    at zero speed."""
    if s1 == s2:
        raise ValueError("lift-off extrapolation needs two different speeds")
    if len(e1["elem"]) != len(e2["elem"]) or np.any(e1["elem"] != e2["elem"]):
        raise ValueError("lift-off extrapolation: the two solutions have different contact elements")
    if np.max(np.hypot(e1["x"] - e2["x"], e1["y"] - e2["y"]), initial=0.0) > CENTROID_TOL:
        raise ValueError("lift-off extrapolation: contact element centroids differ between the two solutions")
    p1, p2 = np.asarray(e1["pres"], float), np.asarray(e2["pres"], float)
    c = (p1 - p2) / (s2 * s2 - s1 * s1)
    p0 = p1 + c * s1 * s1
    falling = c > 0
    if not falling.any():
        return math.inf, -1
    s2_lift = np.full(len(c), np.inf)
    s2_lift[falling] = np.maximum(p0[falling], 0.0) / c[falling]
    i = int(np.argmin(s2_lift))
    return math.sqrt(s2_lift[i]), int(e1["elem"][i])


def linearity_deviation(s1: float, e1: Mapping[str, np.ndarray], s2: float, e2: Mapping[str, np.ndarray],
                        s: float, e: Mapping[str, np.ndarray]) -> float:
    """Largest |p_e(s) - line_e(s)| / p0_e over the elements: how far the solution at speed ``s``
    lies from the per-element lines through (s1, e1) and (s2, e2) (same mesh, all closed), relative
    to each element's zero-speed pressure p0_e of the line."""
    if np.any(e["elem"] != e1["elem"]):
        raise ValueError("linearity check: the solutions have different contact elements")
    p1, p2 = np.asarray(e1["pres"], float), np.asarray(e2["pres"], float)
    c = (p1 - p2) / (s2 * s2 - s1 * s1)
    p0 = p1 + c * s1 * s1
    dev = np.abs(np.asarray(e["pres"], float) - (p0 - c * s * s)) / np.where(p0 > 0, p0, np.nan)
    return float(np.nanmax(dev)) if np.isfinite(dev).any() else math.inf


def liftoff_at_temperature(entries: Sequence[tuple[str, float, Mapping[str, np.ndarray]]]) -> tuple[float, str]:
    """Lift-off speed (rpm) from the load cases solved at one temperature.

    ``entries``: (case name, speed in rpm, contact element arrays) per case, all at the same
    temperature. Rules:

    * the lowest-speed case has an open element: 0 if it is at standstill (the interface is
      open before the rotor turns), else NaN (lift-off is below the lowest solved speed);
    * fewer than two fully closed cases at different speeds: NaN;
    * otherwise :func:`liftoff_from_elements` through the lowest and the highest fully closed
      speed (``inf`` if no element's pressure falls with speed).

    Returns (rpm, basis): the basis names the two cases used or the reason for 0 / NaN. Two
    checks of the closed-contact linearity are flagged in the basis (the value is still
    returned): a case that is open below the extrapolated speed, and, with three or more fully
    closed speeds, a middle one deviating from the two-point line by more than
    :data:`LINEARITY_TOL` (:func:`linearity_deviation`)."""
    rows, seen = [], set()
    for name, rpm, el in sorted(entries, key=lambda t: (t[1], MECH_CASES.index(t[0]) if t[0] in MECH_CASES else 0)):
        if rpm not in seen:
            seen.add(rpm)
            rows.append((name, float(rpm), el))
    if not rows:
        return math.nan, "no case at this temperature"
    name0, rpm0, el0 = rows[0]
    if not fully_closed(el0):
        if rpm0 == 0:
            return 0.0, f"open at {name0}"
        return math.nan, f"open at {name0}, the lowest solved speed"
    closed = [r for r in rows if fully_closed(r[2])]
    if len(closed) < 2:
        return math.nan, "needs two cases with closed contact at different speeds"
    (n1, s1, e1), (n2, s2, e2) = closed[0], closed[-1]
    rpm, _ = liftoff_from_elements(s1, e1, s2, e2)
    basis = f"{n1}+{n2}"
    worst = max(((linearity_deviation(s1, e1, s2, e2, s, e), n) for n, s, e in closed[1:-1]), default=None)
    if worst is not None and worst[0] > LINEARITY_TOL:
        basis += (f" (not linear: {worst[1]} deviates {100 * worst[0]:.2g} % from the line, tolerance "
                  f"{100 * LINEARITY_TOL:g} %)")
    below = [n for n, s, e in rows if not fully_closed(e) and s < rpm]
    if below:
        basis += f" (inconsistent: {', '.join(below)} open below this speed)"
    return rpm, basis


# ----------------------------------------------------------------------------- lift-off by search
SEARCH_MAX_SOLVES = 12   # extra MAPDL runs per temperature, upper bound of liftoff_search
SEARCH_GROWTH = 1.5      # factor between trial speeds while no bracket exists above the solved cases
# A contact whose pressure does not fall with speed (the inf of the two-point line) is checked with
# one trial at this multiple of the fastest solved case before the search returns inf. The fastest
# case is the overspeed (overspeed_factor x mech_rpm, default 1.2), so the check reaches 2.4 x rated
# speed and 4 x the centrifugal load of that case: beyond any operating or spin-test speed (spin
# tests run at about 1.2 x the maximum speed), so a lift-off above it has no design meaning.
SEARCH_NEVER_OPEN_FACTOR = 2.0


def search_applies(p: Params) -> bool:
    """Whether the lift-off speed of ``p`` is found by :func:`liftoff_search` (extra MAPDL runs)
    instead of the two-point extrapolation: a rotor of :data:`NOT_LINEAR_TOPOLOGIES` with
    frictionless sleeve contact and bonded magnets (as the linear path requires; with friction the
    columns stay NaN) and ``liftoff_search_tol`` > 0."""
    return (p.topology in NOT_LINEAR_TOPOLOGIES and p.mu_sleeve == 0 and p.magnet_contact == "bonded"
            and p.liftoff_search_tol > 0)


def _rows(entries):
    """(name, rpm, elems) sorted by speed (then canonical case order), one per speed."""
    rows, seen = [], set()
    for name, rpm, el in sorted(entries, key=lambda t: (t[1], MECH_CASES.index(t[0]) if t[0] in MECH_CASES else 0)):
        if rpm not in seen:
            seen.add(rpm)
            rows.append((name, float(rpm), el))
    return rows


def _quadratic_zero(points) -> float:
    """Lowest speed above the fastest of three closed solutions ``[(s, elems), ...]`` (ascending)
    at which an element's pressure, taken as quadratic in s^2 through its three values, reaches
    zero (the curvature of a contact that is not linear in speed; inf if none falls to zero)."""
    (s0, e0), (s1, e1), (s2, e2) = points
    if np.any(e0["elem"] != e2["elem"]) or np.any(e1["elem"] != e2["elem"]):
        return math.nan
    u0, u1, u2 = s0 * s0, s1 * s1, s2 * s2
    p0, p1, p2 = (np.asarray(e["pres"], float) for e in (e0, e1, e2))
    f01, f12 = (p1 - p0) / (u1 - u0), (p2 - p1) / (u2 - u1)
    a = (f12 - f01) / (u2 - u0)
    # Newton form p(u) = p2 + f12 (u - u2) + a (u - u2)(u - u1); in v = u - u2: p = c + bb v + a v^2
    c, bb = p2, f12 + a * (u2 - u1)
    roots = np.full(len(c), np.inf)
    du = u2 - u0
    lin = np.abs(a) * du * du <= 1e-9 * (np.abs(bb) * du + np.abs(c))
    m = lin & (bb < 0)
    roots[m] = -c[m] / bb[m]
    disc = bb * bb - 4 * a * c
    q = ~lin & (disc >= 0)
    sq = np.sqrt(np.where(q, disc, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        r1, r2 = (-bb - sq) / (2 * a), (-bb + sq) / (2 * a)
    for r in (r1, r2):
        ok = q & (r > 0)
        roots[ok] = np.minimum(roots[ok], r[ok])
    v = float(np.min(roots))
    return math.sqrt(u2 + v) if math.isfinite(v) else math.inf


def pressure_closed(elems: Mapping[str, np.ndarray]) -> bool:
    """Every sleeve contact element carries a positive contact pressure (``CONT,PRES`` > 0): the
    lift-off criterion of the two-point line (an element's pressure reaching zero), whatever the
    contact status of the element's integration points."""
    return bool(len(elems["pres"]) and np.all(np.asarray(elems["pres"], float) > 0))


def liftoff_search(entries: Sequence[tuple[str, float, Mapping[str, np.ndarray]]],
                   solve: Callable[[float], Mapping[str, np.ndarray]], tol: float,
                   max_solves: int = SEARCH_MAX_SOLVES) -> tuple[float, str, list[dict]]:
    """Lift-off speed (rpm) at one temperature, found directly with extra solves.

    The definition is that of the two-point line (:func:`liftoff_from_elements`): the lowest speed
    at which a contact element's pressure reaches zero (:func:`pressure_closed` turns false). Here
    it is bracketed and narrowed with solutions at trial speeds instead of read off a line, for
    contacts that are not linear in speed.

    ``entries`` as for :func:`liftoff_at_temperature` (the solved load cases at this temperature);
    ``solve(rpm)`` runs one more solution at that speed and temperature and returns its contact
    element arrays (it may raise; the search then stops with NaN). Rules:

    * the lowest-speed case not :func:`fully_closed`: 0 at standstill, NaN otherwise, exactly as on
      the linear path (no solves);
    * bracket: the solved speeds are walked upwards; ``lo`` is the last of the run of speeds with
      positive pressure everywhere from the lowest, ``hi`` the first speed after it where an
      element's pressure is zero. With only a closed standstill case there is no speed scale: NaN;
    * trial speeds: from the estimate ``est``, the element-by-element two-point extrapolation
      through the two fastest closed solutions known so far (locally accurate however non-linear
      the contact is over the whole range): ``est + 0.45 tol``, so that a good estimate is
      overshot by less than half the tolerance, and, once a zero-pressure speed ``hi`` lies within
      ``tol`` above it, ``est - 0.45 tol``, which closes the bracket; bisection when the estimate
      lies outside the bracket or the bracket did not halve over the last three trials; without
      an estimate and no zero-pressure speed yet, :data:`SEARCH_GROWTH` x the fastest closed
      speed;
    * a contact that does not open: while no zero-pressure speed is known and the two-point line
      gives ``inf`` both through the slowest and the fastest closed solution (the criterion of the
      linear path: no element's pressure lower at the faster one) and through the two fastest,
      and the quadratic in s^2 through the three fastest has no zero up to ``limit`` =
      :data:`SEARCH_NEVER_OPEN_FACTOR` x the fastest solved case (no element's pressure falls at
      the top of the range, and none curves down to zero below ``limit``), one trial runs at
      ``limit``. If the same holds with it, the result is ``inf`` (one extra solve); otherwise
      (open there, or falling) the search goes on as above;
    * stop when ``hi - lo <= tol``; the result is the midpoint, within ``tol / 2`` of the lift-off
      speed as resolved by the trials.

    At most ``max_solves`` extra solves; if the tolerance is not reached within them the result is
    NaN, with the last bracket in the basis. Like the two-point line, the search assumes one
    transition from closed to open in [lo, hi]; an opening that closes again between two closed
    speeds would not be seen. For ``inf`` it assumes in addition that a pressure which has not
    fallen at any element between the slowest case and ``limit`` (judged at the solved speeds and
    the trial at ``limit``) does not reach zero above ``limit``: a contact whose pressure first
    rises and then falls is found when the fall shows at the solved speeds or at ``limit`` (an
    element's pressure lower at the faster of the two fastest closed solutions, a zero of the
    quadratic through three below ``limit``, or the contact open at ``limit``); a fall that does
    not show there reads as ``inf``. Returns (rpm, basis, trials), trials = one dict per extra solve
    (rpm, closed = :func:`pressure_closed`, p_min, n_open = elements with zero pressure or an
    open status)."""
    rows = _rows(entries)
    trials: list[dict] = []
    if not rows:
        return math.nan, "no case at this temperature", trials
    name0, rpm0, el0 = rows[0]
    if not fully_closed(el0):
        if rpm0 == 0:
            return 0.0, f"open at {name0}", trials
        return math.nan, f"open at {name0}, the lowest solved speed", trials
    if not tol > 0:
        raise ValueError(f"lift-off search tolerance must be > 0 rpm, got {tol}")
    closed = [(rpm0, el0)]   # solutions with positive pressure everywhere, by speed
    hi = None
    for _name, rpm, el in rows[1:]:
        if pressure_closed(el):
            closed.append((rpm, el))
        else:
            hi = rpm
            break
    trying = [math.nan]   # speed of the solve in progress (for the message if it fails)

    def trial(s: float) -> bool:
        trying[0] = s
        el = solve(s)
        trying[0] = math.nan   # the run itself succeeded; a later error is not its failure
        ok = pressure_closed(el)
        pres, st = np.asarray(el["pres"], float), np.asarray(el["status"], float)
        trials.append(dict(rpm=s, closed=ok, p_min=float(pres.min()) if len(pres) else math.nan,
                           n_open=int(np.sum((st < CLOSED_STATUS) | (pres <= 0)))))
        if ok:
            closed.append((s, el))
            closed.sort(key=lambda t: t[0])
        return ok

    def estimate() -> float:
        if len(closed) < 2:
            return math.nan
        (s1, e1), (s2, e2) = closed[-2], closed[-1]
        s, _ = liftoff_from_elements(s1, e1, s2, e2)   # also checks that the meshes agree
        if len(closed) >= 3:
            q = _quadratic_zero(closed[-3:])
            if math.isfinite(q) and q > s2:
                return q
        return s

    def stopped(why: str) -> tuple[float, str, list[dict]]:
        return math.nan, f"not found: search stopped after {len(trials)} extra solve(s), {why}", trials

    limit = SEARCH_NEVER_OPEN_FACTOR * rows[-1][1]   # check speed of a contact that does not open

    def never_falls() -> bool:
        # the inf of the two-point line through the slowest and the fastest closed solution (the criterion
        # of the linear path) and through the two fastest, and no zero of the quadratic through the three
        # fastest below the limit (a zero beyond it, e.g. from the noise of a flat contact, is not chased)
        if any(not math.isinf(liftoff_from_elements(sa, ea, sb, eb)[0])
               for (sa, ea), (sb, eb) in ((closed[0], closed[-1]), (closed[-2], closed[-1]))):
            return False
        return len(closed) < 3 or not _quadratic_zero(closed[-3:]) <= limit

    widths: list[float] = []
    try:
        if hi is None and closed[-1][0] == 0:
            return math.nan, "needs a case with closed contact at a non-zero speed", trials
        while hi is None or hi - closed[-1][0] > tol:
            lo = closed[-1][0]
            flat = hi is None and len(closed) >= 2 and never_falls()
            if flat and lo >= limit:
                return math.inf, (
                    f"inf: contact pressure does not fall with speed (no element's pressure lower at {lo:.1f} "
                    f"rpm than at {closed[0][0]:.1f} rpm or falling between the two fastest closed solutions), "
                    f"positive everywhere up to {lo:.1f} rpm (check speed {limit:.1f} rpm = "
                    f"{SEARCH_NEVER_OPEN_FACTOR:g} x the fastest case), "
                    f"from the cases {'+'.join(r[0] for r in rows)} and {len(trials)} extra solve(s)"), trials
            if len(trials) >= max_solves:
                return stopped(f"still closed at {lo:.1f} rpm" if hi is None else
                               f"closed at {lo:.1f}, zero pressure at {hi:.1f} rpm, tolerance {tol:g} rpm not reached")
            est = estimate()
            if flat:                  # one check at the limit speed before inf
                s = limit
            elif hi is None:
                s = est + 0.45 * tol if math.isfinite(est) and est > lo else SEARCH_GROWTH * lo
                s = max(s, lo + tol)
            else:
                widths.append(hi - lo)
                stalled = len(widths) >= 4 and widths[-1] > 0.5 * widths[-4]
                if stalled or not (math.isfinite(est) and lo < est < hi):
                    s = 0.5 * (lo + hi)
                elif hi - est <= tol:     # zero pressure just above the estimate: close the bracket below it
                    s = max(est - 0.45 * tol, lo + 0.05 * tol)
                else:                     # aim just beyond the estimate, so that a good one ends the search
                    s = min(est + 0.45 * tol, hi - 0.05 * tol)
            if not trial(s):
                hi = s if hi is None else min(hi, s)
    except Exception as e:   # a failed trial run (or an inconsistent one): no value, the point itself stands
        if math.isnan(trying[0]):
            return stopped(f"{type(e).__name__}: {e}")
        return stopped(f"the run at {trying[0]:.1f} rpm failed ({type(e).__name__}: {e})")
    lo = closed[-1][0]
    basis = (f"found by iteration: pressure positive at {lo:.1f} rpm, zero at {hi:.1f} rpm, from the cases "
             f"{'+'.join(r[0] for r in rows)} and {len(trials)} extra solve(s) at trial speeds (tolerance "
             f"{tol:g} rpm)")
    return 0.5 * (lo + hi), basis, trials


def nonlinear_reason(p: Params) -> str:
    """Why the closed-contact model is not linear in the loads ('' if it is): the lift-off
    extrapolation needs frictionless sleeve contact and bonded magnets, and a rotor whose contact
    is linear in speed (not :data:`NOT_LINEAR_TOPOLOGIES`)."""
    why = []
    if p.topology in NOT_LINEAR_TOPOLOGIES:
        why.append(NOT_LINEAR_TOPOLOGIES[p.topology])
    if p.mu_sleeve != 0:
        why.append(f"mu_sleeve = {p.mu_sleeve:g}")
    if p.magnet_contact != "bonded":
        why.append(f"magnet_contact = {p.magnet_contact}")
    return ", ".join(why)


def liftoff_entries(p: Params, results: Mapping[str, Mapping], temperature: float) -> list:
    """(case, rpm, contact element arrays) of the cases solved at ``temperature``, whatever their label."""
    entries = []
    for case, r in results.items():
        rpm_c, t_c = p.case_speed_temperature(case)
        if t_c == temperature and "_contact_elems" in r:
            entries.append((case, rpm_c, r["_contact_elems"]))
    return entries


def liftoff_outputs(p: Params, results: Mapping[str, Mapping], searched: Mapping[str, tuple] | None = None) -> dict:
    """``liftoff_rpm_cold`` / ``_hot`` and their ``liftoff_basis_*`` from the per-case results
    (each with ``_contact_elems``). The cases at the temperature of the label are used, whatever
    their own label: with mech_temperature = mech_temperature_cold both use all cases.

    ``searched`` (label -> (rpm, basis), from :func:`liftoff_search`, see
    :func:`sipm.mapdl.solve_mech_cases`) gives the values of a label directly; it is used for
    the labels it has and only for them, and only the solver passes it (for :func:`search_applies`
    points or when the search is forced)."""
    out = {}
    why = nonlinear_reason(p)
    for label, temp in (("cold", p.mech_temperature_cold), ("hot", p.mech_temperature)):
        if searched and label in searched:
            rpm, basis = searched[label]
        elif why:
            rpm, basis = math.nan, f"not computed: model not linear ({why})"
        else:
            rpm, basis = liftoff_at_temperature(liftoff_entries(p, results, temp))
        out[f"liftoff_rpm_{label}"] = rpm
        out[f"liftoff_basis_{label}"] = basis
    return out


def combine(p: Params, results: Mapping[str, Mapping], wall_s: float | None = None,
            searched: Mapping[str, tuple] | None = None) -> dict:
    """One flat result dict for a point from its per-case result dicts (case -> dict).

    Un-suffixed columns from :func:`base_case` (named in ``mech_case``), every public output of
    every case as ``<output>__<case>``, the :func:`envelope` and the :func:`liftoff_outputs`
    (``searched``: lift-off values found by :func:`liftoff_search`, passed on to it).
    ``mech_wall_s`` becomes ``wall_s`` (the whole point) if given. Private entries: ``_raw``,
    ``_sleeve_nodes``, ``_contact_elems`` of the base case and ``_cases`` (all per-case dicts)."""
    cases = p.mech_case_list
    missing = [c for c in cases if c not in results]
    if missing:
        raise ValueError(f"missing results for load case(s) {missing}")
    base = base_case(cases)
    out = {k: v for k, v in results[base].items() if not k.startswith("_")}
    out["mech_case"] = base
    for case in canonical(cases):
        out.update({case_column(k, case): v for k, v in results[case].items() if not k.startswith("_")})
    out.update(envelope({c: results[c] for c in cases}))
    if any("_contact_elems" in results[c] for c in cases):
        out.update(liftoff_outputs(p, {c: results[c] for c in cases}, searched))
    if wall_s is not None:
        out["mech_wall_s"] = wall_s
    for k in ("_raw", "_sleeve_nodes", "_contact_elems"):
        if k in results[base]:
            out[k] = results[base][k]
    out["_cases"] = dict(results)
    return out


__all__ = ["ENVELOPE", "case_column", "worst_column", "worst_case_column", "base_case", "canonical", "envelope",
           "fully_closed", "pressure_closed", "liftoff_from_elements", "linearity_deviation",
           "liftoff_at_temperature", "liftoff_outputs", "liftoff_entries", "nonlinear_reason",
           "search_applies", "liftoff_search", "combine"]
