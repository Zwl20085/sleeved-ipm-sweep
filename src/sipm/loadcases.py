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
  cases used, or why there is no value): see :func:`liftoff_at_temperature`.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

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


def liftoff_outputs(p: Params, results: Mapping[str, Mapping]) -> dict:
    """``liftoff_rpm_cold`` / ``_hot`` and their ``liftoff_basis_*`` from the per-case results
    (each with ``_contact_elems``). The cases at the temperature of the label are used, whatever
    their own label: with mech_temperature = mech_temperature_cold both use all cases."""
    out = {}
    why = nonlinear_reason(p)
    for label, temp in (("cold", p.mech_temperature_cold), ("hot", p.mech_temperature)):
        if why:
            rpm, basis = math.nan, f"not computed: model not linear ({why})"
        else:
            entries = []
            for case, r in results.items():
                rpm_c, t_c = p.case_speed_temperature(case)
                if t_c == temp and "_contact_elems" in r:
                    entries.append((case, rpm_c, r["_contact_elems"]))
            rpm, basis = liftoff_at_temperature(entries)
        out[f"liftoff_rpm_{label}"] = rpm
        out[f"liftoff_basis_{label}"] = basis
    return out


def combine(p: Params, results: Mapping[str, Mapping], wall_s: float | None = None) -> dict:
    """One flat result dict for a point from its per-case result dicts (case -> dict).

    Un-suffixed columns from :func:`base_case` (named in ``mech_case``), every public output of
    every case as ``<output>__<case>``, the :func:`envelope` and the :func:`liftoff_outputs`.
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
        out.update(liftoff_outputs(p, {c: results[c] for c in cases}))
    if wall_s is not None:
        out["mech_wall_s"] = wall_s
    for k in ("_raw", "_sleeve_nodes", "_contact_elems"):
        if k in results[base]:
            out[k] = results[base][k]
    out["_cases"] = dict(results)
    return out


__all__ = ["ENVELOPE", "case_column", "worst_column", "worst_case_column", "base_case", "canonical", "envelope",
           "fully_closed", "liftoff_from_elements", "linearity_deviation",
           "liftoff_at_temperature", "liftoff_outputs", "nonlinear_reason",
           "combine"]
