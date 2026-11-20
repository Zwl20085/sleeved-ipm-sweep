"""Feasibility checks and derived metrics of a built rotor pole.

Every failed check raises :class:`GeometryError` naming the violated constraint. The checks
that need only the construction points (negative or too short magnet, slot not reaching the
bridge circle, overlapping bridge bands / slot mouths, layer order) are made in
:func:`sipm.geometry.rotor.build_rotor` before the lamination is evaluated; the ones below
need the evaluated lamination:

* every magnet lies in air (it does not overlap steel);
* air regions that must stay separate do (upper/lower slot of a layer where a rib exists,
  layer 1 / layer 2, slots / cavity / holes), and the steel between them, and between them
  and the bore, is at least ``min_web``;
* the steel on the q-axis between adjacent poles is at least ``min_q_web``;
* bridgeless rotor: each OD opening lies inside its slot mouth (lips vs openings).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from .primitives import GeometryError, chain_polyline

if TYPE_CHECKING:  # pragma: no cover
    from .rotor import RotorPole


def _in_chain(chain, pts: np.ndarray) -> np.ndarray:
    from matplotlib.path import Path

    P = chain_polyline(chain, 0.25, 48)
    return Path(np.vstack([P, P[:1]]), closed=True).contains_points(pts)


def in_steel(rp: RotorPole, pts: np.ndarray) -> np.ndarray:
    """True where a point lies in the lamination (sampled boundaries, 0.25 deg / 48 per curve)."""
    pts = np.atleast_2d(pts)
    out = np.zeros(len(pts), bool)
    for f in rp.lamination:
        m = _in_chain(f.outer, pts)
        for h in f.holes:
            m &= ~_in_chain(h, pts)
        out |= m
    return out


def _dist(P: np.ndarray, Q: np.ndarray) -> float:
    """Minimum distance between two closed polylines (0 if they cross)."""
    from .polygons import polyline_distance

    return polyline_distance(P, Q)


def _air_regions(rp: RotorPole) -> list[tuple[str, np.ndarray]]:
    """Hole loops of the lamination with the name of the feature they contain."""
    out = []
    for f in rp.lamination:
        for h in f.holes:
            P = chain_polyline(h, 0.25, 48)
            out.append((_name_hole(rp, h), P))
    return out


def _name_hole(rp: RotorPole, h) -> str:
    names = []
    for m in rp.magnets:
        if _in_chain(h, m.centre[None, :])[0]:
            names.append(f"layer {m.layer} {'upper' if m.side > 0 else 'lower'} slot")
    for nm, q in rp.features.items():
        if _in_chain(h, np.array([q]))[0]:
            names.append(nm.replace("_", " "))
    d_hole = np.array([[rp.r_lam - rp.p.d_hole_depth, 0.0]])
    if _in_chain(h, d_hole)[0]:
        names.append("d-axis hole")
    return " + ".join(names) if names else "air region"


def check_rotor(rp: RotorPole) -> dict:
    """Measure clearances and derived quantities; raise GeometryError on any infeasibility."""
    p = rp.p
    hp = rp.half_pitch
    m: dict = {}
    # ---- magnets in air --------------------------------------------------------------
    for mg in rp.magnets:
        P = chain_polyline(mg.chain, 2.0, 16)
        c = mg.centre
        probe = np.vstack([c[None, :], P + 0.02 * (c - P) / np.linalg.norm(c - P, axis=1)[:, None]])
        if in_steel(rp, probe).any():
            raise GeometryError(f"layer {mg.layer} magnet overlaps the lamination (it is not inside its slot)")
    # ---- bridgeless: lips vs openings ------------------------------------------------------
    if p.is_bridgeless:
        from .rotor import opening_angles

        r_in = rp.r_lam - p.lip_depth
        for L in rp.layers:
            lo, hi = opening_angles(p, L)
            a1 = _wall_angle(L, r_in, +1)
            a2 = _wall_angle(L, r_in, -1)
            if not (a1 - 1e-12 <= lo and hi <= a2 + 1e-12):
                raise GeometryError(f"layer {L.inp.k}: OD opening {math.degrees(lo):.3f}..{math.degrees(hi):.3f} deg "
                                    f"is not inside the slot mouth {math.degrees(a1):.3f}..{math.degrees(a2):.3f} "
                                    f"deg under the lips")
            m[f"pm{L.inp.k}_opening_deg"] = (math.degrees(lo), math.degrees(hi))
    # ---- air regions: separation, webs ---------------------------------------------------
    regions = _air_regions(rp)
    for name, _ in regions:
        same_layer = p.is_bridgeless and name in ("layer 1 upper slot + layer 1 lower slot",
                                                  "layer 2 upper slot + layer 2 lower slot")
        if " + " in name and not same_layer:
            raise GeometryError(f"air regions merge: {name} form one pocket (steel between them is missing)")
    for i in range(len(regions)):
        for j in range(i + 1, len(regions)):
            d = _dist(regions[i][1], regions[j][1])
            if d < p.min_web - 1e-9:
                raise GeometryError(f"steel between the {regions[i][0]} and the {regions[j][0]} is {d:.3f} mm "
                                    f"(min_web = {p.min_web} mm)")
    rmin = min((float(np.hypot(*P.T).min()) for _, P in regions), default=math.inf)
    if not p.is_bridgeless:
        if rmin - rp.r_shaft < p.min_web - 1e-9:
            raise GeometryError(f"a pocket or hole comes within {rmin - rp.r_shaft:.3f} mm of the bore "
                                f"(min_web = {p.min_web} mm)")
    # ---- q-axis web (half width = distance of any air region from the q-axis line) ---------
    def q_dist(P):
        th = np.arctan2(P[:, 1], P[:, 0])
        r = np.hypot(P[:, 0], P[:, 1])
        return (r * np.sin(hp - np.abs(th))).min()

    named = {n: P for n, P in regions}
    qd = min((q_dist(P) for n, P in regions if "slot" in n or "air region" in n), default=math.inf)
    if p.is_bridgeless:
        qd = min(qd, _bridgeless_q_dist(rp))
    m["q_web"] = 2 * float(qd)
    if m["q_web"] < p.min_q_web - 1e-9:
        raise GeometryError(f"steel between the pockets and the q-axis is {m['q_web']:.3f} mm across the pole "
                            f"boundary, below min_q_web = {p.min_q_web} mm")
    # ---- derived quantities ---------------------------------------------------------------
    for L in rp.layers:
        k = L.inp.k
        m[f"pm{k}_length"] = L.length
        up = next(mg for mg in rp.magnets if mg.layer == k and mg.side > 0)
        m[f"pm{k}_magnet_area"] = up.area()
        slot = named.get(f"layer {k} upper slot")
        if slot is not None and not p.is_bridgeless:
            m[f"pm{k}_bridge_min"] = rp.r_lam - float(np.hypot(*slot.T).max())
    if not p.is_bridgeless:
        for k in (1, 2):
            a, b = named.get(f"layer {k} upper slot"), named.get(f"layer {k} lower slot")
            if a is not None and b is not None:
                m[f"pm{k}_rib_min"] = _dist(a, b)
        a, b = named.get("layer 1 upper slot"), named.get("layer 2 upper slot")
        if a is not None and b is not None:
            m["layer_web"] = _dist(a, b)
        if "cavity" in named and "layer 2 upper slot" in named:
            m["pm2_web_min"] = _dist(named["cavity"], named["layer 2 upper slot"])
    m["magnet_area_total"] = rp.magnet_area()
    m["lamination_area"] = float(sum(f.area() for f in rp.lamination))
    m["lamination_pieces"] = len(rp.lamination)
    m["r_lam"] = rp.r_lam
    return m


def _wall_angle(L, r: float, which: int) -> float:
    """Polar angle where slot wall 1 (pole side, which=+1) or 2 (q side, -1) meets radius r."""
    x0 = L.scalars["x1"] if which > 0 else L.scalars["x2"]
    th = L.inp.theta
    b = x0 * math.cos(th)
    disc = b * b - (x0 * x0 - r * r)
    if disc < 0:
        raise GeometryError(f"layer {L.inp.k}: slot wall does not reach radius {r:.3f} mm")
    s = -b + math.sqrt(disc)
    return math.atan2(s * math.sin(th), x0 + s * math.cos(th))


def _bridgeless_q_dist(rp: RotorPole) -> float:
    """Bridgeless: distance from the q-side wall of each upper slot to the q-axis line. The wall
    is straight from the strip's inner corner to the lips (which are steel up to the OD
    openings), so the minimum is at one of its ends."""
    hp = rp.half_pitch
    qn = np.array([-math.sin(hp), math.cos(hp)])  # unit normal of the q-axis line (towards +theta)
    best = math.inf
    for L in rp.layers:
        k = L.inp.k
        u = np.array(L.u)
        n = np.array([-u[1], u[0]])
        a = np.array([L.inp.xpm, 0.0]) + getattr(rp.p, f"pm{k}_slot_start") * u + L.inp.W / 2 * n
        r_lip = rp.r_lam - rp.p.lip_depth
        a_od = _wall_angle(L, r_lip, -1)
        b = r_lip * np.array([math.cos(a_od), math.sin(a_od)])
        best = min(best, float(-(a @ qn)), float(-(b @ qn)))
    return best
