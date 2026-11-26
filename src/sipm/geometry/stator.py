"""Stator pole pitch and winding layout (the stator of the Maxwell models).

Slot (slot frame: slot centre line on +x): a parallel opening of width ``slot_opening``
from the bore to x2 = Rsi cos(asin(Wo / 2 Rsi)) + ``slot_tip_depth``, a taper at
``slot_tip_angle`` to the slot width at x3, parallel sides to the flat bottom at
Rsi + ``slot_depth``; every slot corner, including the two bore corners, is filleted with
``slot_fillet``. Conductors: ``cond_count`` touching rectangles of ``cond_height`` (radial) x
``cond_width`` from ``cond_start_radius``, merged into one coil region per slot that carries
``hairpin_layers`` turns (model: 8 rectangles, Nac = 7).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..params import Params
from .primitives import (
    Arc,
    Chain,
    GeometryError,
    Line,
    Loop,
    chain_area,
    chain_to_loop,
    chain_transformed,
    fillet_chain,
    polygon_chain,
    rotation,
)


@dataclass
class Conductor:
    loop: Loop
    slot: int            # 0-based slot index within the pole pitch
    layer: int           # 1 (one coil region per slot)
    centre_angle: float  # rad, pole frame
    turns: int = 1
    chain: Chain = ()


@dataclass
class StatorPole:
    p: Params
    loop: Loop
    conductors: list
    slot_angles: np.ndarray
    metrics: dict
    chain: Chain = ()


def _pol(r: float, th: float) -> tuple[float, float]:
    return (r * math.cos(th), r * math.sin(th))


def slot_profile(p: Params) -> list[tuple[float, float]]:
    """Unfilleted slot corners in the slot frame, from the upper bore corner round to the lower."""
    ri = p.r_stator_in
    wo, ws = p.slot_opening, p.slot_width
    xb0 = ri * math.cos(math.asin(wo / 2 / ri))
    x2 = xb0 + p.slot_tip_depth
    x3 = x2 + (ws / 2 - wo / 2) / math.tan(math.radians(p.slot_tip_angle))
    xb = ri + p.slot_depth
    if x3 >= xb:
        raise GeometryError("slot taper ends beyond the slot bottom (slot_depth too small)")
    return [(xb0, wo / 2), (x2, wo / 2), (x3, ws / 2), (xb, ws / 2), (xb, -ws / 2), (x3, -ws / 2),
            (x2, -wo / 2), (xb0, -wo / 2)]


def build_stator(p: Params) -> StatorPole:
    p.validate()
    hp = math.pi / p.poles
    nsl = p.slots // p.poles
    tau = 2 * math.pi / p.slots
    ri, ro = p.r_stator_in, p.stator_od / 2
    prof = slot_profile(p)
    xb = prof[3][0]
    yoke = ro - math.hypot(xb, p.slot_width / 2)
    if yoke < 2.0:
        raise GeometryError(f"stator yoke is only {yoke:.3f} mm (slot too deep for stator_od)")
    x3 = prof[2][0]
    tooth_min = 2 * x3 * math.sin(tau / 2) - p.slot_width
    tooth_tip = 2 * ri * math.sin(tau / 2) - p.slot_opening
    if tooth_min < 1.0:
        raise GeometryError(f"tooth width at the slot top is {tooth_min:.3f} mm (< 1 mm)")
    c0, c1 = p.cond_start_radius, p.cond_start_radius + p.cond_count * p.cond_height
    if p.cond_width > p.slot_width + 1e-12:
        raise GeometryError("conductor wider than the slot")
    if c0 < x3 - 1e-12 or c1 > xb + 1e-12:
        raise GeometryError(f"conductors ({c0:.3f}..{c1:.3f} mm) do not fit between the slot taper ({x3:.3f}) "
                            f"and the slot bottom ({xb:.3f} mm)")
    slot_angles = np.array([-hp + (j + 0.5) * tau for j in range(nsl)])
    o = (0.0, 0.0)
    a_open = math.asin(p.slot_opening / 2 / ri)
    segs: list = [Line(_pol(ri, -hp), _pol(ro, -hp)), Arc(o, ro, -hp, 2 * hp), Line(_pol(ro, hp), _pol(ri, hp))]
    corners: list[tuple[float, float]] = []
    cur = hp
    for th in slot_angles[::-1]:
        R = rotation(th)
        pts = [tuple(R @ np.array(q)) for q in prof]
        a_hi = th + a_open
        segs.append(Arc(o, ri, cur, a_hi - cur))
        for q0, q1 in zip(pts[:-1], pts[1:]):
            segs.append(Line(q0, q1))
        corners += pts
        cur = th - a_open
    segs.append(Arc(o, ri, cur, -hp - cur))
    chain: Chain = tuple(segs)
    if p.slot_fillet > 0:
        for q in corners:
            idx = next(i for i, s in enumerate(chain) if math.dist(s.start(), q) < 1e-9)
            chain = fillet_chain(chain, idx, p.slot_fillet)
    conductors = []
    hw = p.cond_width / 2
    rect = polygon_chain([(c0, -hw), (c1, -hw), (c1, hw), (c0, hw)])
    for j, th in enumerate(slot_angles):
        ch = chain_transformed(rect, rotation(th))
        conductors.append(Conductor(chain_to_loop(ch), j, 1, float(th), p.hairpin_layers, ch))
    loop = chain_to_loop(chain)
    metrics = dict(slot_bottom_radius=xb, yoke=yoke, tooth_min=tooth_min, tooth_tip=tooth_tip,
                   stator_area=chain_area(chain), cond_area=chain_area(rect))
    return StatorPole(p=p, loop=loop, conductors=conductors, slot_angles=slot_angles, metrics=metrics, chain=chain)


def winding_layout(p: Params, slot_angles: np.ndarray):
    """Full-pitch integer-slot layout with 60-degree phase belts (illustrative).

    Returns a list of (phase, sign) per slot of the modelled pole pitch and the
    electrical axis angle (deg) of each phase in the same frame as ``slot_angles``.
    The axis is where the phase's own MMF drives flux radially outwards:
    phi = arg(sum s*exp(j*p*theta)) - 90 deg.
    """
    pp = p.pole_pairs
    hp = math.pi / p.poles
    belts = [("A", +1), ("C", -1), ("B", +1), ("A", -1), ("C", +1), ("B", -1)]
    out = []
    for th in slot_angles:
        e = math.degrees(pp * (th + hp)) % 360.0  # 0 at the start of the pole pitch
        out.append(belts[int(e // 60.0) % 6])
    axes = {}
    for ph in "ABC":
        z = sum(s * np.exp(1j * pp * th) for (q, s), th in zip(out, slot_angles) if q == ph)
        axes[ph] = (math.degrees(np.angle(z)) - 90.0)
    return out, axes
