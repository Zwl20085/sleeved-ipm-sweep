"""Rotor pole of the two sleeved topologies, built the way the Maxwell models build it.

Pole frame: shaft axis at the origin, d-axis of the modelled pole on +x, the pole spans
[-pi/p, +pi/p]. The lamination is a constructive expression (see :mod:`.csg`) evaluated
exactly; the operation order is the order of the Maxwell model trees:

sleeved, bridged (``01_ICEM_01_simplerotor``)
    sector - q-axis lightening holes - layer-1 strips + centre-rib band - layer-2 strips
    (each minus its tilted end band) - central cavity + magnet lugs - corner reliefs
    - d-axis hole + outer bridge bands
sleeved, bridgeless (``02_ICEM_01_teslamotor1``)
    as above without cavity and bridge bands; + pole-tip lips with angular openings
    - rib cut (the region between the upper and lower slot of a layer, which separates the
    lamination into pieces)

Where the Maxwell model used a hand-placed helper whose only job was to cover something
(fixed angular spans of the bridge bands and lips, fixed extents of the rib band and tilted
end band, hard-coded prism coordinates that cut the bridgeless ribs) the helper is derived
from the geometry instead; see docs/methods.md for the list and the rib-cut rule.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..params import Params
from . import csg
from .construction import LayerInput, LayerPoints, layer_points
from .primitives import (
    MIRROR_X,
    Chain,
    GeometryError,
    Hermite,
    Line,
    Loop,
    Pt,
    annular_sector,
    chain_area,
    chain_ccw,
    chain_mirror_x,
    chain_reversed,
    chain_to_loop,
    chain_transformed,
    circle_chain,
    fillet_chain,
    offset_line_segment,
    polygon_chain,
    rotation,
)

FAR = 3.0          # slot strips run this many lamination radii beyond their inner end (through the OD)
BAND_LONG = 40.0   # mm, half length of the tilted end band of the layer-2 slot (covers the slot)
OPEN_PAD = 1.0     # mm, radial overshoot of the opening cutters beyond the lip

# Magnet retention features of the sleeved models (fixed polygons, mm). Each polygon is given
# in a frame whose origin is a magnet corner of the nominal (unshrunk) magnet and whose x axis
# points to a second corner; y is x turned +90 deg. Corner names: P8/P9 inner end (pole side /
# q side), P11/P10 outer end. "W/2" stands for half the slot width (the model has the number
# 2.15295.. / 1.82450.., i.e. half of its RD_Wpm). Lugs are added to the lamination, reliefs
# removed from it.
LUGS: dict[str, tuple[int, str, str, tuple]] = {
    "layer2_inner_lug": (2, "P9", "P10", ((0.0, -0.6), (0.0, "-W/2"), (-2.0, 0.8), (-0.8, 0.4))),
    "layer2_outer_lug": (2, "P10", "P11", ((0.4, 0.0), (1.6, 0.0), (-0.6, 2.2), (-0.6, 0.8))),
    "layer1_inner_lug": (1, "P9", "P8", ((0.6, 0.0), ("W/2", 0.0), (-0.4, -1.8), (-0.6, -0.8))),
    "layer1_outer_lug": (1, "P10", "P11", ((0.5, 0.0), (1.5, 0.0), (-1.0, 2.0), (-1.0, 1.0))),
}
RELIEFS: dict[str, tuple[int, str, str, tuple]] = {
    "layer2_inner_relief": (2, "P8", "P11", ((1.0, 0.0), (0.0, -0.2), (-0.4, 0.0))),
    "layer2_outer_relief": (2, "P11", "P10", ((0.0, 0.8), (-0.2, 0.0), (0.0, -0.2))),
    "layer1_inner_relief": (1, "P8", "P11", ((0.6, 0.0), (0.0, -0.1), (-0.1, 0.0))),
    "layer1_outer_relief": (1, "P11", "P8", ((-0.2, 0.0), (0.0, 0.2), (0.8, 0.0))),
}
# Central cavity (layer-2 construction points FB, FBd; see _cavity_chain): fixed numbers of the
# model's FB_5..FB_8 curves, mm.
CAVITY = dict(tip=1.0, step=2.0, tangent_tip=3.0, tangent_tip_steep=6.0, tangent_fb=2.0, tangent_fbd=6.0)


# ======================================================================================
# data classes
# ======================================================================================
@dataclass
class Magnet:
    chain: Chain                 # exact boundary (lines + corner arcs), counter-clockwise
    magnetisation: np.ndarray    # unit vector, pole frame
    layer: int
    side: int                    # +1 upper (y > 0), -1 lower
    axis_angle: float            # rad, angle of the long axis to +x
    length: float                # derived magnet length (nominal, before clearances)
    thickness: float             # actual thickness (slot width - side clearance)
    loop: Loop = field(init=False)

    def __post_init__(self):
        self.loop = chain_to_loop(self.chain)

    @property
    def centre(self) -> np.ndarray:
        P = self.loop.discretize(2.0)
        return P.mean(axis=0)

    def area(self) -> float:
        return chain_area(self.chain)


@dataclass
class LaminationFace:
    """One lamination piece: exact chains and the solver loops (Hermite curves as biarcs)."""

    outer: Chain
    holes: tuple

    def area(self) -> float:
        return chain_area(self.outer) + sum(chain_area(h) for h in self.holes)

    def solver_loops(self) -> tuple[Loop, list[Loop]]:
        return chain_to_loop(self.outer), [chain_to_loop(h).ccw() for h in self.holes]


@dataclass
class RotorPole:
    p: Params
    r_lam: float
    r_shaft: float
    r_sleeve_out: float
    half_pitch: float                    # rad
    lamination: list                     # list[LaminationFace], full pole
    magnets: list                        # list[Magnet], full pole
    sleeve: Loop
    layers: list                         # list[LayerPoints]
    features: dict                       # feature name -> representative point (checks, figures)
    metrics: dict
    _expr_builder: object = None
    _half: list | None = None

    @property
    def sector(self) -> Loop:
        """Full-pole annular sector (shaft to lamination OD) before any hole."""
        return chain_to_loop(annular_sector(self.r_shaft, self.r_lam, -self.half_pitch, self.half_pitch))

    def half_lamination(self) -> list[LaminationFace]:
        """Lamination in the upper half pole 0 <= theta <= pi/p (mechanical model), with the
        Hermite curves already replaced by biarcs."""
        if self._half is None:
            self._half = _evaluate(self._expr_builder, 0.0, self.half_pitch, solver=True)
        return self._half

    def half_magnets(self) -> list[Magnet]:
        return [m for m in self.magnets if m.side == +1]

    def magnet_area(self) -> float:
        return sum(m.area() for m in self.magnets)


# ======================================================================================
# helpers
# ======================================================================================
def _pol(r: float, th: float) -> Pt:
    return (r * math.cos(th), r * math.sin(th))


def layer_inputs(p: Params) -> tuple[LayerInput, LayerInput]:
    def one(k: int, fb: float, tilt: float) -> LayerInput:
        return LayerInput(k=k, xpm=getattr(p, f"pm{k}_d_crossing"),
                          theta=math.radians(getattr(p, f"pm{k}_angle")) / 2,
                          W=getattr(p, f"pm{k}_width"), dpm=getattr(p, f"pm{k}_end_gap"),
                          K=getattr(p, f"pm{k}_length_ratio"), FB=fb,
                          INNO=getattr(p, f"pm{k}_inner_offset"), tINNO=tilt,
                          WO=getattr(p, f"pm{k}_bridge_thickness"))
    # layer 1 ends on the flat rib line (no tilt); layer 2 ends on the tilted web
    return one(1, p.pm1_rib_half_width, 0.0), one(2, p.pm2_web_width, math.radians(p.pm2_inner_tilt))


def construction_points(p: Params) -> list[LayerPoints]:
    i1, i2 = layer_inputs(p)
    return [layer_points(i1, p.r_lam), layer_points(i2, p.r_lam)]


def _frame_poly(L: LayerPoints, origin: str, xto: str, pts: tuple) -> Chain:
    o = np.array(L.p(origin))
    xd = np.array(L.p(xto)) - o
    xd /= np.linalg.norm(xd)
    yd = np.array([-xd[1], xd[0]])
    half = L.inp.W / 2

    def val(v):
        if isinstance(v, str):
            return -half if v.startswith("-") else half
        return v

    out = [tuple(o + val(a) * xd + val(b) * yd) for a, b in pts]
    return polygon_chain(out)


def _strip(L: LayerPoints, start: float) -> Chain:
    """Upper slot strip of width W along the magnet axis from ``start`` (beyond the d-axis
    crossing) through the OD."""
    u = np.array(L.u)
    n = np.array([-u[1], u[0]])
    o = np.array([L.inp.xpm, 0.0])
    far = FAR * L.r_lam
    w = L.inp.W / 2
    return polygon_chain([o + start * u - w * n, o + far * u - w * n, o + far * u + w * n, o + start * u + w * n])


def _rib_band(half: float, r_lam: float) -> Chain:
    return polygon_chain([(0.0, -half), (r_lam + 1.0, -half), (r_lam + 1.0, half), (0.0, half)])


def _end_band(L: LayerPoints) -> Chain:
    """Band under the tilted inner end of the layer-2 slot (Maxwell Rectangle12): height
    FB + INNO, rotated by -tINNO about the d-axis crossing and moved by d2 along x."""
    h = L.inp.FB + L.inp.INNO
    R = rotation(-L.inp.tINNO)
    base = np.array([L.inp.xpm + L.scalars["d2"], 0.0])
    pts = [base + R @ np.array(q) for q in ((-BAND_LONG, 0.0), (BAND_LONG, 0.0), (BAND_LONG, h), (-BAND_LONG, h))]
    return polygon_chain(pts)


def _q_hole(p: Params, at: float) -> Chain:
    """Lightening hole centred on the q-axis at angle ``at``: inner arc, two lines to the apex,
    three corner fillets."""
    a = math.radians(p.q_hole_half_angle)
    ri, ra = p.q_hole_inner_radius, p.q_hole_apex_radius
    if ra <= ri:
        raise GeometryError("q_hole_apex_radius must exceed q_hole_inner_radius")
    pa, pb, apex = _pol(ri, -a), _pol(ri, a), (ra, 0.0)
    from .primitives import Arc
    ch: Chain = (Line(pa, apex), Line(apex, pb), Arc((0.0, 0.0), ri, a, -2 * a))
    if p.q_hole_fillet > 0:
        for idx in (2, 1, 0):
            ch = fillet_chain(ch, idx, p.q_hole_fillet)
    return chain_transformed(ch, rotation(at))


def _d_hole(p: Params, r_lam: float) -> Chain:
    return circle_chain((r_lam - p.d_hole_depth, 0.0), p.d_hole_diameter / 2)


def _cavity_chain(L2: LayerPoints) -> Chain:
    """Central cavity on the d-axis inside the layer-2 V (model FB_5..FB_8 and mirror)."""
    ti = L2.inp.tINNO
    xFB, yFB = L2.p("FB")
    xFBd, yFBd = L2.p("FBd")
    Ya = yFBd - (L2.scalars["dFB"] - L2.inp.INNO) / math.cos(ti) - CAVITY["step"] * math.tan(ti)
    t_tip = CAVITY["tangent_tip_steep"] if ti > math.radians(30) else CAVITY["tangent_tip"]
    fb5 = Hermite((xFB + CAVITY["tip"], 0.0), (0.0, t_tip), (xFB, yFB),
                  (-CAVITY["tangent_fb"] * math.cos(ti), CAVITY["tangent_fb"] * math.sin(ti)))
    q1 = (xFBd + CAVITY["step"], Ya)
    fb7 = Hermite((xFBd, Ya * 3 / 6), (0.0, 0.5 * L2.inp.INNO), q1,
                  (CAVITY["tangent_fbd"] * math.cos(ti), -CAVITY["tangent_fbd"] * math.sin(ti)))
    upper = (fb5, Line((xFB, yFB), q1), fb7.reversed())       # (xFB+1,0) -> ... -> (xFBd, Ya/2)
    lower = chain_reversed(chain_mirror_x_open(upper))           # (xFBd, -Ya/2) -> ... -> (xFB+1, 0)
    ch = upper + (Line((xFBd, Ya * 3 / 6), (xFBd, -Ya * 3 / 6)),) + lower
    return chain_ccw(ch)


def chain_mirror_x_open(segs) -> tuple:
    """Mirror an open chain about the x axis without reversing it."""
    return tuple(s.transformed(MIRROR_X) for s in segs)


def _magnet_chain(p: Params, L: LayerPoints) -> Chain:
    P = L.P
    ch = polygon_chain([P["P9"], P["P8"], P["P11"], P["P10"]])  # P10->P9 = q-side face, P9->P8 = inner end
    if chain_area(ch) <= 0:
        raise GeometryError(f"layer {L.inp.k}: magnet length is not positive")
    if p.pm_side_clearance > 0:
        ch = offset_line_segment(ch, 3, p.pm_side_clearance)
    if p.pm_end_clearance > 0:
        ch = offset_line_segment(ch, 0, p.pm_end_clearance)
    if p.pm_corner_radius > 0:
        for idx in (3, 2, 1, 0):
            ch = fillet_chain(ch, idx, p.pm_corner_radius)
    return ch


def opening_angles(p: Params, L: LayerPoints) -> tuple[float, float]:
    k = L.inp.k
    hi = L.scalars["opening_reference"] - math.radians(getattr(p, f"pm{k}_opening_shift"))
    return hi - math.radians(getattr(p, f"pm{k}_opening_angle")), hi


def _rib_cut(p: Params, L: LayerPoints) -> Chain:
    """Bridgeless rotor: air region joining the upper and lower slot of a layer (derived
    replacement of the model's hand-placed cutting prisms). Convex hull of the slot strip's
    inner-end corners, the slot's inner pole-side corner P5 and those inner magnet corners
    (P8, P9) that lie before the strip's inner end, together with their mirror images."""
    k = L.inp.k
    u = np.array(L.u)
    n = np.array([-u[1], u[0]])
    o = np.array([L.inp.xpm, 0.0])
    s0 = getattr(p, f"pm{k}_slot_start")
    w = L.inp.W / 2
    pts = [o + s0 * u + w * n, o + s0 * u - w * n, np.array(L.p("P5"))]
    for c in ("P8", "P9"):
        q = np.array(L.p(c))
        if (q - o) @ u < s0:
            pts.append(q)
    pts += [np.array([q[0], -q[1]]) for q in pts]
    hull = _convex_hull([tuple(q) for q in pts])
    return polygon_chain(hull)


def _convex_hull(points: list[Pt]) -> list[Pt]:
    pts = sorted(set(points))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[Pt] = []
    for q in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], q) <= 1e-12:
            lower.pop()
        lower.append(q)
    upper: list[Pt] = []
    for q in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], q) <= 1e-12:
            upper.pop()
        upper.append(q)
    return lower[:-1] + upper[:-1]


# ======================================================================================
# expression builders
# ======================================================================================
@dataclass
class _Build:
    """Everything needed to evaluate the lamination for any angular wedge."""

    p: Params
    layers: list
    r_shaft: float
    hp: float
    curved_holes: list   # exact chains with Hermite curves, subtracted last (they touch nothing)
    features: dict
    split_angle: float | None = None


def _sleeved_expr(b: _Build, base: csg.Expr) -> csg.Expr:
    p, (L1, L2), hp = b.p, b.layers, b.hp
    Leaf = csg.Leaf
    lam = base
    for at in (-hp, hp):
        lam = lam - Leaf(_q_hole(p, at), "q_hole")
    s1 = _strip(L1, p.pm1_slot_start)
    lam = lam - Leaf(s1, "slot1") - Leaf(chain_mirror_x(s1), "slot1_lower")
    lam = lam | (Leaf(_rib_band(p.pm1_rib_half_width, b.p.r_lam), "rib1") & base)
    s2, e2 = _strip(L2, p.pm2_slot_start), _end_band(L2)
    lam = lam - (Leaf(s2, "slot2") - Leaf(e2, "end2"))
    lam = lam - (Leaf(chain_mirror_x(s2), "slot2_lower") - Leaf(chain_mirror_x(e2), "end2_lower"))
    for name, (k, o, x, pts) in LUGS.items():
        ch = _frame_poly(b.layers[k - 1], o, x, pts)
        lam = lam | ((Leaf(ch, name) | Leaf(chain_mirror_x(ch), name + "_lower")) & base)
    for name, (k, o, x, pts) in RELIEFS.items():
        ch = _frame_poly(b.layers[k - 1], o, x, pts)
        lam = lam - Leaf(ch, name) - Leaf(chain_mirror_x(ch), name + "_lower")
    lam = lam - Leaf(_d_hole(p, p.r_lam), "d_hole")
    R = p.r_lam
    if not p.is_bridgeless:
        ts = b.split_angle
        band1 = annular_sector(R - p.pm1_bridge_thickness, R, -ts, ts)
        bands = Leaf(band1, "bridge1")
        for a0, a1 in ((ts, hp), (-hp, -ts)):
            bands = bands | Leaf(annular_sector(R - p.pm2_bridge_thickness, R, a0, a1), "bridge2")
        lam = lam | (bands & base)
    else:
        lips: csg.Expr = Leaf(annular_sector(R - p.lip_depth, R, -hp, hp), "lips")
        for L in (L1, L2):
            lo, hi = opening_angles(p, L)
            for sgn in (1, -1):
                a0, a1 = (lo, hi) if sgn > 0 else (-hi, -lo)
                lips = lips - Leaf(annular_sector(R - p.lip_depth - OPEN_PAD, R + OPEN_PAD, a0, a1), "opening")
        lam = lam | (lips & base)
        for L in (L1, L2):
            lam = lam - Leaf(_rib_cut(p, L), f"rib_cut{L.inp.k}")
    return lam


def _evaluate(b: _Build, a0: float, a1: float, solver: bool) -> list[LaminationFace]:
    base = csg.Leaf(annular_sector(b.r_shaft, b.p.r_lam, a0, a1), "sector")
    expr = _sleeved_expr(b, base)
    if solver:
        for ch in b.curved_holes:
            expr = expr - csg.Leaf(_loop_chain(chain_to_loop(ch)), "curved")
        faces = csg.evaluate(expr)
        return [LaminationFace(f.outer, f.holes) for f in faces]
    faces = csg.evaluate(expr)
    out = []
    pending = list(b.curved_holes)
    for f in faces:
        holes = list(f.holes)
        for ch in list(pending):
            q = ch[0].start()
            if csg._crossings(f.outer, q) % 2 == 1:
                holes.append(chain_reversed(ch))
                pending.remove(ch)
        out.append(LaminationFace(f.outer, tuple(holes)))
    if pending:
        raise GeometryError("a pocket or the central cavity lies outside the lamination")
    return out


def _loop_chain(loop: Loop) -> Chain:
    """Exact chain of a solver Loop (lines and minor arcs)."""
    from .primitives import Arc, arc_from
    out = []
    n = len(loop)
    for i in range(n):
        p0, p1, c = tuple(loop.pts[i]), tuple(loop.pts[(i + 1) % n]), loop.arcs[i]
        if c is None:
            out.append(Line(p0, p1))
        else:
            _, d = Loop._sweep(np.array(p0), np.array(p1), c)
            a = arc_from(tuple(c), p0, p1, ccw=d > 0)
            out.append(Arc(a.c, a.r, a.a0, a.da))
    return tuple(out)


# ======================================================================================
# public
# ======================================================================================
def build_rotor(p: Params) -> RotorPole:
    """Build one rotor pole (full pole, pole frame) and check its feasibility."""
    from .checks import check_rotor

    p.validate()
    hp = math.pi / p.poles
    r_sh = p.shaft_diameter / 2
    layers = construction_points(p)
    L1, L2 = layers
    if L1.inp.xpm <= L2.inp.xpm:
        raise GeometryError("layer order: layer 1 must cross the d-axis outside layer 2 "
                            "(pm1_d_crossing > pm2_d_crossing)")
    # the lightening holes' inner arc (centred on the shaft axis) must stay min_web outside the
    # bore; checked here because a bore that touches or crosses that arc leaves no valid outline
    web = p.q_hole_inner_radius - r_sh
    if web < p.min_web - 1e-9:
        raise GeometryError(f"the q-axis lightening holes come within {web:.3f} mm of the bore "
                            f"(q_hole_inner_radius {p.q_hole_inner_radius} mm - shaft radius {r_sh} mm; "
                            f"min_web = {p.min_web} mm)")
    for L in layers:
        lmin = 2 * p.pm_corner_radius + 0.5
        if L.scalars["free_length"] <= 0:
            raise GeometryError(f"layer {L.inp.k}: the slot does not reach the bridge circle beyond the magnet's "
                                f"inner end (free slot length {L.scalars['free_length']:.3f} mm <= 0)")
        if L.length < lmin:
            raise GeometryError(f"layer {L.inp.k}: derived magnet length {L.length:.3f} mm is shorter than "
                                f"{lmin:.3f} mm (2 x pm_corner_radius + 0.5 mm)")
    curved: list = []
    features: dict = {}
    split = None
    if not p.is_bridgeless:
        a41 = math.atan2(L1.p("P4")[1], L1.p("P4")[0])
        a32 = math.atan2(L2.p("P3")[1], L2.p("P3")[0])
        a42 = math.atan2(L2.p("P4")[1], L2.p("P4")[0])
        if a41 >= a32:
            raise GeometryError("bridge bands: the layer-1 and layer-2 slot mouths overlap at the bridge "
                                f"circles ({math.degrees(a41):.3f} >= {math.degrees(a32):.3f} deg)")
        if a42 >= hp:
            raise GeometryError("bridge band does not cover the slot: the layer-2 slot meets its bridge "
                                f"circle beyond the q-axis ({math.degrees(a42):.3f} deg)")
        split = 0.5 * (a41 + a32)
        curved.append(_cavity_chain(L2))
        features["cavity"] = (0.5 * (L2.p("FB")[0] + L2.p("FBd")[0]), 0.0)
    b = _Build(p=p, layers=layers, r_shaft=r_sh, hp=hp, curved_holes=curved, features=features, split_angle=split)
    faces = _evaluate(b, -hp, hp, solver=False)
    magnets = []
    for L in layers:
        ch = _magnet_chain(p, L)
        th = L.inp.theta
        nvec = np.array([math.sin(th), -math.cos(th)])
        thick = L.inp.W - p.pm_side_clearance
        magnets.append(Magnet(ch, nvec, L.inp.k, +1, th, L.length, thick))
        magnets.append(Magnet(chain_mirror_x(ch), np.array([nvec[0], -nvec[1]]), L.inp.k, -1, -th, L.length, thick))
    sleeve = chain_to_loop(annular_sector(p.r_lam, p.r_sleeve_out, -hp, hp))
    rp = RotorPole(p=p, r_lam=p.r_lam, r_shaft=r_sh, r_sleeve_out=p.r_sleeve_out, half_pitch=hp,
                   lamination=faces, magnets=magnets, sleeve=sleeve, layers=layers, features=features,
                   metrics={}, _expr_builder=b)
    rp.metrics = check_rotor(rp)
    return rp
