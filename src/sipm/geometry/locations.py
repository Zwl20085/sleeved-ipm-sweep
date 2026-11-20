"""Named load-carrying sections of the lamination, where the headline stress is read.

The global lamination stress maximum of these rotors sits at sharp corners of the Maxwell
outline (slot ends, lug corners, rib-cut corners). In a linear-elastic model such corners are
singular: the peak grows without limit as the mesh is refined. The design quantity is
therefore the peak von Mises stress over a few named sections that carry the centrifugal load,
each taken over the section's span minus ``stress_corner_margin`` at both ends, so that the
singular end corners are left out and the value converges with the mesh (docs/methods.md,
section 3).

Sleeved bridge-retained rotor (upper half pole):

* ``bridge1``, ``bridge2`` - the outer bridge of each layer: the steel between the slot's
  closing arc (radius r_lam - pmk_bridge_thickness) and the OD, over the slot mouth;
* ``rib1`` - the layer-1 centre rib: the steel between the flat inner end of the layer-1 slot
  (y = pm1_rib_half_width + pm1_inner_offset) and the d-axis;
* ``web2`` - the layer-2 web: the band of thickness pm2_web_width between the tilted inner end
  of the layer-2 slot and the central cavity;
* ``qhole`` - the steel along the three corner fillets of the q-axis lightening hole (a band
  ``QHOLE_BAND`` deep), where the hoop load of the q-axis body concentrates. The fillets are
  smooth, so this section needs no corner margin.

Sleeved bridgeless rotor: no steel joins the three lamination pieces; the sleeve and the
(bonded) magnets hold them together. Its load-carrying sections are

* ``qhole`` - as above, in the q-axis body, which carries the hoop load of the pole;
* ``lip1_pole``, ``lip1_q``, ``lip2_pole``, ``lip2_q`` - the pole-tip lips that cover each slot
  mouth under the OD (layer, pole side or q side of the opening), cantilevered from the slot
  walls to the OD openings and loaded by their own centrifugal force and the sleeve contact
  pressure.

The other stress peaks of the bridgeless rotor sit at sharp corners (rib cut, slot ends, lug
corners next to the bonded magnets) and are left to the mesh-dependent global maximum.

Each section is an exact region in the pole frame (an annular sector or a band along a
straight edge), found from the lamination's own boundary segments, so it follows the geometry
for any parameter set. Lengths in mm, angles in rad.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .primitives import GeometryError

TOL = 1e-6        # mm, node-on-boundary tolerance and geometric matching tolerance
QHOLE_BAND = 1.0  # mm, depth of the steel band along the q-hole fillets


@dataclass(frozen=True)
class Section:
    """A named region. ``kind`` 'sector': r0 <= r <= r1 and a0 <= theta <= a1 about the centre
    ``c`` (shaft axis by default; theta counter-clockwise from a0). ``kind`` 'band':
    o + s d + t n with s0 <= s <= s1 and 0 <= t <= depth (d unit vector along the section, n
    unit normal into the steel). Several regions may share a name (one section)."""

    name: str
    kind: str
    r0: float = 0.0
    r1: float = 0.0
    a0: float = 0.0
    a1: float = 0.0
    o: tuple = (0.0, 0.0)
    d: tuple = (1.0, 0.0)
    n: tuple = (0.0, 1.0)
    s0: float = 0.0
    s1: float = 0.0
    depth: float = 0.0
    c: tuple = (0.0, 0.0)

    def contains(self, xy: np.ndarray, tol: float = TOL) -> np.ndarray:
        """Boolean mask of the points (N x 2) inside the section (boundary included)."""
        xy = np.atleast_2d(np.asarray(xy, float))
        if self.kind == "sector":
            q = xy - np.asarray(self.c)
            r = np.hypot(q[:, 0], q[:, 1])
            ta = tol / max(self.r0, 1.0)
            da = np.mod(np.arctan2(q[:, 1], q[:, 0]) - (self.a0 - ta), 2 * math.pi)
            return (r >= self.r0 - tol) & (r <= self.r1 + tol) & (da <= self.a1 - self.a0 + 2 * ta)
        q = xy - np.asarray(self.o)
        s = q @ np.asarray(self.d)
        t = q @ np.asarray(self.n)
        return (s >= self.s0 - tol) & (s <= self.s1 + tol) & (t >= -tol) & (t <= self.depth + tol)

    @property
    def span(self) -> float:
        """Length of the section along its edge (mm), after the corner margins."""
        if self.kind == "sector":
            return self.r0 * (self.a1 - self.a0)
        return self.s1 - self.s0


def _segments(rp):
    """Every boundary segment of the lamination in the upper half pole."""
    for face in rp.lamination:
        for ch in (face.outer, *face.holes):
            for s in ch:
                yield s


def _qhole_sections(rp) -> list[Section]:
    """Bands along the corner fillets of the q-axis lightening hole (upper half pole)."""
    rf = rp.p.q_hole_fillet
    if rf <= 0:
        return []
    out = []
    for s in _segments(rp):
        if s.kind != "arc" or abs(s.r - rf) > TOL or math.hypot(*s.c) < rp.p.q_hole_inner_radius - rf - TOL:
            continue
        if math.atan2(s.c[1], s.c[0]) <= 0:
            continue
        a0, a1 = sorted((s.a0, s.a0 + s.da))
        out.append(Section("qhole", "sector", r0=rf, r1=rf + QHOLE_BAND, a0=a0, a1=a1, c=tuple(map(float, s.c))))
    if not out:
        raise GeometryError("stress section qhole: the fillets of the q-axis lightening hole were not found")
    return out


def _arc_spans(rp, radius: float) -> list[tuple[float, float]]:
    """Angular spans (a0 < a1, upper half) of boundary arcs centred on the axis with this radius."""
    out = []
    for s in _segments(rp):
        if s.kind != "arc" or abs(s.r - radius) > TOL or math.hypot(*s.c) > TOL:
            continue
        a0, a1 = sorted((s.a0, s.a0 + s.da))
        if a1 > 0:
            out.append((max(a0, 0.0), a1))
    return out


def _span_containing(spans, lo: float, hi: float, name: str) -> tuple[float, float]:
    """The single span that lies within [lo, hi] (the slot mouth of one layer)."""
    hits = [sp for sp in spans if sp[0] >= lo - 1e-9 and sp[1] <= hi + 1e-9]
    if len(hits) != 1:
        raise GeometryError(f"stress section {name}: expected one boundary arc in "
                            f"{math.degrees(lo):.3f}..{math.degrees(hi):.3f} deg, found {len(hits)}")
    return hits[0]


def _sector(name: str, r_in: float, r_out: float, a0: float, a1: float, margin: float) -> Section:
    m = margin / r_in
    if a1 - a0 <= 2 * m:
        raise GeometryError(f"stress section {name}: span {r_in * (a1 - a0):.3f} mm is not longer than twice "
                            f"stress_corner_margin ({margin} mm)")
    return Section(name, "sector", r0=r_in, r1=r_out, a0=a0 + m, a1=a1 - m)


def _end_edge(rp, d: np.ndarray, a: np.ndarray, b: np.ndarray, name: str) -> tuple[np.ndarray, float]:
    """The slot's straight inner-end edge: the longest boundary line segment parallel to ``d``
    whose midpoint lies within 0.1 mm of the line through the construction points a, b (the
    model's end edge is not exactly through them) and between them along d. Returns (start
    point, length) with the edge running along +d from the start point."""
    n = np.array([-d[1], d[0]])
    lo, hi = sorted(((a - a) @ d, (b - a) @ d))
    best = None
    for seg in _segments(rp):
        if seg.kind != "line":
            continue
        p0, p1 = np.array(seg.start()), np.array(seg.end())
        v = p1 - p0
        if abs(v @ n) > 1e-9 * max(1.0, np.linalg.norm(v)):
            continue
        mid = 0.5 * (p0 + p1)
        if abs((mid - a) @ n) > 0.1 or not lo - 1.0 <= (mid - a) @ d <= hi + 1.0:
            continue
        length = float(np.linalg.norm(v))
        if best is None or length > best[1]:
            best = (p0 if v @ d > 0 else p1, length)
    if best is None:
        raise GeometryError(f"stress section {name}: the slot's straight inner end was not found")
    return best


def _band(name: str, o, d, n, s0: float, s1: float, depth: float, margin: float) -> Section:
    if s1 - s0 <= 2 * margin:
        raise GeometryError(f"stress section {name}: span {s1 - s0:.3f} mm is not longer than twice "
                            f"stress_corner_margin ({margin} mm)")
    return Section(name, "band", o=tuple(map(float, o)), d=tuple(map(float, d)), n=tuple(map(float, n)),
                   s0=s0 + margin, s1=s1 - margin, depth=depth)


def stress_sections(rp) -> list[Section]:
    """The named sections of a built rotor (see the module docstring)."""
    p = rp.p
    m = p.stress_corner_margin
    L1, L2 = rp.layers
    out: list[Section] = []
    if not p.is_bridgeless:
        for L in (L1, L2):
            k = L.inp.k
            R = rp.r_lam - getattr(p, f"pm{k}_bridge_thickness")
            a3 = math.atan2(L.p("P3")[1], L.p("P3")[0])
            a4 = math.atan2(L.p("P4")[1], L.p("P4")[0])
            a0, a1 = _span_containing(_arc_spans(rp, R), a3 - 1e-9, a4 + 1e-9, f"bridge{k}")
            out.append(_sector(f"bridge{k}", R, rp.r_lam, a0, a1, m))
        # layer-1 centre rib: below the flat inner end of the layer-1 slot, down to the d-axis
        d = np.array([1.0, 0.0])
        o, length = _end_edge(rp, d, np.array(L1.p("P6")), np.array(L1.p("P5")), "rib1")
        out.append(_band("rib1", o, d, np.array([0.0, -1.0]), 0.0, length, float(o[1]), m))
        # layer-2 web: below the tilted inner end of the layer-2 slot, thickness pm2_web_width
        ti = L2.inp.tINNO
        d = np.array([math.cos(ti), -math.sin(ti)])            # along the end edge (end band rotated by -tilt)
        n = np.array([-math.sin(ti), -math.cos(ti)])           # into the web, towards the cavity
        o, length = _end_edge(rp, d, np.array(L2.p("P6")), np.array(L2.p("P5")), "web2")
        out.append(_band("web2", o, d, n, 0.0, length, p.pm2_web_width, m))
        out += _qhole_sections(rp)
    else:
        from .rotor import opening_angles

        r_in = rp.r_lam - p.lip_depth
        spans = _arc_spans(rp, r_in)
        for L in (L1, L2):
            k = L.inp.k
            lo, hi = opening_angles(p, L)
            pole = [sp for sp in spans if abs(sp[1] - lo) < 1e-9]
            q = [sp for sp in spans if abs(sp[0] - hi) < 1e-9]
            if len(pole) != 1 or len(q) != 1:
                raise GeometryError(f"stress section: lips of layer {k} not found next to its OD opening")
            out.append(_sector(f"lip{k}_pole", r_in, rp.r_lam, *pole[0], m))
            out.append(_sector(f"lip{k}_q", r_in, rp.r_lam, *q[0], m))
        out += _qhole_sections(rp)
    return out


def section_peaks(sections: list[Section], xy: np.ndarray, seqv: np.ndarray) -> dict[str, dict]:
    """Peak von Mises stress and its location for each section from nodal values.

    Regions that share a name form one section. Returns {name: dict(value, x, y, nodes)} in
    the order of first appearance; a section that contains no node raises."""
    masks: dict[str, np.ndarray] = {}
    for s in sections:
        m = s.contains(xy)
        masks[s.name] = masks[s.name] | m if s.name in masks else m
    out = {}
    for name, mask in masks.items():
        if not mask.any():
            raise GeometryError(f"stress section {name}: no lamination node inside (mesh too coarse?)")
        i = np.flatnonzero(mask)[np.argmax(seqv[mask])]
        out[name] = dict(value=float(seqv[i]), x=float(xy[i, 0]), y=float(xy[i, 1]), nodes=int(mask.sum()))
    return out
