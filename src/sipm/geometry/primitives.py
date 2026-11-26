"""Exact 2-D boundary segments and the solver-facing :class:`Loop`.

Two boundary representations are used:

* an *exact chain*: a closed tuple of :class:`Line`, :class:`Arc` and :class:`Hermite`
  segments, each running from ``start()`` to ``end()``. Every construction in the package is
  done on exact chains, and the reference fixtures from the Maxwell models are compared with
  them;
* a :class:`Loop`: a closed chain of straight lines and minor circular arcs only, the form the
  Maxwell and MAPDL deck writers consume. :func:`chain_to_loop` converts an exact chain:
  lines and arcs are copied exactly, a cubic Hermite curve (the curves of the central
  cavity) is replaced by tangent-continuous biarcs whose
  distance from the true curve is at most ``tol`` (default :data:`BIARC_TOL`), identically
  for both solvers.

Lengths in mm, angles in rad.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np

TWO_PI = 2.0 * math.pi
BIARC_TOL = 5e-4          # mm, maximum distance of the biarc replacement from a Hermite curve
ARC_SPLIT = math.radians(90.0)  # longest arc piece in a Loop (Loop arcs must be minor arcs)

Pt = tuple[float, float]


class GeometryError(ValueError):
    """Raised for infeasible geometry; the message names the violated constraint."""


def _p(v) -> Pt:
    return (float(v[0]), float(v[1]))


# ======================================================================================
# exact segments
# ======================================================================================
@dataclass(frozen=True)
class Line:
    p0: Pt
    p1: Pt
    kind = "line"

    def start(self) -> Pt:
        return self.p0

    def end(self) -> Pt:
        return self.p1

    def point(self, u: float) -> Pt:
        return (self.p0[0] + u * (self.p1[0] - self.p0[0]), self.p0[1] + u * (self.p1[1] - self.p0[1]))

    def length(self) -> float:
        return math.hypot(self.p1[0] - self.p0[0], self.p1[1] - self.p0[1])

    def tangent(self, u: float = 0.0) -> Pt:
        L = self.length()
        return ((self.p1[0] - self.p0[0]) / L, (self.p1[1] - self.p0[1]) / L)

    def reversed(self) -> Line:
        return Line(self.p1, self.p0)

    def transformed(self, M: np.ndarray, b: Sequence[float] = (0.0, 0.0)) -> Line:
        return Line(_affine(M, b, self.p0), _affine(M, b, self.p1))

    def area_term(self) -> float:
        """Integral of x dy - y dx along the segment (twice the signed-area contribution)."""
        return self.p0[0] * self.p1[1] - self.p1[0] * self.p0[1]

    def sample(self, n: int = 2) -> list[Pt]:
        return [self.point(k / (n - 1)) for k in range(n)]


@dataclass(frozen=True)
class Arc:
    """Circular arc: centre ``c``, radius ``r``, start angle ``a0``, signed sweep ``da``
    (counter-clockwise positive; ``|da| = 2 pi`` is a full circle)."""

    c: Pt
    r: float
    a0: float
    da: float
    kind = "arc"

    def point(self, u: float) -> Pt:
        a = self.a0 + u * self.da
        return (self.c[0] + self.r * math.cos(a), self.c[1] + self.r * math.sin(a))

    def start(self) -> Pt:
        return self.point(0.0)

    def end(self) -> Pt:
        return self.point(1.0)

    def length(self) -> float:
        return abs(self.da) * self.r

    def tangent(self, u: float = 0.0) -> Pt:
        a = self.a0 + u * self.da
        s = 1.0 if self.da > 0 else -1.0
        return (-math.sin(a) * s, math.cos(a) * s)

    def reversed(self) -> Arc:
        return Arc(self.c, self.r, self.a0 + self.da, -self.da)

    def transformed(self, M: np.ndarray, b: Sequence[float] = (0.0, 0.0)) -> Arc:
        c = _affine(M, b, self.c)
        p0 = _affine(M, b, self.start())
        det = float(np.linalg.det(M))
        a0 = math.atan2(p0[1] - c[1], p0[0] - c[0])
        return Arc(c, self.r, a0, self.da if det > 0 else -self.da)

    def area_term(self) -> float:
        cx, cy, r = self.c[0], self.c[1], self.r
        a0, a1 = self.a0, self.a0 + self.da
        return r * r * (a1 - a0) + cx * r * (math.sin(a1) - math.sin(a0)) - cy * r * (math.cos(a1) - math.cos(a0))

    def contains_angle(self, a: float, tol: float = 1e-12) -> float | None:
        """Parameter u in [0, 1] of polar angle ``a`` on this arc, or None if outside."""
        rel = (a - self.a0) if self.da > 0 else (self.a0 - a)
        rel %= TWO_PI
        span = abs(self.da)
        if rel <= span + tol:
            return min(rel / span, 1.0)
        if rel >= TWO_PI - tol:  # just before the start
            return 0.0
        return None

    def sample(self, n: int = 2) -> list[Pt]:
        return [self.point(k / (n - 1)) for k in range(n)]


@dataclass(frozen=True)
class Hermite:
    """Cubic Hermite curve p(t) = h00 pa + h10 ta + h01 pb + h11 tb, t in [0, 1] (the
    ``EquationBasedCurve`` form of the Maxwell models). ``forward=False`` runs it from pb to pa."""

    pa: Pt
    ta: Pt
    pb: Pt
    tb: Pt
    forward: bool = True
    kind = "hermite"

    def at_t(self, t: float) -> Pt:
        t2, t3 = t * t, t * t * t
        h00, h10, h01, h11 = 2 * t3 - 3 * t2 + 1, t3 - 2 * t2 + t, -2 * t3 + 3 * t2, t3 - t2
        return (h00 * self.pa[0] + h10 * self.ta[0] + h01 * self.pb[0] + h11 * self.tb[0],
                h00 * self.pa[1] + h10 * self.ta[1] + h01 * self.pb[1] + h11 * self.tb[1])

    def deriv_t(self, t: float) -> Pt:
        t2 = t * t
        d00, d10, d01, d11 = 6 * t2 - 6 * t, 3 * t2 - 4 * t + 1, -6 * t2 + 6 * t, 3 * t2 - 2 * t
        return (d00 * self.pa[0] + d10 * self.ta[0] + d01 * self.pb[0] + d11 * self.tb[0],
                d00 * self.pa[1] + d10 * self.ta[1] + d01 * self.pb[1] + d11 * self.tb[1])

    def _t(self, u: float) -> float:
        return u if self.forward else 1.0 - u

    def point(self, u: float) -> Pt:
        return self.at_t(self._t(u))

    def start(self) -> Pt:
        return self.point(0.0)

    def end(self) -> Pt:
        return self.point(1.0)

    def tangent(self, u: float = 0.0) -> Pt:
        d = self.deriv_t(self._t(u))
        L = math.hypot(*d)
        s = 1.0 if self.forward else -1.0
        return (s * d[0] / L, s * d[1] / L)

    def reversed(self) -> Hermite:
        return Hermite(self.pa, self.ta, self.pb, self.tb, not self.forward)

    def transformed(self, M: np.ndarray, b: Sequence[float] = (0.0, 0.0)) -> Hermite:
        return Hermite(_affine(M, b, self.pa), _linear(M, self.ta), _affine(M, b, self.pb),
                       _linear(M, self.tb), self.forward)

    def length(self, n: int = 400) -> float:
        pts = self.sample(n + 1)
        return float(sum(math.dist(pts[k], pts[k + 1]) for k in range(n)))

    def area_term(self) -> float:
        # x(t) y'(t) - y(t) x'(t) is a polynomial of degree 4 -> 3-point Gauss-Legendre is exact
        xs, ws = np.polynomial.legendre.leggauss(3)
        tot = 0.0
        for x, w in zip(xs, ws):
            t = 0.5 * (x + 1.0)
            p, d = self.at_t(t), self.deriv_t(t)
            tot += 0.5 * w * (p[0] * d[1] - p[1] * d[0])
        return tot if self.forward else -tot

    def sample(self, n: int = 2) -> list[Pt]:
        return [self.point(k / (n - 1)) for k in range(n)]


Segment = Line | Arc | Hermite
Chain = tuple  # tuple[Segment, ...], closed: seg[i].end() == seg[i+1].start()


def _affine(M: np.ndarray, b: Sequence[float], p: Sequence[float]) -> Pt:
    return (float(M[0, 0] * p[0] + M[0, 1] * p[1] + b[0]), float(M[1, 0] * p[0] + M[1, 1] * p[1] + b[1]))


def _linear(M: np.ndarray, v: Sequence[float]) -> Pt:
    return (float(M[0, 0] * v[0] + M[0, 1] * v[1]), float(M[1, 0] * v[0] + M[1, 1] * v[1]))


MIRROR_X = np.diag([1.0, -1.0])


def rotation(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]])


# ---- chain helpers ------------------------------------------------------------------------
def chain_area(chain: Iterable[Segment]) -> float:
    """Signed area (counter-clockwise positive)."""
    return 0.5 * sum(s.area_term() for s in chain)


def chain_reversed(chain: Sequence[Segment]) -> Chain:
    return tuple(s.reversed() for s in reversed(chain))


def chain_ccw(chain: Sequence[Segment]) -> Chain:
    return tuple(chain) if chain_area(chain) > 0 else chain_reversed(chain)


def chain_transformed(chain: Sequence[Segment], M: np.ndarray, b: Sequence[float] = (0.0, 0.0)) -> Chain:
    """Apply x -> M x + b; a reflection reverses the chain so its orientation is kept."""
    out = tuple(s.transformed(M, b) for s in chain)
    return chain_reversed(out) if np.linalg.det(M) < 0 else out


def chain_mirror_x(chain: Sequence[Segment]) -> Chain:
    return chain_transformed(chain, MIRROR_X)


def chain_polyline(chain: Sequence[Segment], step_deg: float = 1.0, curve_n: int = 64) -> np.ndarray:
    """Dense polyline (without repeated end point) for plotting and sampled checks."""
    out: list[Pt] = []
    for s in chain:
        if s.kind == "line":
            out.append(s.p0)
        elif s.kind == "arc":
            n = max(2, int(math.ceil(abs(math.degrees(s.da)) / step_deg)) + 1)
            out.extend(s.sample(n)[:-1])
        else:
            out.extend(s.sample(curve_n + 1)[:-1])
    return np.array(out, dtype=float)


def chain_closure_error(chain: Sequence[Segment]) -> float:
    n = len(chain)
    return max(math.dist(chain[i].end(), chain[(i + 1) % n].start()) for i in range(n))


def polygon_chain(points: Sequence[Sequence[float]]) -> Chain:
    pts = [_p(q) for q in points]
    return tuple(Line(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts)))


def arc_from(c: Sequence[float], p_start: Sequence[float], p_end: Sequence[float], ccw: bool) -> Arc:
    """Arc around ``c`` from ``p_start`` to ``p_end`` in the given direction (radius from p_start)."""
    a0 = math.atan2(p_start[1] - c[1], p_start[0] - c[0])
    a1 = math.atan2(p_end[1] - c[1], p_end[0] - c[0])
    da = (a1 - a0) % TWO_PI if ccw else -((a0 - a1) % TWO_PI)
    return Arc(_p(c), math.hypot(p_start[0] - c[0], p_start[1] - c[1]), a0, da)


def circle_chain(c: Sequence[float], r: float) -> Chain:
    """Full circle, counter-clockwise, starting at angle 0."""
    return (Arc(_p(c), float(r), 0.0, TWO_PI),)


def annular_sector(r_in: float, r_out: float, a0: float, a1: float) -> Chain:
    """Region r_in <= r <= r_out, a0 <= theta <= a1 (counter-clockwise chain); r_in = 0 gives a
    circular sector."""
    o = (0.0, 0.0)
    if r_in <= 0:
        p1, p2 = _pol(r_out, a0), _pol(r_out, a1)
        return (Line(o, p1), Arc(o, r_out, a0, a1 - a0), Line(p2, o))
    return (Line(_pol(r_in, a0), _pol(r_out, a0)), Arc(o, r_out, a0, a1 - a0),
            Line(_pol(r_out, a1), _pol(r_in, a1)), Arc(o, r_in, a1, a0 - a1))


def _pol(r: float, th: float) -> Pt:
    return (r * math.cos(th), r * math.sin(th))


# ======================================================================================
# fillets on exact chains
# ======================================================================================
def fillet_chain(chain: Sequence[Segment], idx: int, r: float) -> Chain:
    """Round the corner at the start of segment ``idx`` (between segment idx-1 and idx) with
    radius ``r``; both segments must be lines or arcs. The fillet lies inside the corner
    angle (< 180 deg) on the side of the region, as Maxwell's vertex fillet does."""
    n = len(chain)
    a, b = chain[(idx - 1) % n], chain[idx]
    v = b.start()
    ta, tb = a.tangent(1.0), b.tangent(0.0)
    turn = ta[0] * tb[1] - ta[1] * tb[0]          # > 0: left turn (convex corner of a CCW region)
    if abs(turn) < 1e-12:
        raise GeometryError(f"no corner to fillet at ({v[0]:.4f}, {v[1]:.4f})")
    side = 1.0 if turn > 0 else -1.0              # fillet centre is on the left of both segments if convex
    c = _fillet_centre(a, b, r, side, v)
    ta_pt = _foot(a, c)
    tb_pt = _foot(b, c)
    new_a = _trim(a, end=ta_pt)
    new_b = _trim(b, start=tb_pt)
    arc = arc_from(c, ta_pt, tb_pt, ccw=side > 0)
    if abs(arc.da) > math.pi:
        raise GeometryError(f"fillet radius {r} mm does not fit at ({v[0]:.4f}, {v[1]:.4f})")
    for s, ref in ((new_a, a), (new_b, b)):
        if s.length() <= 1e-9 or s.length() > ref.length() + 1e-9:
            raise GeometryError(f"fillet radius {r} mm does not fit at corner ({v[0]:.4f}, {v[1]:.4f}) mm")
    out = list(chain)
    out[(idx - 1) % n] = new_a
    out[idx] = new_b
    out.insert(idx, arc)
    return tuple(out)


def _offset_carrier(s: Segment, d: float, side: float):
    """Carrier of segment ``s`` offset by ``d`` to its left (side=+1) or right (side=-1)."""
    if s.kind == "line":
        t = s.tangent()
        nrm = (-t[1] * side, t[0] * side)
        return ("line", (s.p0[0] + d * nrm[0], s.p0[1] + d * nrm[1]), t)
    if s.kind == "arc":
        # left of a CCW arc is towards the centre
        towards_centre = (s.da > 0) == (side > 0)
        return ("circle", s.c, s.r - d if towards_centre else s.r + d)
    raise GeometryError("fillets are only supported between lines and arcs")


def _fillet_centre(a: Segment, b: Segment, r: float, side: float, near: Pt) -> Pt:
    ca, cb = _offset_carrier(a, r, side), _offset_carrier(b, r, side)
    cands = _carrier_intersections(ca, cb)
    if not cands:
        raise GeometryError(f"fillet radius {r} mm does not fit at ({near[0]:.4f}, {near[1]:.4f})")
    return min(cands, key=lambda q: math.dist(q, near))


def _carrier_intersections(ca, cb) -> list[Pt]:
    if ca[0] == "line" and cb[0] == "line":
        q = _line_line(ca[1], ca[2], cb[1], cb[2])
        return [] if q is None else [q]
    if ca[0] == "circle" and cb[0] == "line":
        ca, cb = cb, ca
    if ca[0] == "line" and cb[0] == "circle":
        return [pt for pt, _ in _line_circle_full(ca[1], ca[2], cb[1], cb[2])]
    return _circle_circle(ca[1], ca[2], cb[1], cb[2])


def _line_line(p: Pt, d: Pt, q: Pt, e: Pt) -> Pt | None:
    den = d[0] * e[1] - d[1] * e[0]
    if abs(den) < 1e-15:
        return None
    t = ((q[0] - p[0]) * e[1] - (q[1] - p[1]) * e[0]) / den
    return (p[0] + t * d[0], p[1] + t * d[1])


def _line_circle_full(p: Pt, d: Pt, c: Pt, r: float) -> list[tuple[Pt, float]]:
    """Points of the infinite line p + s d (d unit) on the circle; returns (point, s)."""
    fx, fy = p[0] - c[0], p[1] - c[1]
    b = fx * d[0] + fy * d[1]
    cc = fx * fx + fy * fy - r * r
    disc = b * b - cc
    if disc < -1e-12 * max(1.0, r * r):
        return []
    disc = max(disc, 0.0)
    sq = math.sqrt(disc)
    out = []
    for s in ((-b - sq, -b + sq) if sq > 0 else (-b,)):
        out.append(((p[0] + s * d[0], p[1] + s * d[1]), s))
    return out


def _circle_circle(c1: Pt, r1: float, c2: Pt, r2: float) -> list[Pt]:
    dx, dy = c2[0] - c1[0], c2[1] - c1[1]
    d = math.hypot(dx, dy)
    if d < 1e-15:
        return []
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    h2 = r1 * r1 - a * a
    if h2 < -1e-12 * max(1.0, r1 * r1):
        return []
    h = math.sqrt(max(h2, 0.0))
    mx, my = c1[0] + a * dx / d, c1[1] + a * dy / d
    if h == 0:
        return [(mx, my)]
    return [(mx - h * dy / d, my + h * dx / d), (mx + h * dy / d, my - h * dx / d)]


def _foot(s: Segment, p: Pt) -> Pt:
    if s.kind == "line":
        t = s.tangent()
        k = (p[0] - s.p0[0]) * t[0] + (p[1] - s.p0[1]) * t[1]
        return (s.p0[0] + k * t[0], s.p0[1] + k * t[1])
    dx, dy = p[0] - s.c[0], p[1] - s.c[1]
    L = math.hypot(dx, dy)
    return (s.c[0] + s.r * dx / L, s.c[1] + s.r * dy / L)


def _trim(s: Segment, start: Pt | None = None, end: Pt | None = None) -> Segment:
    p0 = start if start is not None else s.start()
    p1 = end if end is not None else s.end()
    if s.kind == "line":
        return Line(_p(p0), _p(p1))
    return arc_from(s.c, p0, p1, ccw=s.da > 0)


def offset_line_segment(chain: Sequence[Segment], idx: int, d: float) -> Chain:
    """Move straight segment ``idx`` of a polygon chain by ``d`` along its left normal (into a
    counter-clockwise region for d > 0) and re-intersect it with its straight neighbours
    (Maxwell ``MoveEdges``)."""
    n = len(chain)
    s = chain[idx]
    if s.kind != "line" or chain[idx - 1].kind != "line" or chain[(idx + 1) % n].kind != "line":
        raise GeometryError("offset_line_segment needs straight neighbours")
    t = s.tangent()
    nrm = (-t[1], t[0])
    p = (s.p0[0] + d * nrm[0], s.p0[1] + d * nrm[1])
    prev, nxt = chain[idx - 1], chain[(idx + 1) % n]
    q0 = _line_line(prev.p0, prev.tangent(), p, t)
    q1 = _line_line(p, t, nxt.p0, nxt.tangent())
    if q0 is None or q1 is None:
        raise GeometryError("offset edge is parallel to a neighbour")
    out = list(chain)
    out[idx - 1] = Line(prev.p0, q0)
    out[idx] = Line(q0, q1)
    out[(idx + 1) % n] = Line(q1, nxt.p1)
    return tuple(out)


# ======================================================================================
# biarc replacement of Hermite curves
# ======================================================================================
def _biarc(p0: Pt, t0: Pt, p1: Pt, t1: Pt) -> list[Segment]:
    """Equal-tangent-length biarc from p0 (unit tangent t0) to p1 (unit tangent t1)."""
    vx, vy = p1[0] - p0[0], p1[1] - p0[1]
    tx, ty = t0[0] + t1[0], t0[1] + t1[1]
    vv = vx * vx + vy * vy
    vt = vx * tx + vy * ty
    tt = tx * tx + ty * ty
    a = tt - 4.0
    if abs(a) < 1e-12:
        if abs(vt) < 1e-15:
            return [Line(p0, p1)]
        dd = vv / (2.0 * vt)
    else:
        disc = vt * vt - a * vv
        dd = (vt + math.sqrt(max(disc, 0.0))) / a
        if dd <= 0:
            dd = (vt - math.sqrt(max(disc, 0.0))) / a
    j = (0.5 * (p0[0] + dd * t0[0] + p1[0] - dd * t1[0]), 0.5 * (p0[1] + dd * t0[1] + p1[1] - dd * t1[1]))
    return [_arc_tangent(p0, t0, j), _arc_tangent(p1, (-t1[0], -t1[1]), j, reverse=True)]


def _arc_tangent(p: Pt, t: Pt, q: Pt, reverse: bool = False) -> Segment:
    """Arc leaving ``p`` with tangent ``t`` and ending at ``q`` (line if collinear). With
    ``reverse`` the arc is returned running from q to p."""
    nx, ny = -t[1], t[0]
    wx, wy = p[0] - q[0], p[1] - q[1]
    den = 2.0 * (nx * wx + ny * wy)
    ww = wx * wx + wy * wy
    if abs(den) < 1e-12 * math.sqrt(ww) or ww / abs(den) > 1e6:
        seg: Segment = Line(p, q)
    else:
        s = -ww / den
        c = (p[0] + s * nx, p[1] + s * ny)
        seg = arc_from(c, p, q, ccw=s > 0)
    return seg.reversed() if reverse else seg


def _seg_distance(seg: Segment, q: Pt) -> float:
    if seg.kind == "line":
        dx, dy = seg.p1[0] - seg.p0[0], seg.p1[1] - seg.p0[1]
        L2 = dx * dx + dy * dy
        u = 0.0 if L2 == 0 else max(0.0, min(1.0, ((q[0] - seg.p0[0]) * dx + (q[1] - seg.p0[1]) * dy) / L2))
        return math.dist(q, (seg.p0[0] + u * dx, seg.p0[1] + u * dy))
    a = math.atan2(q[1] - seg.c[1], q[0] - seg.c[0])
    if seg.contains_angle(a, 1e-9) is not None:
        return abs(math.hypot(q[0] - seg.c[0], q[1] - seg.c[1]) - seg.r)
    return min(math.dist(q, seg.start()), math.dist(q, seg.end()))


def hermite_to_biarcs(h: Hermite, tol: float = BIARC_TOL, depth: int = 0) -> list[Segment]:
    """Replace a Hermite curve by biarcs within ``tol`` mm (recursive bisection in t)."""
    return _biarc_range(h, 0.0, 1.0, tol, 0) if h.forward else \
        [s.reversed() for s in reversed(_biarc_range(h, 0.0, 1.0, tol, 0))]


def _unit(v: Pt) -> Pt:
    L = math.hypot(*v)
    return (v[0] / L, v[1] / L)


def _biarc_range(h: Hermite, ta: float, tb: float, tol: float, depth: int) -> list[Segment]:
    p0, p1 = h.at_t(ta), h.at_t(tb)
    segs = _biarc(p0, _unit(h.deriv_t(ta)), p1, _unit(h.deriv_t(tb)))
    err = 0.0
    for k in range(1, 16):
        q = h.at_t(ta + (tb - ta) * k / 16)
        err = max(err, min(_seg_distance(s, q) for s in segs))
    if (err <= tol and all(s.kind == "line" or abs(s.da) < math.radians(120) for s in segs)) or depth > 18:
        if err > tol:
            raise GeometryError(f"biarc replacement of a Hermite curve misses tolerance ({err:.2e} mm)")
        return segs
    tm = 0.5 * (ta + tb)
    return _biarc_range(h, ta, tm, tol, depth + 1) + _biarc_range(h, tm, tb, tol, depth + 1)


def hermite_max_deviation(h: Hermite, segs: Sequence[Segment], n: int = 2000) -> float:
    """Largest distance of ``n`` points of the true curve from the replacement chain."""
    return max(min(_seg_distance(s, h.at_t(k / n)) for s in segs) for k in range(n + 1))


# ======================================================================================
# solver-facing loop
# ======================================================================================
@dataclass
class Loop:
    """Closed boundary of lines and minor arcs. Segment i runs from pts[i] to pts[(i+1) % n];
    arcs[i] is None for a straight line or the (x, y) centre of a minor circular arc."""

    pts: np.ndarray
    arcs: list = field(default_factory=list)

    def __post_init__(self):
        self.pts = np.asarray(self.pts, dtype=float).reshape(-1, 2)
        if not self.arcs:
            self.arcs = [None] * len(self.pts)
        if len(self.arcs) != len(self.pts):
            raise ValueError("arcs must have one entry per segment")
        self.arcs = [None if c is None else np.asarray(c, float) for c in self.arcs]
        for i, c in enumerate(self.arcs):
            if c is not None:
                p, q = self.pts[i], self.pts[(i + 1) % len(self.pts)]
                r1, r2 = np.hypot(*(p - c)), np.hypot(*(q - c))
                if abs(r1 - r2) > 1e-6 * max(1.0, r1):
                    raise ValueError(f"arc {i}: end points not equidistant from centre ({r1} vs {r2})")

    def __len__(self):
        return len(self.pts)

    def segments(self):
        n = len(self.pts)
        for i in range(n):
            yield self.pts[i], self.pts[(i + 1) % n], self.arcs[i]

    @staticmethod
    def _sweep(p, q, c):
        a0 = math.atan2(p[1] - c[1], p[0] - c[0])
        a1 = math.atan2(q[1] - c[1], q[0] - c[0])
        d = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi  # minor arc in (-pi, pi]
        return a0, d

    def arc_mid(self, i):
        p, q, c = self.pts[i], self.pts[(i + 1) % len(self.pts)], self.arcs[i]
        a0, d = self._sweep(p, q, c)
        r = np.hypot(*(p - c))
        return c + r * np.array([math.cos(a0 + d / 2), math.sin(a0 + d / 2)])

    def signed_area(self) -> float:
        a = 0.0
        for p, q, c in self.segments():
            a += 0.5 * (p[0] * q[1] - q[0] * p[1])
            if c is not None:
                _, d = self._sweep(p, q, c)
                r = np.hypot(*(p - c))
                a += 0.5 * r * r * (d - math.sin(d))
        return a

    def area(self) -> float:
        return abs(self.signed_area())

    def discretize(self, max_step_deg: float = 1.0) -> np.ndarray:
        """Polyline (without repeated end point) approximating the loop."""
        out = []
        for p, q, c in self.segments():
            out.append(p)
            if c is not None:
                a0, d = self._sweep(p, q, c)
                r = np.hypot(*(p - c))
                n = max(2, int(math.ceil(abs(math.degrees(d)) / max_step_deg)))
                for k in range(1, n):
                    t = a0 + d * k / n
                    out.append(c + r * np.array([math.cos(t), math.sin(t)]))
        return np.array(out)

    def transformed(self, rot_deg=0.0, mirror_x=False) -> Loop:
        """Mirror about the x axis (optional) then rotate about the origin."""
        th = math.radians(rot_deg)
        R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        M = np.diag([1.0, -1.0]) if mirror_x else np.eye(2)
        T = R @ M
        pts = self.pts @ T.T
        arcs = [None if c is None else T @ c for c in self.arcs]
        if mirror_x:  # keep orientation (counter-clockwise) after a reflection
            pts, arcs = _reverse(pts, arcs)
        return Loop(pts, arcs)

    def ccw(self) -> Loop:
        if self.signed_area() < 0:
            p, a = _reverse(self.pts, self.arcs)
            return Loop(p, a)
        return self


def _reverse(pts, arcs):
    n = len(pts)
    newpts = pts[::-1].copy()
    # segment i of the reversed loop runs newpts[i] -> newpts[i+1] = pts[n-1-i] -> pts[n-2-i],
    # which is original segment (n-2-i) traversed backwards.
    newarcs = [arcs[(n - 2 - i) % n] for i in range(n)]
    return newpts, newarcs


def polygon(points) -> Loop:
    return Loop(np.asarray(points, float))


def chain_to_loop(chain: Sequence[Segment], tol: float = BIARC_TOL) -> Loop:
    """Solver form of an exact chain (lines and minor arcs; Hermite curves -> biarcs)."""
    pts: list[Pt] = []
    arcs: list = []
    for s in chain:
        pieces: list[Segment]
        if s.kind == "hermite":
            pieces = hermite_to_biarcs(s, tol)
        elif s.kind == "arc" and abs(s.da) > ARC_SPLIT + 1e-12:
            n = int(math.ceil(abs(s.da) / ARC_SPLIT))
            pieces = [Arc(s.c, s.r, s.a0 + s.da * k / n, s.da / n) for k in range(n)]
        else:
            pieces = [s]
        for q in pieces:
            pts.append(q.start())
            arcs.append(None if q.kind == "line" else q.c)
    return Loop(np.array(pts), arcs)
