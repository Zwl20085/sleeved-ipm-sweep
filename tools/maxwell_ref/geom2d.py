"""Exact 2-D boundary geometry: lines, circular arcs and parametric curves.

Coordinates are in mm.  A *sheet* is a set of faces (outer loop CCW, holes CW),
a *wire* a set of segment chains.  Every segment carries

* ``tag``   - a provenance name stable across re-evaluation (which primitive
              edge it comes from),
* ``ids``   - AEDT edge IDs attached to it,
* ``v0/v1`` - AEDT vertex IDs at its start / end.

Boolean operations are computed with shapely on densely sampled loops and the
exact boundary is then *recovered*: every output edge is mapped back to the
exact source segment it lies on, consecutive pieces are merged and every
junction is recomputed as the exact intersection (or exact shared end point)
of the two exact source carriers.  The result is again made only of exact
lines/arcs/curves; the sampling only decides topology.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable, Sequence

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.polygon import orient

TWO_PI = 2.0 * math.pi
ARC_GRID = 14400          # global angular sampling grid (0.025 deg)
CURVE_SAMPLES = 600       # samples per parametric curve
SNAP = 2e-6               # mm: junction snapping to an exact end point
ONLINE = 1e-7             # mm: output edge lies on a sampled source chord
SLIVER_AREA = 1e-7        # mm^2: faces smaller than this are dropped (logged)
GRID = 1e-9               # mm: sampled coordinates and shapely results are snapped to this grid


def _q(v: float) -> float:
    return round(v / GRID) * GRID


def _key(p) -> tuple:
    return (round(p[0] / GRID), round(p[1] / GRID))


class GeomError(RuntimeError):
    pass


Pt = tuple


def _pt(x: float, y: float) -> Pt:
    return (float(x), float(y))


def dist(a: Pt, b: Pt) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Affine:
    """p' = M @ p + t"""
    a: float = 1.0
    b: float = 0.0
    c: float = 0.0
    d: float = 1.0
    tx: float = 0.0
    ty: float = 0.0

    def __call__(self, p: Pt) -> Pt:
        return _pt(self.a * p[0] + self.b * p[1] + self.tx, self.c * p[0] + self.d * p[1] + self.ty)

    def vec(self, v: Pt) -> Pt:
        return _pt(self.a * v[0] + self.b * v[1], self.c * v[0] + self.d * v[1])

    @property
    def det(self) -> float:
        return self.a * self.d - self.b * self.c

    def then(self, other: "Affine") -> "Affine":
        """other o self"""
        return Affine(other.a * self.a + other.b * self.c, other.a * self.b + other.b * self.d,
                      other.c * self.a + other.d * self.c, other.c * self.b + other.d * self.d,
                      other.a * self.tx + other.b * self.ty + other.tx,
                      other.c * self.tx + other.d * self.ty + other.ty)

    def inverse(self) -> "Affine":
        det = self.det
        a, b, c, d = self.d / det, -self.b / det, -self.c / det, self.a / det
        return Affine(a, b, c, d, -(a * self.tx + b * self.ty), -(c * self.tx + d * self.ty))

    def as_list(self) -> list:
        return [[self.a, self.b, self.tx], [self.c, self.d, self.ty]]

    @staticmethod
    def rotation(angle: float, center: Pt = (0.0, 0.0)) -> "Affine":
        c, s = math.cos(angle), math.sin(angle)
        return Affine(c, -s, s, c, center[0] - c * center[0] + s * center[1], center[1] - s * center[0] - c * center[1])

    @staticmethod
    def translation(dx: float, dy: float) -> "Affine":
        return Affine(tx=dx, ty=dy)

    @staticmethod
    def mirror(base: Pt, normal: Pt) -> "Affine":
        nx, ny = normal
        n = math.hypot(nx, ny)
        if n == 0:
            raise GeomError("mirror normal of zero length")
        nx, ny = nx / n, ny / n
        a, b, c, d = 1 - 2 * nx * nx, -2 * nx * ny, -2 * nx * ny, 1 - 2 * ny * ny
        return Affine(a, b, c, d, base[0] - (a * base[0] + b * base[1]), base[1] - (c * base[0] + d * base[1]))


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Seg:
    tag: str
    ids: frozenset = frozenset()
    v0: frozenset = frozenset()
    v1: frozenset = frozenset()

    def with_meta(self, **kw) -> "Seg":
        return replace(self, **kw)


@dataclass(frozen=True)
class Line(Seg):
    p0: Pt = (0.0, 0.0)
    p1: Pt = (0.0, 0.0)

    kind = "line"

    def start(self) -> Pt:
        return self.p0

    def end(self) -> Pt:
        return self.p1

    def at(self, u: float) -> Pt:
        return _pt(self.p0[0] + u * (self.p1[0] - self.p0[0]), self.p0[1] + u * (self.p1[1] - self.p0[1]))

    def length(self) -> float:
        return dist(self.p0, self.p1)

    def sample(self) -> list[tuple[Pt, float]]:
        return [(self.p0, 0.0)]

    def transform(self, T: Affine) -> "Line":
        return replace(self, p0=T(self.p0), p1=T(self.p1))

    def reversed(self) -> "Line":
        return replace(self, p0=self.p1, p1=self.p0, v0=self.v1, v1=self.v0)

    def project(self, p: Pt) -> float:
        dx, dy = self.p1[0] - self.p0[0], self.p1[1] - self.p0[1]
        L2 = dx * dx + dy * dy
        return ((p[0] - self.p0[0]) * dx + (p[1] - self.p0[1]) * dy) / L2

    def sub(self, p_start: Pt, p_end: Pt) -> "Line":
        return replace(self, p0=p_start, p1=p_end, v0=frozenset(), v1=frozenset())

    def area_term(self) -> float:  # contribution to 2*signed area (shoelace)
        return self.p0[0] * self.p1[1] - self.p1[0] * self.p0[1]

    def tangent_out(self) -> Pt:
        L = self.length()
        return _pt((self.p1[0] - self.p0[0]) / L, (self.p1[1] - self.p0[1]) / L)

    tangent_in = tangent_out

    def describe(self) -> dict:
        return {"type": "line", "p0": list(self.p0), "p1": list(self.p1)}


@dataclass(frozen=True)
class Arc(Seg):
    c: Pt = (0.0, 0.0)
    r: float = 1.0
    a0: float = 0.0
    da: float = 0.0   # signed sweep, CCW positive

    kind = "arc"

    def at(self, u: float) -> Pt:
        a = self.a0 + u * self.da
        return _pt(self.c[0] + self.r * math.cos(a), self.c[1] + self.r * math.sin(a))

    def start(self) -> Pt:
        return self.at(0.0)

    def end(self) -> Pt:
        return self.at(1.0)

    @property
    def closed(self) -> bool:
        return abs(abs(self.da) - TWO_PI) < 1e-12

    def length(self) -> float:
        return abs(self.da) * self.r

    def sample(self) -> list[tuple[Pt, float]]:
        step = TWO_PI / ARC_GRID
        out = [(self.at(0.0), 0.0)]
        a_lo, a_hi = sorted((self.a0, self.a0 + self.da))
        k0, k1 = math.floor(a_lo / step) + 1, math.ceil(a_hi / step) - 1
        ks = range(k0, k1 + 1) if self.da > 0 else range(k1, k0 - 1, -1)
        for k in ks:
            a = k * step
            u = (a - self.a0) / self.da
            if 1e-9 < u < 1 - 1e-9:
                out.append((_pt(self.c[0] + self.r * math.cos(a), self.c[1] + self.r * math.sin(a)), u))
        return out

    def transform(self, T: Affine) -> "Arc":
        det = T.det
        if abs(abs(det) - 1.0) > 1e-12:
            raise GeomError("non-rigid transform of an arc")
        c = T(self.c)
        p0 = T(self.at(0.0))
        a0 = math.atan2(p0[1] - c[1], p0[0] - c[0])
        da = self.da if det > 0 else -self.da
        return replace(self, c=c, a0=a0, da=da)

    def reversed(self) -> "Arc":
        return replace(self, a0=self.a0 + self.da, da=-self.da, v0=self.v1, v1=self.v0)

    def project(self, p: Pt) -> float:
        a = math.atan2(p[1] - self.c[1], p[0] - self.c[0])
        rel = a - self.a0
        if self.da > 0:
            rel = rel % TWO_PI
            if rel > self.da + (TWO_PI - self.da) / 2:
                rel -= TWO_PI
        else:
            rel = -((-rel) % TWO_PI)
            if rel < self.da - (TWO_PI + self.da) / 2:
                rel += TWO_PI
        return rel / self.da

    def angle_of(self, p: Pt) -> float:
        return self.a0 + self.project(p) * self.da

    def sub(self, p_start: Pt, p_end: Pt, reverse: bool = False, full: bool = False) -> "Arc":
        a_s = math.atan2(p_start[1] - self.c[1], p_start[0] - self.c[0])
        a_e = math.atan2(p_end[1] - self.c[1], p_end[0] - self.c[0])
        sgn = -1.0 if (self.da < 0) != reverse else 1.0
        sweep = (a_e - a_s) % TWO_PI if sgn > 0 else -((a_s - a_e) % TWO_PI)
        if full:
            sweep = sgn * TWO_PI
        elif abs(sweep) < 1e-15:
            sweep = 0.0
        return replace(self, a0=a_s, da=sweep, v0=frozenset(), v1=frozenset())

    def area_term(self) -> float:
        # integral of (x dy - y dx) along the arc
        cx, cy, r, a0, a1 = self.c[0], self.c[1], self.r, self.a0, self.a0 + self.da
        return (r * r * (a1 - a0) + cx * r * (math.sin(a1) - math.sin(a0))
                - cy * r * (math.cos(a1) - math.cos(a0)))

    def tangent_at(self, u: float) -> Pt:
        a = self.a0 + u * self.da
        s = 1.0 if self.da > 0 else -1.0
        return _pt(-math.sin(a) * s, math.cos(a) * s)

    def tangent_in(self) -> Pt:
        return self.tangent_at(1.0)

    def tangent_out(self) -> Pt:
        return self.tangent_at(0.0)

    def describe(self) -> dict:
        return {"type": "arc", "center": list(self.c), "radius": self.r,
                "start_angle_deg": math.degrees(self.a0), "sweep_deg": math.degrees(self.da),
                "p0": list(self.at(0.0)), "pmid": list(self.at(0.5)), "p1": list(self.at(1.0))}


@dataclass(frozen=True)
class Curve(Seg):
    """Parametric curve p(t) = T(f(t)), t in [t0, t1] (direction t0 -> t1)."""
    f: Callable = None  # type: ignore[assignment]
    t0: float = 0.0
    t1: float = 1.0
    T: Affine = Affine()
    desc: tuple = ()  # (("x", expr), ("y", expr), ...) exact description

    kind = "curve"

    def point_t(self, t: float) -> Pt:
        return self.T(self.f(t))

    def at(self, u: float) -> Pt:
        return self.point_t(self.t0 + u * (self.t1 - self.t0))

    def start(self) -> Pt:
        return self.at(0.0)

    def end(self) -> Pt:
        return self.at(1.0)

    def sample(self) -> list[tuple[Pt, float]]:
        return [(self.at(k / CURVE_SAMPLES), k / CURVE_SAMPLES) for k in range(CURVE_SAMPLES)]

    def length(self) -> float:
        pts = [self.at(k / 2000) for k in range(2001)]
        return sum(dist(pts[k], pts[k + 1]) for k in range(2000))

    def transform(self, T: Affine) -> "Curve":
        return replace(self, T=self.T.then(T))

    def reversed(self) -> "Curve":
        return replace(self, t0=self.t1, t1=self.t0, v0=self.v1, v1=self.v0)

    def project(self, p: Pt) -> float:
        us = np.linspace(0.0, 1.0, 2001)
        pts = np.array([self.at(u) for u in us])
        k = int(np.argmin(np.hypot(pts[:, 0] - p[0], pts[:, 1] - p[1])))
        lo, hi = us[max(k - 1, 0)], us[min(k + 1, 2000)]
        for _ in range(60):  # golden-section refinement
            m1, m2 = lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618
            if dist(self.at(m1), p) < dist(self.at(m2), p):
                hi = m2
            else:
                lo = m1
        return 0.5 * (lo + hi)

    def sub(self, p_start: Pt, p_end: Pt, reverse: bool = False) -> "Curve":
        us, ue = self.project(p_start), self.project(p_end)
        ts, te = self.t0 + us * (self.t1 - self.t0), self.t0 + ue * (self.t1 - self.t0)
        return replace(self, t0=ts, t1=te, v0=frozenset(), v1=frozenset())

    def area_term(self) -> float:
        # Gauss-Legendre on 64 panels x 8 points of x dy - y dx
        xs, ws = np.polynomial.legendre.leggauss(8)
        tot = 0.0
        n = 64
        h = 1e-7
        for k in range(n):
            a, b = k / n, (k + 1) / n
            for x, w in zip(xs, ws):
                u = 0.5 * (b - a) * x + 0.5 * (a + b)
                p = self.at(u)
                pa, pb = self.at(max(u - h, 0.0)), self.at(min(u + h, 1.0))
                du = min(u + h, 1.0) - max(u - h, 0.0)
                dx, dy = (pb[0] - pa[0]) / du, (pb[1] - pa[1]) / du
                tot += 0.5 * (b - a) * w * (p[0] * dy - p[1] * dx)
        return tot

    def tangent_in(self) -> Pt:
        a, b = self.at(1.0 - 1e-6), self.at(1.0)
        L = dist(a, b)
        return _pt((b[0] - a[0]) / L, (b[1] - a[1]) / L)

    def tangent_out(self) -> Pt:
        a, b = self.at(0.0), self.at(1e-6)
        L = dist(a, b)
        return _pt((b[0] - a[0]) / L, (b[1] - a[1]) / L)

    def describe(self) -> dict:
        d = {"type": "curve", "t_start": self.t0, "t_end": self.t1, "transform": self.T.as_list(),
             "p0": list(self.start()), "pmid": list(self.at(0.5)), "p1": list(self.end()),
             "samples_uniform_t_41": [list(self.at(k / 40)) for k in range(41)]}
        d.update({k: v for k, v in self.desc})
        return d


Loop = tuple  # tuple[Seg, ...]


@dataclass(frozen=True)
class Face:
    outer: Loop
    holes: tuple = ()

    def loops(self) -> list:
        return [self.outer, *self.holes]


@dataclass(frozen=True)
class Sheet:
    faces: tuple

    kind = "sheet"

    def transform(self, T: Affine) -> "Sheet":
        flip = T.det < 0
        out = []
        for f in self.faces:
            loops = [tuple(s.transform(T) for s in lp) for lp in f.loops()]
            if flip:  # keep outer CCW / holes CW
                loops = [reverse_loop(lp) for lp in loops]
            out.append(Face(loops[0], tuple(loops[1:])))
        return Sheet(tuple(out))

    def segments(self) -> list:
        return [s for f in self.faces for lp in f.loops() for s in lp]

    def loops(self) -> list:
        return [lp for f in self.faces for lp in f.loops()]


@dataclass(frozen=True)
class Wire:
    chains: tuple  # tuple[tuple[Seg,...], ...]

    kind = "wire"

    def transform(self, T: Affine) -> "Wire":
        return Wire(tuple(tuple(s.transform(T) for s in ch) for ch in self.chains))

    def segments(self) -> list:
        return [s for ch in self.chains for s in ch]

    def loops(self) -> list:
        return list(self.chains)


def reverse_loop(lp: Sequence[Seg]) -> Loop:
    return tuple(s.reversed() for s in reversed(lp))


def loop_area(lp: Sequence[Seg]) -> float:
    return 0.5 * sum(s.area_term() for s in lp)


def face_area(f: Face) -> float:
    return sum(loop_area(lp) for lp in f.loops())


def sheet_area(sh: Sheet) -> float:
    return sum(face_area(f) for f in sh.faces)


def loop_closure_error(lp: Sequence[Seg]) -> float:
    return max(dist(lp[i].end(), lp[(i + 1) % len(lp)].start()) for i in range(len(lp)))


# ---------------------------------------------------------------------------
# carriers and intersections
def _line_line(p, q, r, s) -> list:
    d1 = (q[0] - p[0], q[1] - p[1])
    d2 = (s[0] - r[0], s[1] - r[1])
    den = d1[0] * d2[1] - d1[1] * d2[0]
    scale = math.hypot(*d1) * math.hypot(*d2)
    if abs(den) < 1e-14 * scale:
        return []
    t = ((r[0] - p[0]) * d2[1] - (r[1] - p[1]) * d2[0]) / den
    return [_pt(p[0] + t * d1[0], p[1] + t * d1[1])]


def _line_circle(p, q, c, r) -> list:
    dx, dy = q[0] - p[0], q[1] - p[1]
    L = math.hypot(dx, dy)
    ux, uy = dx / L, dy / L
    # foot of perpendicular from c
    t = (c[0] - p[0]) * ux + (c[1] - p[1]) * uy
    fx, fy = p[0] + t * ux, p[1] + t * uy
    h2 = r * r - ((c[0] - fx) ** 2 + (c[1] - fy) ** 2)
    if h2 < -1e-12 * r * r:
        return []
    h = math.sqrt(max(h2, 0.0))
    return [_pt(fx + h * ux, fy + h * uy), _pt(fx - h * ux, fy - h * uy)]


def _circle_circle(c1, r1, c2, r2) -> list:
    d = dist(c1, c2)
    if d < 1e-14:
        return []
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    h2 = r1 * r1 - a * a
    if h2 < -1e-12 * r1 * r1:
        return []
    h = math.sqrt(max(h2, 0.0))
    ex, ey = (c2[0] - c1[0]) / d, (c2[1] - c1[1]) / d
    mx, my = c1[0] + a * ex, c1[1] + a * ey
    return [_pt(mx - h * ey, my + h * ex), _pt(mx + h * ey, my - h * ex)]


def _implicit(s: Seg):
    """F(p) with F=0 on the carrier, for line or circle."""
    if isinstance(s, Line):
        t = s.tangent_out()
        return lambda p: (p[0] - s.p0[0]) * t[1] - (p[1] - s.p0[1]) * t[0]
    if isinstance(s, Arc):
        return lambda p: dist(p, s.c) - s.r
    raise GeomError("no implicit form")


def intersect_near(a: Seg, b: Seg, near: Pt) -> Pt:
    """Exact intersection of the carriers of a and b closest to ``near``."""
    if isinstance(a, Curve) or isinstance(b, Curve):
        return _intersect_curve(a, b, near)
    if isinstance(a, Line) and isinstance(b, Line):
        c = _line_line(a.p0, a.p1, b.p0, b.p1)
    elif isinstance(a, Line) and isinstance(b, Arc):
        c = _line_circle(a.p0, a.p1, b.c, b.r)
    elif isinstance(a, Arc) and isinstance(b, Line):
        c = _line_circle(b.p0, b.p1, a.c, a.r)
    else:
        c = _circle_circle(a.c, a.r, b.c, b.r)
    if not c:
        raise GeomError(f"carriers of {a.tag} and {b.tag} do not intersect near {near}")
    return min(c, key=lambda p: dist(p, near))


def _intersect_curve(a: Seg, b: Seg, near: Pt) -> Pt:
    if isinstance(a, Curve) and isinstance(b, Curve):
        u, v = a.project(near), b.project(near)
        for _ in range(50):
            pa, pb = a.at(u), b.at(v)
            h = 1e-7
            ja = np.subtract(a.at(u + h), pa) / h
            jb = np.subtract(b.at(v + h), pb) / h
            J = np.array([[ja[0], -jb[0]], [ja[1], -jb[1]]])
            rhs = -np.array([pa[0] - pb[0], pa[1] - pb[1]])
            du, dv = np.linalg.solve(J, rhs)
            u, v = u + du, v + dv
            if abs(du) + abs(dv) < 1e-15:
                break
        return a.at(u)
    cur, other = (a, b) if isinstance(a, Curve) else (b, a)
    F = _implicit(other)
    u = cur.project(near)
    for _ in range(60):
        f0 = F(cur.at(u))
        h = 1e-8
        df = (F(cur.at(u + h)) - f0) / h
        if df == 0:
            break
        step = f0 / df
        u -= step
        if abs(step) < 1e-15:
            break
    return cur.at(u)


# ---------------------------------------------------------------------------
# sampling and recovery
@dataclass
class _Src:
    segs: list = field(default_factory=list)   # exact source segments
    chords: list = field(default_factory=list)  # (x0, y0, x1, y1, src_index)


def _sample_loop(lp: Sequence[Seg], src: _Src) -> list:
    coords = []
    for s in lp:
        k = len(src.segs)
        src.segs.append(s)
        pts = [(_q(p[0]), _q(p[1])) for p, _ in s.sample()]
        e = s.end()
        nxt = (_q(e[0]), _q(e[1]))
        for i, p in enumerate(pts):
            q = pts[i + 1] if i + 1 < len(pts) else nxt
            src.chords.append((p[0], p[1], q[0], q[1], k))
        coords.extend(pts)
    return coords


def to_shapely(sh: Sheet, src: _Src):
    polys = []
    for f in sh.faces:
        outer = _sample_loop(f.outer, src)
        holes = [_sample_loop(h, src) for h in f.holes]
        polys.append(Polygon(outer, holes))
    geom = shapely.union_all(polys, grid_size=GRID) if len(polys) > 1 else shapely.set_precision(polys[0], GRID)
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    return geom


class _ChordIndex:
    def __init__(self, src: _Src, cell: float = 0.25):
        self.src = src
        self.cell = cell
        self.grid: dict = {}
        arr = np.array([c[:4] for c in src.chords]) if src.chords else np.zeros((0, 4))
        self.arr = arr
        self.long: list = []
        for i, (x0, y0, x1, y1, _) in enumerate(src.chords):
            if (abs(x1 - x0) + abs(y1 - y0)) > 16 * cell:
                self.long.append(i)  # long chords are always tested
                continue
            for gx in range(math.floor(min(x0, x1) / cell), math.floor(max(x0, x1) / cell) + 1):
                for gy in range(math.floor(min(y0, y1) / cell), math.floor(max(y0, y1) / cell) + 1):
                    self.grid.setdefault((gx, gy), []).append(i)

    def carrier_class(self) -> list:
        """Map source index -> lowest index of a source with the same carrier
        (identical line or circle); coincident boundaries from different
        operands are then recovered as one piece."""
        if getattr(self, "_canon", None) is None:
            segs = self.src.segs
            canon = list(range(len(segs)))
            for j in range(len(segs)):
                for i in range(j):
                    if canon[i] == i and _carrier_equal(segs[i], segs[j]):
                        canon[j] = i
                        break
            self._canon = canon
        return self._canon

    @staticmethod
    def _pd(p, x0, y0, x1, y1) -> float:
        dx, dy = x1 - x0, y1 - y0
        L2 = dx * dx + dy * dy
        if L2 == 0:
            return math.hypot(p[0] - x0, p[1] - y0)
        t = max(0.0, min(1.0, ((p[0] - x0) * dx + (p[1] - y0) * dy) / L2))
        return math.hypot(p[0] - x0 - t * dx, p[1] - y0 - t * dy)

    def _incident(self) -> dict:
        if getattr(self, "_inc", None) is None:
            inc: dict = {}
            for i, (x0, y0, x1, y1, _) in enumerate(self.src.chords):
                inc.setdefault(_key((x0, y0)), []).append(i)
                inc.setdefault(_key((x1, y1)), []).append(i)
            self._inc = inc
        return self._inc

    def source_of(self, p: Pt, q: Pt) -> int:
        m = _pt(0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]))
        inc = self._incident()
        ip, iq = inc.get(_key(p)), inc.get(_key(q))
        if ip is not None or iq is not None:
            cand = sorted(set(ip or []) | set(iq or []))
            hit = self._best(cand, p, q, m)
            if hit is not None:
                return hit
        gx, gy = math.floor(m[0] / self.cell), math.floor(m[1] / self.cell)
        cand = sorted({i for dx in (-1, 0, 1) for dy in (-1, 0, 1) for i in self.grid.get((gx + dx, gy + dy), [])}
                      | set(self.long))
        hit = self._best(cand, p, q, m)
        if hit is None:
            raise GeomError(f"output edge {p}->{q} matches no source chord")
        return hit

    def _best(self, cand, p, q, m):
        best = None
        for i in cand:
            x0, y0, x1, y1, k = self.src.chords[i]
            e = max(self._pd(p, x0, y0, x1, y1), self._pd(q, x0, y0, x1, y1), self._pd(m, x0, y0, x1, y1))
            if e < ONLINE and (best is None or k < best[0]):
                best = (k, e)
        return None if best is None else best[0]


def _recover_ring(coords: list, index: _ChordIndex, log: list) -> Loop:
    pts = [tuple(c) for c in coords[:-1]]
    n = len(pts)
    segs = index.src.segs
    canon = index.carrier_class()
    raw = [index.source_of(pts[i], pts[(i + 1) % n]) for i in range(n)]
    srcs = [canon[k] for k in raw]
    # rotate so that the ring starts at a source change
    start = next((i for i in range(n) if srcs[i] != srcs[i - 1]), None)
    if start is None:  # single source around the whole ring (full circle)
        s = segs[srcs[0]]
        if isinstance(s, Arc) and s.closed:
            sign = 1.0 if Polygon(pts).exterior.is_ccw else -1.0
            if s.da * sign < 0:
                s = s.reversed()
            return (s,)
        raise GeomError(f"ring made of a single open source {s.tag}")
    pts = pts[start:] + pts[:start]
    srcs = srcs[start:] + srcs[:start]
    raw = raw[start:] + raw[:start]
    groups = []  # [src, first_pt, last_pt, n_edges, raw sources]
    for i in range(n):
        if groups and groups[-1][0] == srcs[i]:
            groups[-1][2] = pts[(i + 1) % n]
            groups[-1][3] += 1
            groups[-1][4].add(raw[i])
        else:
            groups.append([srcs[i], pts[i], pts[(i + 1) % n], 1, {raw[i]}])
    # exact junctions
    m = len(groups)
    members: dict = {}
    for k, c in enumerate(canon):
        members.setdefault(c, []).append(segs[k])
    junction = []
    for j in range(m):
        a, b = segs[groups[j - 1][0]], segs[groups[j][0]]
        p = groups[j][1]
        cands = [q for s in members[groups[j - 1][0]] + members[groups[j][0]] for q in (s.start(), s.end())]
        snap = min(cands, key=lambda c: dist(c, p))
        if dist(snap, p) < SNAP:
            junction.append(snap)
        else:
            x = intersect_near(a, b, p)
            if dist(x, p) > 1e-3:
                raise GeomError(f"junction of {a.tag}/{b.tag} at {p} recomputed {dist(x, p):.3g} mm away")
            junction.append(x)
    out = []
    for j in range(m):
        s = segs[groups[j][0]]
        p0, p1 = junction[j], junction[(j + 1) % m]
        if dist(p0, p1) < 1e-12:
            log.append(f"dropped zero-length piece of {s.tag}")
            continue
        piece = _sub_of(s, p0, p1, groups[j][1], groups[j][2])
        contrib = [segs[k] for k in sorted(groups[j][4])]
        # merged coincident edges keep one ID; AEDT retains the highest one
        # (see 'MergedEdges RetainedID' in Unite identities, e.g. coil unite)
        allids = frozenset().union(*(c.ids for c in contrib))
        ids = frozenset({max(allids)}) if allids else frozenset()
        # vertex IDs survive where the piece ends exactly on an original end point
        ends = [(q, vs) for c in contrib for q, vs in ((c.start(), c.v0), (c.end(), c.v1))]
        v0 = frozenset().union(*(vs for q, vs in ends if dist(q, p0) < SNAP))
        v1 = frozenset().union(*(vs for q, vs in ends if dist(q, p1) < SNAP))
        out.append(piece.with_meta(ids=ids, v0=v0, v1=v1))
    return tuple(_merge_collinear(out))


def _sub_of(s: Seg, p0: Pt, p1: Pt, a0: Pt, a1: Pt) -> Seg:
    if isinstance(s, Line):
        return s.sub(p0, p1)
    if isinstance(s, Arc):
        u0, u1 = s.project(a0), s.project(a1)
        rev = u1 < u0
        if s.closed:  # direction from the sampled path
            rev = False
            cross = (a0[0] - s.c[0]) * (a1[1] - s.c[1]) - (a0[1] - s.c[1]) * (a1[0] - s.c[0])
            rev = (cross > 0) != (s.da > 0) if abs(cross) > 1e-12 else False
        return s.sub(p0, p1, reverse=rev)
    if isinstance(s, Curve):
        return s.sub(p0, p1)
    raise GeomError(f"unknown segment {type(s)}")


def _carrier_equal(a: Seg, b: Seg) -> bool:
    """Same infinite line / same circle (direction ignored)."""
    if isinstance(a, Line) and isinstance(b, Line):
        ta = a.tangent_out()
        tb = b.tangent_out()
        cross = ta[0] * tb[1] - ta[1] * tb[0]
        off0 = abs((b.p0[0] - a.p0[0]) * ta[1] - (b.p0[1] - a.p0[1]) * ta[0])
        off1 = abs((b.p1[0] - a.p0[0]) * ta[1] - (b.p1[1] - a.p0[1]) * ta[0])
        return abs(cross) < 1e-11 and off0 < 1e-9 and off1 < 1e-9
    if isinstance(a, Arc) and isinstance(b, Arc):
        return dist(a.c, b.c) < 1e-9 and abs(a.r - b.r) < 1e-9
    return False


def _same_carrier(a: Seg, b: Seg) -> bool:
    if isinstance(a, Line) and isinstance(b, Line):
        ta, tb = a.tangent_out(), b.tangent_out()
        cross = ta[0] * tb[1] - ta[1] * tb[0]
        off = abs((b.p1[0] - a.p0[0]) * ta[1] - (b.p1[1] - a.p0[1]) * ta[0])
        return abs(cross) < 1e-12 and ta[0] * tb[0] + ta[1] * tb[1] > 0 and off < 1e-9
    if isinstance(a, Arc) and isinstance(b, Arc):
        return dist(a.c, b.c) < 1e-9 and abs(a.r - b.r) < 1e-9 and (a.da > 0) == (b.da > 0)
    return False


def _origin(s: Seg) -> str:
    return s.tag.split(".")[0]


def _merge_collinear(segs: list, across_ops_only: bool = False) -> list:
    def ok(a, b):
        return (_same_carrier(a, b) and dist(a.end(), b.start()) < 1e-9
                and (not across_ops_only or _origin(a) != _origin(b)))
    out: list = []
    for s in segs:
        if out and ok(out[-1], s):
            prev = out.pop()
            out.append(_join(prev, s))
        else:
            out.append(s)
    if len(out) > 1 and ok(out[-1], out[0]):
        first = out.pop(0)
        out[-1] = _join(out[-1], first)
    return out


def _join(a: Seg, b: Seg) -> Seg:
    tag = a.tag if a.tag <= b.tag else b.tag
    ids = frozenset({max(a.ids | b.ids)}) if (a.ids | b.ids) else frozenset()
    if isinstance(a, Line):
        return replace(a, p1=b.p1, tag=tag, ids=ids, v1=b.v1)
    return replace(a, da=a.da + b.da, tag=tag, ids=ids, v1=b.v1)


def from_shapely(geom, src: _Src, log: list) -> Sheet:
    index = _ChordIndex(src)
    faces = []
    polys = list(geom.geoms) if isinstance(geom, (MultiPolygon, shapely.GeometryCollection)) else [geom]
    for p in polys:
        if not isinstance(p, Polygon) or p.is_empty:
            continue
        if p.area < SLIVER_AREA:
            log.append(f"dropped sliver face of area {p.area:.3g} mm^2")
            continue
        p = orient(p, 1.0)
        outer = _recover_ring(list(p.exterior.coords), index, log)
        holes = []
        for h in p.interiors:
            if Polygon(h).area < SLIVER_AREA:
                log.append(f"dropped sliver hole of area {Polygon(h).area:.3g} mm^2")
                continue
            holes.append(_recover_ring(list(h.coords), index, log))
        faces.append(Face(outer, tuple(holes)))
    faces.sort(key=lambda f: (-abs(face_area(f)), f.outer[0].start()))
    return Sheet(tuple(faces))


def boolean(op: str, a: Sheet, tools: Sequence[Sheet], log: list) -> Sheet:
    src = _Src()
    ga = to_shapely(a, src)
    gts = [to_shapely(t, src) for t in tools]
    if op == "union":
        g = shapely.union_all([ga, *gts], grid_size=GRID)
    elif op == "difference":
        g = shapely.difference(ga, shapely.union_all(gts, grid_size=GRID), grid_size=GRID) if gts else ga
    elif op == "intersection":
        g = ga
        for t in gts:
            g = shapely.intersection(g, t, grid_size=GRID)
    else:
        raise GeomError(f"unknown boolean {op}")
    return from_shapely(g, src, log)


def half_plane(p: Pt, direction: Pt, tag: str, size: float = 5000.0) -> Sheet:
    """Sheet = half plane on the LEFT of the directed line through p along direction."""
    L = math.hypot(*direction)
    ux, uy = direction[0] / L, direction[1] / L
    nx, ny = -uy, ux
    a = _pt(p[0] - size * ux, p[1] - size * uy)
    b = _pt(p[0] + size * ux, p[1] + size * uy)
    c = _pt(b[0] + size * nx, b[1] + size * ny)
    d = _pt(a[0] + size * nx, a[1] + size * ny)
    loop = (Line(tag=tag, p0=a, p1=b), Line(tag="farbox", p0=b, p1=c), Line(tag="farbox", p0=c, p1=d),
            Line(tag="farbox", p0=d, p1=a))
    return Sheet((Face(loop),))


# ---------------------------------------------------------------------------
# wires -> faces
def chain_loops(segs: Sequence[Seg], tol: float = 1e-7) -> list:
    """Chain segments into closed loops (orientation fixed later)."""
    pool = list(segs)
    loops = []
    while pool:
        cur = [pool.pop(0)]
        if isinstance(cur[0], Arc) and cur[0].closed:
            loops.append(tuple(cur))
            continue
        guard = 0
        while dist(cur[-1].end(), cur[0].start()) > tol:
            guard += 1
            if guard > 10000:
                raise GeomError("chaining did not terminate")
            e = cur[-1].end()
            hit = None
            for i, s in enumerate(pool):
                if dist(s.start(), e) <= tol:
                    hit = (i, s)
                    break
                if dist(s.end(), e) <= tol:
                    hit = (i, s.reversed())
                    break
            if hit is None:
                raise GeomError(f"open wire: no segment continues from {e} (tags {[s.tag for s in cur]})")
            pool.pop(hit[0])
            cur.append(hit[1])
        loops.append(tuple(_snap_chain(cur)))
    return loops


def _snap_chain(cur: list) -> list:
    """Make consecutive end/start points identical (they agree within tolerance)."""
    out = list(cur)
    for i in range(len(out)):
        a, b = out[i], out[(i + 1) % len(out)]
        e = a.end()
        if isinstance(b, Line) and b.p0 != e:
            out[(i + 1) % len(out)] = replace(b, p0=e)
        elif isinstance(a, Line) and isinstance(b, Arc) and dist(a.p1, b.start()) > 0:
            out[i] = replace(a, p1=b.start())
    return out


def cover(w: Wire) -> Sheet:
    loops = chain_loops(w.segments())
    faces = []
    for lp in loops:
        # AEDT merges consecutive edges on one line/circle when wires are united
        lp = tuple(_merge_collinear(list(lp), across_ops_only=True)) if len(lp) > 1 else lp
        if loop_area(lp) < 0:
            lp = reverse_loop(lp)
        faces.append(Face(lp))
    return Sheet(tuple(faces))


# ---------------------------------------------------------------------------
# vertex / edge lookup
def find_vertex(body, vid: int) -> list:
    """[(loop_index, seg_index)] where vertex = start of seg (end of previous)."""
    out = []
    for li, lp in enumerate(body.loops()):
        n = len(lp)
        closed = isinstance(body, Sheet) or dist(lp[-1].end(), lp[0].start()) < 1e-9
        for i in range(n):
            prev = lp[i - 1] if (i > 0 or closed) else None
            ids = set(lp[i].v0) | (set(prev.v1) if prev is not None else set())
            if vid in ids:
                out.append((li, i))
        if not closed and vid in lp[-1].v1:
            out.append((li, n))
    return out


def find_edge(body, eid: int) -> list:
    return [(li, i) for li, lp in enumerate(body.loops()) for i, s in enumerate(lp) if eid in s.ids]


def replace_loop(body, li: int, new_loop: Loop):
    if isinstance(body, Wire):
        chains = list(body.chains)
        chains[li] = new_loop
        return Wire(tuple(chains))
    k = 0
    faces = []
    for f in body.faces:
        loops = f.loops()
        for j in range(len(loops)):
            if k == li:
                loops[j] = new_loop
            k += 1
        faces.append(Face(loops[0], tuple(loops[1:])))
    return Sheet(tuple(faces))


# ---------------------------------------------------------------------------
# local operations
def _left_normal(t: Pt) -> Pt:
    return _pt(-t[1], t[0])


def _offset_carrier(s: Seg, d: float) -> Seg:
    """Carrier of s shifted by d to the LEFT of its direction."""
    if isinstance(s, Line):
        n = _left_normal(s.tangent_out())
        return replace(s, p0=_pt(s.p0[0] + d * n[0], s.p0[1] + d * n[1]),
                       p1=_pt(s.p1[0] + d * n[0], s.p1[1] + d * n[1]))
    if isinstance(s, Arc):
        r = s.r - d if s.da > 0 else s.r + d
        if r <= 0:
            raise GeomError(f"offset of arc {s.tag} by {d} gives non-positive radius")
        return replace(s, r=r)
    raise GeomError(f"offset of {s.kind} {s.tag} not supported")


def _foot(s: Seg, p: Pt) -> Pt:
    if isinstance(s, Line):
        u = s.project(p)
        return s.at(u)
    if isinstance(s, Arc):
        a = math.atan2(p[1] - s.c[1], p[0] - s.c[0])
        return _pt(s.c[0] + s.r * math.cos(a), s.c[1] + s.r * math.sin(a))
    raise GeomError("foot on curve not supported")


def _trim_end(s: Seg, p: Pt) -> Seg:
    if isinstance(s, Line):
        return replace(s, p1=p)
    if isinstance(s, Arc):
        a = s.angle_of(p)
        return replace(s, da=a - s.a0)
    raise GeomError("trim of curve not supported")


def _trim_start(s: Seg, p: Pt) -> Seg:
    if isinstance(s, Line):
        return replace(s, p0=p)
    if isinstance(s, Arc):
        a = s.angle_of(p)
        return replace(s, a0=a, da=s.a0 + s.da - a)
    raise GeomError("trim of curve not supported")


def fillet_vertex(lp: Loop, i: int, r: float, tag: str) -> tuple:
    """Round the corner at the start of lp[i] with radius r (face on the left)."""
    A, B = lp[i - 1], lp[i]
    V = B.start()
    ta, tb = A.tangent_in(), B.tangent_out()
    turn = ta[0] * tb[1] - ta[1] * tb[0]
    if abs(turn) < 1e-12:
        raise GeomError(f"fillet at a tangent (smooth) vertex {V} between {A.tag} and {B.tag}")
    side = r if turn > 0 else -r  # convex corner: centre inside the face (left); reflex: outside
    oa, ob = _offset_carrier(A, side), _offset_carrier(B, side)
    cands = []
    for C in _all_intersections(oa, ob):
        ta, tb = _foot(A, C), _foot(B, C)
        ua, ub = A.project(ta), B.project(tb)
        if -1e-9 <= ua <= 1 + 1e-9 and -1e-9 <= ub <= 1 + 1e-9:
            cands.append((dist(C, V), C, ta, tb))
    if not cands:
        raise GeomError(f"fillet r={r} at vertex {V} between {A.tag} and {B.tag} does not fit")
    _, C, ta, tb = min(cands, key=lambda c: c[0])
    a0 = math.atan2(ta[1] - C[1], ta[0] - C[0])
    a1 = math.atan2(tb[1] - C[1], tb[0] - C[0])
    da = (a1 - a0 + math.pi) % TWO_PI - math.pi
    arc = Arc(tag=tag, c=C, r=r, a0=a0, da=da)
    newA = _trim_end(A, ta).with_meta(v1=frozenset())
    newB = _trim_start(B, tb).with_meta(v0=frozenset())
    out = list(lp)
    out[i - 1] = newA
    out[i] = newB
    out.insert(i, arc)
    return tuple(out)


def _all_intersections(a: Seg, b: Seg) -> list:
    if isinstance(a, Line) and isinstance(b, Line):
        return _line_line(a.p0, a.p1, b.p0, b.p1)
    if isinstance(a, Line) and isinstance(b, Arc):
        return _line_circle(a.p0, a.p1, b.c, b.r)
    if isinstance(a, Arc) and isinstance(b, Line):
        return _line_circle(b.p0, b.p1, a.c, a.r)
    if isinstance(a, Arc) and isinstance(b, Arc):
        return _circle_circle(a.c, a.r, b.c, b.r)
    raise GeomError("intersection with curve carriers not supported here")


def offset_edge(lp: Loop, i: int, d_out: float, closed: bool = True) -> tuple:
    """Move edge lp[i] by d_out along its outward normal (face on the left), re-trim neighbours."""
    E = lp[i]
    E2 = _offset_carrier(E, -d_out)
    n = len(lp)
    out = list(lp)
    if closed or i > 0:
        A = lp[i - 1]
        P = min(_all_intersections(A, E2), key=lambda p: dist(p, E.start()))
        out[i - 1] = _trim_end(A, P)
        E2 = _trim_start(E2, P)
    if closed or i < n - 1:
        B = lp[(i + 1) % n]
        Q = min(_all_intersections(B, E2), key=lambda p: dist(p, E.end()))
        out[(i + 1) % n] = _trim_start(B, Q)
        E2 = _trim_end(E2, Q)
    out[i] = E2
    return tuple(out)


def body_vertices(body) -> list:
    pts = []
    for lp in body.loops():
        for s in lp:
            pts.append(s.start())
            pts.append(s.end())
    return pts


def sample_loop_xy(lp: Sequence[Seg], max_step_deg: float = 1.0, curve_n: int = 100) -> np.ndarray:
    out = []
    for s in lp:
        if isinstance(s, Line):
            out.append(s.p0)
        elif isinstance(s, Arc):
            n = max(2, int(math.ceil(abs(math.degrees(s.da)) / max_step_deg)))
            out.extend(s.at(k / n) for k in range(n))
        else:
            out.extend(s.at(k / curve_n) for k in range(curve_n))
    out.append(lp[-1].end())
    return np.array(out)


def self_intersections(sh: Sheet) -> list:
    """Return descriptions of invalidities of the sampled sheet (empty if valid)."""
    issues = []
    for k, f in enumerate(sh.faces):
        rings = [sample_loop_xy(lp, 0.5) for lp in f.loops()]
        try:
            poly = Polygon(rings[0], rings[1:])
        except ValueError as e:
            issues.append(f"face {k}: {e}")
            continue
        if not poly.is_valid:
            issues.append(f"face {k}: {shapely.validation.explain_validity(poly)}")
    return issues
