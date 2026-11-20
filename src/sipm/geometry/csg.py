"""Exact Boolean operations on regions bounded by lines and circular arcs.

The rotor laminations are built the way the Maxwell models build them: a sequence of
unite / subtract operations on simple primitives (sector, slots, rib bands, lugs, bridge
bands, holes). :func:`evaluate` computes the boundary of such a constructive expression
exactly:

1. every primitive edge is split at every exact intersection with every other edge
   (line-line, line-circle, circle-circle in closed form);
2. each piece is kept if the expression's membership differs on its two sides (tested
   :data:`PROBE` mm to the left and right of its midpoint) and oriented with the material on
   its left;
3. the pieces are chained into loops (outer loops counter-clockwise, holes clockwise),
   consecutive pieces on the same line or circle are merged, and holes are assigned to the
   outer loop that contains them.

Vertices are therefore exact intersections of exact carriers (to floating point), so a
topology change caused by a parameter (an edge disappearing, two pockets merging) is handled
without special cases.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .primitives import (
    TWO_PI,
    Arc,
    GeometryError,
    Line,
    Pt,
    Segment,
    _circle_circle,
    arc_from,
    chain_area,
)

PROBE = 1e-6       # mm, side offset of the membership probes
SNAP = 1e-9        # mm, intersection points closer than this to an end point are snapped to it
JOIN = 1e-7        # mm, end points closer than this are the same vertex when chaining
MIN_PIECE = 1e-9   # mm, pieces shorter than this are dropped


# ======================================================================================
# expression tree
# ======================================================================================
class Expr:
    def contains(self, p: Pt) -> bool:  # pragma: no cover - interface
        raise NotImplementedError

    def leaves(self) -> list[Leaf]:  # pragma: no cover - interface
        raise NotImplementedError

    def __or__(self, other: Expr) -> Expr:
        return Op("union", self, other)

    def __sub__(self, other: Expr) -> Expr:
        return Op("diff", self, other)

    def __and__(self, other: Expr) -> Expr:
        return Op("inter", self, other)


@dataclass(frozen=True)
class Op(Expr):
    kind: str
    a: Expr
    b: Expr

    def contains(self, p: Pt) -> bool:
        if self.kind == "union":
            return self.a.contains(p) or self.b.contains(p)
        if self.kind == "diff":
            return self.a.contains(p) and not self.b.contains(p)
        return self.a.contains(p) and self.b.contains(p)

    def leaves(self) -> list[Leaf]:
        return self.a.leaves() + self.b.leaves()


class Leaf(Expr):
    """A primitive region bounded by one closed chain of lines and arcs (any orientation)."""

    def __init__(self, chain: Sequence[Segment], name: str = ""):
        for s in chain:
            if s.kind not in ("line", "arc"):
                raise GeometryError(f"Boolean primitive {name!r} has a {s.kind} edge (lines/arcs only)")
        self.chain = tuple(chain)
        self.name = name
        pts = [q for s in self.chain for q in _bbox_points(s)]
        self.bbox = (min(q[0] for q in pts), min(q[1] for q in pts), max(q[0] for q in pts), max(q[1] for q in pts))

    def contains(self, p: Pt) -> bool:
        x0, y0, x1, y1 = self.bbox
        if p[0] < x0 or p[0] > x1 or p[1] < y0 or p[1] > y1:
            return False
        return _crossings(self.chain, p) % 2 == 1

    def leaves(self) -> list[Leaf]:
        return [self]


def union_all(exprs: Sequence[Expr]) -> Expr:
    out = exprs[0]
    for e in exprs[1:]:
        out = out | e
    return out


def _bbox_points(s: Segment) -> list[Pt]:
    if s.kind == "line":
        return [s.p0, s.p1]
    pts = [s.start(), s.end()]
    for k in range(4):
        a = k * math.pi / 2
        if s.contains_angle(a) is not None:
            pts.append((s.c[0] + s.r * math.cos(a), s.c[1] + s.r * math.sin(a)))
    return pts


# ---- point in region (ray to +x, half-open rule) ---------------------------------------
def _crossings(chain: Sequence[Segment], p: Pt) -> int:
    px, py = p
    n = 0
    for s in chain:
        if s.kind == "line":
            (x0, y0), (x1, y1) = s.p0, s.p1
            if (y0 > py) != (y1 > py):
                x = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
                if x > px:
                    n += 1
        else:
            n += _arc_crossings(s, px, py)
    return n


def _arc_crossings(s: Arc, px: float, py: float) -> int:
    """Crossings of the ray y = py, x > px with the arc, counted on y-monotone sub-arcs."""
    lo, hi = (s.a0, s.a0 + s.da) if s.da > 0 else (s.a0 + s.da, s.a0)
    cuts = [lo]
    k = math.ceil((lo - math.pi / 2) / math.pi)
    while True:
        a = math.pi / 2 + k * math.pi
        if a >= hi:
            break
        if a > lo:
            cuts.append(a)
        k += 1
    cuts.append(hi)
    n = 0
    cx, cy, r = s.c[0], s.c[1], s.r
    for a0, a1 in zip(cuts[:-1], cuts[1:]):
        y0, y1 = cy + r * math.sin(a0), cy + r * math.sin(a1)
        if (y0 > py) != (y1 > py):
            am = 0.5 * (a0 + a1)
            dx = math.sqrt(max(r * r - (py - cy) ** 2, 0.0))
            x = cx + (dx if math.cos(am) > 0 else -dx)
            if x > px:
                n += 1
    return n


# ======================================================================================
# splitting
# ======================================================================================
def _seg_param(s: Segment, q: Pt) -> float | None:
    """Parameter of point q on segment s (None if q is not on s within SNAP)."""
    if s.kind == "line":
        dx, dy = s.p1[0] - s.p0[0], s.p1[1] - s.p0[1]
        L2 = dx * dx + dy * dy
        u = ((q[0] - s.p0[0]) * dx + (q[1] - s.p0[1]) * dy) / L2
        L = math.sqrt(L2)
        if u < -SNAP / L or u > 1 + SNAP / L:
            return None
        foot = (s.p0[0] + u * dx, s.p0[1] + u * dy)
        if math.dist(foot, q) > SNAP:
            return None
        return min(max(u, 0.0), 1.0)
    d = math.hypot(q[0] - s.c[0], q[1] - s.c[1])
    if abs(d - s.r) > SNAP:
        return None
    return s.contains_angle(math.atan2(q[1] - s.c[1], q[0] - s.c[0]), SNAP / s.r)


def _intersections(a: Segment, b: Segment) -> list[Pt]:
    """Exact intersection points of two segments plus end points of one lying on the other."""
    pts: list[Pt] = []
    for q in (b.start(), b.end()):
        if _seg_param(a, q) is not None:
            pts.append(q)
    for q in (a.start(), a.end()):
        if _seg_param(b, q) is not None:
            pts.append(q)
    if a.kind == "line" and b.kind == "line":
        dx1, dy1 = a.p1[0] - a.p0[0], a.p1[1] - a.p0[1]
        dx2, dy2 = b.p1[0] - b.p0[0], b.p1[1] - b.p0[1]
        den = dx1 * dy2 - dy1 * dx2
        if abs(den) > 1e-14 * math.hypot(dx1, dy1) * math.hypot(dx2, dy2):
            t = ((b.p0[0] - a.p0[0]) * dy2 - (b.p0[1] - a.p0[1]) * dx2) / den
            q = (a.p0[0] + t * dx1, a.p0[1] + t * dy1)
            pts.append(q)
    elif a.kind == "arc" and b.kind == "arc":
        if math.dist(a.c, b.c) > 1e-12 or abs(a.r - b.r) > 1e-12:
            pts.extend(_circle_circle(a.c, a.r, b.c, b.r))
    else:
        ln, ar = (a, b) if a.kind == "line" else (b, a)
        dx, dy = ln.p1[0] - ln.p0[0], ln.p1[1] - ln.p0[1]
        L = math.hypot(dx, dy)
        fx, fy = ln.p0[0] - ar.c[0], ln.p0[1] - ar.c[1]
        bb = (fx * dx + fy * dy) / L
        cc = fx * fx + fy * fy - ar.r * ar.r
        disc = bb * bb - cc
        if disc >= -1e-12 * ar.r:
            sq = math.sqrt(max(disc, 0.0))
            for s in {-bb - sq, -bb + sq}:
                pts.append((ln.p0[0] + s * dx / L, ln.p0[1] + s * dy / L))
    return [q for q in pts if _seg_param(a, q) is not None and _seg_param(b, q) is not None]


def _split(s: Segment, pts: list[Pt]) -> list[Segment]:
    closed = s.kind == "arc" and abs(abs(s.da) - TWO_PI) < 1e-12
    params = [(0.0, s.start())]
    for q in pts:
        u = _seg_param(s, q)
        if u is not None and math.dist(q, s.start()) > SNAP and math.dist(q, s.end()) > SNAP:
            params.append((u, q))
    params.sort(key=lambda t: t[0])
    # collapse points that coincide (within SNAP); end points stay exact
    clean: list[tuple[float, Pt]] = [params[0]]
    for u, q in params[1:]:
        if math.dist(q, clean[-1][1]) > SNAP:
            clean.append((u, q))
    clean.append((1.0, s.end()))
    if closed and len(clean) == 2:
        return [s]
    out = []
    for (_, q0), (_, q1) in zip(clean[:-1], clean[1:]):
        if math.dist(q0, q1) <= SNAP:
            continue
        piece = Line(q0, q1) if s.kind == "line" else arc_from(s.c, q0, q1, ccw=s.da > 0)
        if s.kind == "arc" and abs(piece.da) > abs(s.da) + 1e-9:  # wrapped the wrong way
            continue
        if piece.length() > MIN_PIECE:
            out.append(piece)
    return out


# ======================================================================================
# evaluation
# ======================================================================================
@dataclass(frozen=True)
class Face:
    """One connected region: counter-clockwise outer chain and clockwise hole chains."""

    outer: tuple
    holes: tuple = ()

    def area(self) -> float:
        return chain_area(self.outer) + sum(chain_area(h) for h in self.holes)

    def loops(self) -> list[tuple]:
        return [self.outer, *self.holes]


def evaluate(expr: Expr, extra_edges: Sequence[Segment] = ()) -> list[Face]:
    """Faces of the region described by ``expr`` (exact lines and arcs)."""
    segs: list[Segment] = [s for lf in expr.leaves() for s in lf.chain] + list(extra_edges)
    boxes = [_seg_box(s) for s in segs]
    cuts: list[list[Pt]] = [[] for _ in segs]
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            if not _boxes_touch(boxes[i], boxes[j]):
                continue
            for q in _intersections(segs[i], segs[j]):
                cuts[i].append(q)
                cuts[j].append(q)
    pieces: list[Segment] = []
    for s, c in zip(segs, cuts):
        pieces.extend(_split(s, c))
    kept: list[Segment] = []
    for pc in pieces:
        m = pc.point(0.5)
        t = pc.tangent(0.5)
        left = (m[0] - PROBE * t[1], m[1] + PROBE * t[0])
        right = (m[0] + PROBE * t[1], m[1] - PROBE * t[0])
        il, ir = expr.contains(left), expr.contains(right)
        if il == ir:
            continue
        kept.append(pc if il else pc.reversed())
    kept = _dedup(kept)
    loops = _chain(kept)
    loops = [_merge(lp) for lp in loops]
    return _faces(loops)


def _seg_box(s: Segment) -> tuple[float, float, float, float]:
    pts = _bbox_points(s)
    return (min(q[0] for q in pts) - SNAP, min(q[1] for q in pts) - SNAP,
            max(q[0] for q in pts) + SNAP, max(q[1] for q in pts) + SNAP)


def _boxes_touch(a, b) -> bool:
    return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])


def _dedup(segs: list[Segment]) -> list[Segment]:
    out: list[Segment] = []
    seen: list[tuple[Pt, Pt, Pt]] = []
    for s in segs:
        key = (s.start(), s.end(), s.point(0.5))
        if any(math.dist(key[0], k[0]) < JOIN and math.dist(key[1], k[1]) < JOIN and math.dist(key[2], k[2]) < JOIN
               for k in seen):
            continue
        seen.append(key)
        out.append(s)
    return out


def _chain(segs: list[Segment]) -> list[list[Segment]]:
    remaining = list(range(len(segs)))
    starts = [s.start() for s in segs]
    used = [False] * len(segs)
    loops: list[list[Segment]] = []
    for i0 in remaining:
        if used[i0]:
            continue
        loop = [segs[i0]]
        used[i0] = True
        first = segs[i0].start()
        cur = segs[i0]
        for _ in range(len(segs) + 1):
            end = cur.end()
            if math.dist(end, first) < JOIN and len(loop) > 0:
                break
            cands = [k for k in range(len(segs)) if not used[k] and math.dist(starts[k], end) < JOIN]
            if not cands:
                raise GeometryError(f"boundary does not close at ({end[0]:.6f}, {end[1]:.6f})")
            if len(cands) > 1:
                cands.sort(key=lambda k: _cw_turn(cur, segs[k]))
            k = cands[0]
            used[k] = True
            cur = segs[k]
            loop.append(cur)
        loops.append(loop)
    return loops


def _cw_turn(a: Segment, b: Segment) -> float:
    """Clockwise angle from the reversed incoming direction to the outgoing direction."""
    ti, to = a.tangent(1.0), b.tangent(0.0)
    back = math.atan2(-ti[1], -ti[0])
    out = math.atan2(to[1], to[0])
    return (back - out) % TWO_PI


def _same_carrier(a: Segment, b: Segment) -> bool:
    if a.kind != b.kind:
        return False
    if a.kind == "line":
        ta, tb = a.tangent(), b.tangent()
        if abs(ta[0] * tb[1] - ta[1] * tb[0]) > 1e-12 or ta[0] * tb[0] + ta[1] * tb[1] < 0:
            return False
        dx, dy = b.p1[0] - a.p0[0], b.p1[1] - a.p0[1]
        return abs(ta[0] * dy - ta[1] * dx) < 1e-9
    return math.dist(a.c, b.c) < 1e-9 and abs(a.r - b.r) < 1e-9 and (a.da > 0) == (b.da > 0)


def _join(a: Segment, b: Segment) -> Segment:
    if a.kind == "line":
        return Line(a.p0, b.p1)
    da = a.da + b.da
    if abs(abs(da) - TWO_PI) < 1e-9:
        return Arc(a.c, a.r, a.a0, math.copysign(TWO_PI, da))
    return Arc(a.c, a.r, a.a0, da)


def _merge(loop: list[Segment]) -> tuple:
    out = list(loop)
    changed = True
    while changed and len(out) > 1:
        changed = False
        for i in range(len(out)):
            a, b = out[i], out[(i + 1) % len(out)]
            if a is b:
                continue
            if _same_carrier(a, b):
                j = (i + 1) % len(out)
                out[i] = _join(a, b)
                del out[j]
                changed = True
                break
    if len(out) == 1 and out[0].kind == "arc" and abs(abs(out[0].da) - TWO_PI) < 1e-9:
        a = out[0]
        out = [Arc(a.c, a.r, 0.0, a.da)]
    return tuple(out)


def _faces(loops: list[tuple]) -> list[Face]:
    outers = [lp for lp in loops if chain_area(lp) > 0]
    holes = [lp for lp in loops if chain_area(lp) <= 0]
    assigned: dict[int, list[tuple]] = {i: [] for i in range(len(outers))}
    for h in holes:
        q = h[0].point(0.5)
        t = h[0].tangent(0.5)
        probe = (q[0] - 10 * PROBE * t[1], q[1] + 10 * PROBE * t[0])  # material side of the hole edge
        owners = [i for i, o in enumerate(outers) if _crossings(o, probe) % 2 == 1]
        if not owners:
            raise GeometryError("a hole boundary lies outside every face")
        best = min(owners, key=lambda i: chain_area(outers[i]))
        assigned[best].append(h)
    faces = [Face(o, tuple(assigned[i])) for i, o in enumerate(outers)]
    return sorted(faces, key=lambda f: -f.area())
