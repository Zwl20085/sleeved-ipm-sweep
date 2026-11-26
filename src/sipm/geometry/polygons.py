"""Sampled-polygon helpers (distances, point-in-polygon, self-intersection) and the Loop fillet."""

from __future__ import annotations

import math

import numpy as np

from .primitives import GeometryError, Loop


def fillet(loop: Loop, idx: int, r: float) -> Loop:
    """Round the corner at vertex ``idx`` (both adjacent segments must be lines) with
    radius ``r``. The arc is placed inside the angle (< 180 deg) formed at the vertex."""
    n = len(loop)
    if r <= 0:
        return loop
    if loop.arcs[idx - 1 if idx > 0 else n - 1] is not None or loop.arcs[idx] is not None:
        raise GeometryError("fillet requested at a corner adjacent to an arc")
    v = loop.pts[idx]
    pprev = loop.pts[(idx - 1) % n]
    pnext = loop.pts[(idx + 1) % n]
    d1 = (pprev - v) / np.linalg.norm(pprev - v)
    d2 = (pnext - v) / np.linalg.norm(pnext - v)
    cosang = float(np.clip(d1 @ d2, -1.0, 1.0))
    ang = math.acos(cosang)
    if ang < 1e-6 or abs(ang - math.pi) < 1e-6:
        return loop
    half = ang / 2
    t = r / math.tan(half)
    l1 = np.linalg.norm(pprev - v)
    l2 = np.linalg.norm(pnext - v)
    # fillets are applied one after another, so l1/l2 are what earlier fillets left over
    if t > l1 - 0.02 or t > l2 - 0.02:
        raise GeometryError(
            f"fillet radius {r} mm does not fit at corner ({v[0]:.3f}, {v[1]:.3f}) mm "
            f"(needs {t:.3f} mm of straight edge, has {min(l1, l2):.3f} mm)")
    t1 = v + t * d1
    t2 = v + t * d2
    bis = (d1 + d2) / np.linalg.norm(d1 + d2)
    c = v + (r / math.sin(half)) * bis
    pts = list(loop.pts)
    arcs = list(loop.arcs)
    pts[idx:idx + 1] = [t1, t2]
    arcs[idx:idx] = [c]  # segment t1->t2 is the arc; segment t2->next keeps old arcs[idx] (line)
    return Loop(np.array(pts), arcs)


# ======================================================================================
# helpers
# ======================================================================================
def _circle_line_param(p0, d, radius):
    """Largest s with |p0 + s d| = radius (d unit vector)."""
    b = p0 @ d
    c = p0 @ p0 - radius * radius
    disc = b * b - c
    if disc < 0:
        raise GeometryError("construction line misses the target circle")
    return -b + math.sqrt(disc)


def _pol(r, th):
    return np.array([r * math.cos(th), r * math.sin(th)])


def seg_seg_distance(p1, p2, q1, q2) -> float:
    def pt_seg(p, a, b):
        ab = b - a
        L = ab @ ab
        t = 0.0 if L == 0 else float(np.clip((p - a) @ ab / L, 0, 1))
        return np.linalg.norm(p - (a + t * ab))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    d1, d2 = cross(q1, q2, p1), cross(q1, q2, p2)
    d3, d4 = cross(p1, p2, q1), cross(p1, p2, q2)
    if (d1 * d2 < 0) and (d3 * d4 < 0):
        return 0.0
    return min(pt_seg(p1, q1, q2), pt_seg(p2, q1, q2), pt_seg(q1, p1, p2), pt_seg(q2, p1, p2))


def _pts_to_segs(P, A, B):
    """Distance matrix (len(P), len(A)) from points P to segments A->B."""
    AB = B - A
    L = np.einsum("ij,ij->i", AB, AB)
    L = np.where(L == 0, 1.0, L)
    AP = P[:, None, :] - A[None, :, :]
    t = np.clip(np.einsum("pij,ij->pi", AP, AB) / L[None, :], 0.0, 1.0)
    proj = A[None, :, :] + t[..., None] * AB[None, :, :]
    return np.linalg.norm(P[:, None, :] - proj, axis=2)


def _segs_cross(A, B, C, D):
    """Boolean matrix: proper crossing of segments A->B (rows) with C->D (cols)."""
    def cross(o, a, b):
        return (a[..., 0] - o[..., 0]) * (b[..., 1] - o[..., 1]) - (a[..., 1] - o[..., 1]) * (b[..., 0] - o[..., 0])

    A_, B_ = A[:, None, :], B[:, None, :]
    C_, D_ = C[None, :, :], D[None, :, :]
    d1, d2 = cross(C_, D_, A_), cross(C_, D_, B_)
    d3, d4 = cross(A_, B_, C_), cross(A_, B_, D_)
    return (d1 * d2 < 0) & (d3 * d4 < 0)


def polyline_distance(P, Q) -> float:
    """Minimum distance between two closed polylines (0 if they cross)."""
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    Pa, Pb = P, np.roll(P, -1, axis=0)
    Qa, Qb = Q, np.roll(Q, -1, axis=0)
    if _segs_cross(Pa, Pb, Qa, Qb).any():
        return 0.0
    return float(min(_pts_to_segs(P, Qa, Qb).min(), _pts_to_segs(Q, Pa, Pb).min()))


def points_in_polygon(points, poly, tol=1e-9) -> np.ndarray:
    """Even-odd test; points within ``tol`` of an edge count as inside."""
    from matplotlib.path import Path

    points = np.asarray(points, float)
    poly = np.asarray(poly, float)
    path = Path(np.vstack([poly, poly[:1]]), closed=True)
    inside = path.contains_points(points, radius=0.0)
    if tol > 0:
        d = _pts_to_segs(points, poly, np.roll(poly, -1, axis=0)).min(axis=1)
        inside |= d <= tol
    return inside


def self_intersects(poly) -> bool:
    poly = np.asarray(poly, float)
    n = len(poly)
    A, B = poly, np.roll(poly, -1, axis=0)
    X = _segs_cross(A, B, A, B)
    idx = np.arange(n)
    near = (np.abs(idx[:, None] - idx[None, :]) <= 1) | (np.abs(idx[:, None] - idx[None, :]) == n - 1)
    return bool((X & ~near).any())
