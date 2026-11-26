import itertools
import math

import numpy as np
import pytest

from sipm.geometry import (
    BIARC_TOL,
    GeometryError,
    Hermite,
    Line,
    Loop,
    build_rotor,
    build_stator,
    chain_area,
    csg,
    fillet,
    points_in_polygon,
    winding_layout,
)
from sipm.geometry.checks import in_steel
from sipm.geometry.primitives import (
    annular_sector,
    chain_mirror_x,
    chain_polyline,
    chain_to_loop,
    chain_transformed,
    circle_chain,
    hermite_max_deviation,
    hermite_to_biarcs,
    polygon_chain,
)
from sipm.params import TOPOLOGIES, Params, default_params


def _pol(r, a):
    return np.array([r * math.cos(a), r * math.sin(a)])


# ======================================================================================
# primitives
# ======================================================================================
def test_loop_area_exact_for_arcs():
    r, a = 10.0, math.pi / 4
    sector = Loop(np.array([[0, 0], _pol(r, 0), _pol(r, a)]), [None, (0, 0), None])
    assert math.isclose(sector.area(), 0.5 * r * r * a, rel_tol=1e-12)
    ann = Loop(np.array([_pol(5, 0), _pol(r, 0), _pol(r, a), _pol(5, a)]), [None, (0, 0), None, (0, 0)])
    assert math.isclose(ann.area(), 0.5 * (r * r - 25) * a, rel_tol=1e-12)
    P = ann.discretize(0.05)
    x, y = P[:, 0], P[:, 1]
    assert math.isclose(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))), ann.area(), rel_tol=1e-6)


def test_loop_mirror_keeps_area_and_orientation():
    lp = Loop(np.array([_pol(5, 0.1), _pol(9, 0.1), _pol(9, 0.5), _pol(5, 0.5)]), [None, (0, 0), None, (0, 0)]).ccw()
    m = lp.transformed(mirror_x=True)
    assert math.isclose(m.area(), lp.area(), rel_tol=1e-12)
    assert m.signed_area() > 0


def test_fillet_tangent_and_area():
    sq = Loop(np.array([[0, 0], [4, 0], [4, 4], [0, 4]]))
    f = fillet(sq, 2, 1.0)
    assert math.isclose(f.area(), 16 - (1 - math.pi / 4), rel_tol=1e-12)
    i = next(k for k, c in enumerate(f.arcs) if c is not None)
    assert np.allclose(f.arcs[i], [3, 3])
    assert np.allclose(f.pts[i], [4, 3]) and np.allclose(f.pts[i + 1], [3, 4])
    with pytest.raises(GeometryError):
        fillet(sq, 2, 3.99)


def test_exact_chain_area_and_solver_loop_agree():
    ch = annular_sector(30.0, 70.0, -0.3, 0.4)
    assert math.isclose(chain_area(ch), 0.5 * (70 ** 2 - 30 ** 2) * 0.7, rel_tol=1e-14)
    assert math.isclose(chain_to_loop(ch).area(), chain_area(ch), rel_tol=1e-12)
    circ = circle_chain((5.0, 1.0), 2.0)
    assert math.isclose(chain_to_loop(circ).area(), 4 * math.pi, rel_tol=1e-12)  # split into minor arcs


def test_hermite_area_is_exact_and_biarcs_meet_tolerance():
    h = Hermite((0.0, 0.0), (3.0, 0.0), (2.0, 1.0), (0.0, 4.0))
    ch = (h, Line((2.0, 1.0), (0.0, 0.0)))
    P = chain_polyline(ch, 1.0, 20000)
    x, y = P[:, 0], P[:, 1]
    poly = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    assert math.isclose(chain_area(ch), poly, rel_tol=1e-8)
    segs = hermite_to_biarcs(h)
    assert math.dist(segs[0].start(), h.start()) < 1e-12 and math.dist(segs[-1].end(), h.end()) < 1e-12
    assert hermite_max_deviation(h, segs) <= BIARC_TOL
    for a, b in zip(segs[:-1], segs[1:]):  # tangent-continuous
        ta, tb = a.tangent(1.0), b.tangent(0.0)
        assert abs(ta[0] * tb[1] - ta[1] * tb[0]) < 1e-9


# ======================================================================================
# exact Booleans
# ======================================================================================
def test_csg_union_difference_and_hole():
    a = csg.Leaf(polygon_chain([(0, 0), (2, 0), (2, 2), (0, 2)]))
    b = csg.Leaf(polygon_chain([(1, 1), (3, 1), (3, 3), (1, 3)]))
    u = csg.evaluate(a | b)
    assert len(u) == 1 and math.isclose(u[0].area(), 7.0, rel_tol=1e-14)
    assert len(u[0].outer) == 8
    d = csg.evaluate(a - b)
    assert math.isclose(d[0].area(), 3.0, rel_tol=1e-14)
    h = csg.evaluate(a - csg.Leaf(circle_chain((1.0, 1.0), 0.5)))
    assert len(h) == 1 and len(h[0].holes) == 1
    assert math.isclose(h[0].area(), 4 - math.pi / 4, rel_tol=1e-13)
    two = csg.evaluate(a - csg.Leaf(polygon_chain([(0.9, -1), (1.1, -1), (1.1, 3), (0.9, 3)])))
    assert len(two) == 2 and math.isclose(sum(f.area() for f in two), 3.6, rel_tol=1e-13)


def test_csg_arc_line_intersection_exact():
    disc = csg.Leaf(circle_chain((0.0, 0.0), 1.0))
    half = csg.Leaf(polygon_chain([(-2, 0), (2, 0), (2, 2), (-2, 2)]))
    f = csg.evaluate(disc & half)
    assert math.isclose(f[0].area(), math.pi / 2, rel_tol=1e-14)
    pts = {tuple(round(v, 12) for v in s.start()) for s in f[0].outer}
    assert pts == {(1.0, 0.0), (-1.0, 0.0)}


# ======================================================================================
# rotor
# ======================================================================================
def _faces_mirror_equal(faces, M) -> float:
    """Worst vertex deviation between the faces and their image under M (vertex sets)."""
    a = np.array([s.start() for f in faces for lp in [f.outer, *f.holes] for s in lp])
    b = np.array([s.start() for f in faces for lp in [f.outer, *f.holes] for s in chain_transformed(lp, M)])
    d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min(1)
    return float(d.max())


@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_rotor_symmetric_about_d_and_q_axes(topo):
    """The mechanical model uses half a pole with symmetry on both cut lines; this needs the
    pole to be symmetric about the d-axis (y = 0) and the q-axis (theta = pi/p)."""
    rp = build_rotor(default_params(topo))
    assert _faces_mirror_equal(rp.lamination, np.diag([1.0, -1.0])) < 1e-9
    # q-axis: membership of random points near the q-axis equals that of their mirror image
    # about the q-axis, mapped back into this pole by the pole pitch (periodicity)
    hp = rp.half_pitch
    rng = np.random.default_rng(1)
    r = rng.uniform(rp.r_shaft, rp.r_lam, 4000)
    th = rng.uniform(hp - 0.25, hp, 4000)
    pts = np.c_[r * np.cos(th), r * np.sin(th)]
    back = np.c_[r * np.cos(hp - th - hp), r * np.sin(hp - th - hp)]  # angle 2 hp - th, minus 2 hp
    assert (in_steel(rp, pts) == in_steel(rp, back)).all()
    for mg in rp.magnets:
        twin = [m for m in rp.magnets if m.layer == mg.layer and m.side == -mg.side][0]
        assert math.isclose(mg.area(), twin.area(), rel_tol=1e-12)


@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_rotor_regions_consistent(topo):
    p = default_params(topo)
    rp = build_rotor(p)
    m = rp.metrics
    # every magnet in air (not overlapping steel), magnets disjoint
    for mg in rp.magnets:
        assert not in_steel(rp, mg.centre[None, :]).any()
    for i, a in enumerate(rp.magnets):
        for b in rp.magnets[i + 1:]:
            assert not points_in_polygon(a.loop.discretize(5), b.loop.discretize(5), tol=0).any()
    # magnet area = (width - side clearance) x (length - end clearance) - corner rounding
    rc = p.pm_corner_radius
    for mg in rp.magnets:
        L = rp.layers[mg.layer - 1]
        expect = (L.inp.W - p.pm_side_clearance) * (L.length - p.pm_end_clearance) - (4 - math.pi) * rc * rc
        assert math.isclose(mg.area(), expect, rel_tol=1e-12)
    # lamination + air inside the sector = sector area (holes are inside the faces)
    assert 0 < m["lamination_area"] < chain_area(annular_sector(rp.r_shaft, rp.r_lam, -rp.half_pitch,
                                                                  rp.half_pitch))
    # half pole is half of the full pole
    assert math.isclose(sum(f.area() for f in rp.half_lamination()), m["lamination_area"] / 2, rel_tol=2e-5)
    if topo != "sleeved_bridgeless":
        assert math.isclose(m["pm1_rib_min"], 2 * p.pm1_rib_half_width, rel_tol=1e-6)
        assert m["pm2_bridge_min"] > 0 and m["pm1_bridge_min"] > 0
    assert m["lamination_pieces"] == (3 if topo == "sleeved_bridgeless" else 1)


def test_bridged_bridges_equal_parameters():
    p = default_params("sleeved_bridged")
    m = build_rotor(p).metrics
    assert math.isclose(m["pm1_bridge_min"], p.pm1_bridge_thickness, abs_tol=1e-9)
    assert math.isclose(m["pm2_bridge_min"], p.pm2_bridge_thickness, abs_tol=1e-9)
    assert math.isclose(m["pm1_rib_min"], 2 * p.pm1_rib_half_width, abs_tol=1e-6)


@pytest.mark.parametrize("t", [0.5, 1.0, 1.5, 2.0])
def test_airgap_constant_across_sleeve_thickness(t):
    for topo in ("sleeved_bridged", "sleeved_bridgeless"):
        p = default_params(topo, sleeve_thickness=t)
        rp, sp = build_rotor(p), build_stator(p)
        r_sl = np.hypot(*rp.sleeve.pts.T).max()
        r_bore = np.hypot(*sp.loop.pts.T).min()
        assert math.isclose(r_bore - r_sl, 1.32, abs_tol=1e-9)
        assert math.isclose(np.hypot(*rp.sector.pts.T).max(), r_sl - t, abs_tol=1e-9)


def test_thicker_sleeve_keeps_layers_and_shortens_magnets():
    """Maxwell construction: the layers stay put when the sleeve thickens; the bridges stay
    the same thickness below the new OD, so the free slot length and the magnet shrink."""
    rps = [build_rotor(default_params("sleeved_bridged", sleeve_thickness=t)) for t in (0.5, 1.0, 1.5, 2.0)]
    for k in (0, 1):
        P8 = [rp.layers[k].p("P8") for rp in rps]
        assert max(math.dist(P8[0], q) for q in P8) < 1e-12
        lengths = [rp.layers[k].length for rp in rps]
        assert all(a > b for a, b in zip(lengths, lengths[1:]))


def test_removed_conventional_topology_gives_a_clear_error():
    for make in (lambda: default_params("conventional"), lambda: Params(topology="conventional")):
        with pytest.raises(ValueError, match="conventional.*removed from the package"):
            make()


def test_bridgeless_pieces_and_openings():
    p = default_params("sleeved_bridgeless")
    rp = build_rotor(p)
    assert len(rp.lamination) == 3
    od = [s for f in rp.lamination for s in f.outer if s.kind == "arc" and abs(s.r - rp.r_lam) < 1e-9]
    covered = sum(abs(s.da) for s in od)
    opened = 2 * math.radians(p.pm1_opening_angle + p.pm2_opening_angle)
    assert math.isclose(covered, 2 * rp.half_pitch - opened, rel_tol=1e-9)


# ---- feasibility checks: one per constraint --------------------------------------------
@pytest.mark.parametrize("topo, bad, match", [
    ("sleeved_bridged", dict(pm1_length_ratio=0.1), "shorter than"),                          # too-short magnet
    ("sleeved_bridged", dict(pm1_bridge_thickness=14.0), "does not reach the bridge circle"),  # negative magnet
    ("sleeved_bridged", dict(pm2_d_crossing=58.0, pm1_d_crossing=60.0), "overlap|merge|steel between"),
    ("sleeved_bridged", dict(pm2_angle=140.0), "beyond the q-axis|q-axis"),                    # band vs slot
    ("sleeved_bridged", dict(pm1_d_crossing=40.0), "layer order"),                              # layer order
    ("sleeved_bridged", dict(d_hole_depth=10.0), "d-axis hole"),                                # pocket vs hole
    ("sleeved_bridged", dict(pm2_web_width=0.2, pm2_slot_start=0.0), "cavity"),                 # pocket vs cavity
    ("sleeved_bridged", dict(shaft_diameter=94.0), "bore"),                                     # pocket vs shaft
    ("sleeved_bridged", dict(shaft_diameter=84.0), "lightening holes come within 0.000 mm of the bore"),
    ("sleeved_bridgeless", dict(shaft_diameter=83.5), "lightening holes come within 0.250 mm of the bore"),
    ("sleeved_bridged", dict(pm2_slot_start=9.0), "overlaps the lamination"),                   # magnet vs steel
    ("sleeved_bridgeless", dict(pm1_opening_shift=3.0), "not inside the slot mouth"),           # lips vs openings
    ("sleeved_bridgeless", dict(pm1_opening_shift=-1.5), "not inside the slot mouth"),          # ... q side
    ("sleeved_bridgeless", dict(pm2_opening_angle=4.5), "not inside the slot mouth"),           # ... too wide
    ("sleeved_bridged", dict(min_q_web=10.0), "q-axis"),                                        # q-axis web
])
def test_infeasible_geometry_rejected(topo, bad, match):
    with pytest.raises(GeometryError, match=match):
        build_rotor(default_params(topo, **bad))


SWEEP_CORNERS = {
    "sleeved_bridged": dict(sleeve_thickness=(0.5, 2.0), pm1_bridge_thickness=(0.5, 1.5),
                            pm2_bridge_thickness=(0.5, 1.5), pm1_rib_half_width=(0.5, 1.5), pm2_web_width=(0.5, 1.5)),
    "sleeved_bridgeless": dict(sleeve_thickness=(0.5, 2.0), lip_depth=(0.5, 2.0), pm1_bridge_thickness=(0.5, 1.5),
                               pm2_bridge_thickness=(0.5, 1.5)),
}
# Corners of the original sweep ranges that are infeasible, with the constraint they violate
# (see docs/methods.md, "Feasible sweep ranges").
INFEASIBLE_CORNERS: dict = {}


@pytest.mark.parametrize("topo", list(SWEEP_CORNERS))
def test_every_sweep_corner_builds_or_names_the_violated_constraint(topo):
    names = list(SWEEP_CORNERS[topo])
    failed = {}
    for combo in itertools.product(*SWEEP_CORNERS[topo].values()):
        kw = dict(zip(names, combo))
        try:
            build_rotor(default_params(topo, **kw))
        except GeometryError as e:
            assert str(e)
            failed[tuple(sorted(kw.items()))] = str(e)
    assert set(failed) == set(INFEASIBLE_CORNERS.get(topo, {})), failed


# ======================================================================================
# stator
# ======================================================================================
def test_infeasible_stator_rejected():
    with pytest.raises((GeometryError, ValueError)):
        build_stator(default_params("sleeved_bridged", slot_width=9.0))
    with pytest.raises(GeometryError):
        build_stator(default_params("sleeved_bridged", cond_height=2.0))


def test_stator_conductors_inside_slots():
    p = default_params("sleeved_bridged")
    sp = build_stator(p)
    assert len(sp.conductors) == p.slots // p.poles
    stator_poly = sp.loop.discretize(0.2)
    for c in sp.conductors:
        assert not points_in_polygon(c.loop.pts, stator_poly, tol=0).any()
        assert math.isclose(c.loop.area(), p.cond_width * p.cond_height * p.cond_count, rel_tol=1e-9)
        assert c.turns == p.hairpin_layers
    ri, ro = p.stator_id / 2, p.stator_od / 2
    full = 0.5 * (ro ** 2 - ri ** 2) * 2 * math.pi / p.poles
    assert sp.loop.area() < full


def test_winding_balanced_three_phase():
    p = default_params("sleeved_bridged")
    layout, axes = winding_layout(p, build_stator(p).slot_angles)
    assert [ph for ph, _ in layout] == ["A", "A", "C", "C", "B", "B"]
    assert math.isclose((axes["B"] - axes["A"]) % 360, 120, abs_tol=1e-9)
    assert math.isclose((axes["C"] - axes["A"]) % 360, 240, abs_tol=1e-9)


def test_lower_magnets_are_mirror_images():
    rp = build_rotor(default_params("sleeved_bridged"))
    for k in (1, 2):
        up = [m for m in rp.magnets if m.layer == k and m.side > 0][0]
        lo = [m for m in rp.magnets if m.layer == k and m.side < 0][0]
        a = sorted(tuple(np.round(s.start(), 12)) for s in chain_mirror_x(up.chain))
        b = sorted(tuple(np.round(s.start(), 12)) for s in lo.chain)
        assert a == b
        assert np.allclose(lo.magnetisation, [up.magnetisation[0], -up.magnetisation[1]])


def test_cavity_curves_as_biarcs_within_tolerance():
    """Both solvers get the Hermite curves of the central cavity as the same tangent-continuous
    biarcs, at most BIARC_TOL from the true curve."""
    rp = build_rotor(default_params("sleeved_bridged"))
    curves = [s for f in rp.lamination for h in f.holes for s in h if s.kind == "hermite"]
    assert len(curves) == 4  # two per half of the cavity
    worst = max(hermite_max_deviation(h, hermite_to_biarcs(h), n=400) for h in curves)
    assert worst <= BIARC_TOL
