"""Layered sleeve and winding prestress: parameters, membrane homogenisation, the disc + N rings
closed form (with hoop eigenstrain) and the sequential-winding conversion, deck text, per-layer
failure and interface outputs, CLI and example configuration."""

import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from sipm import analytic as an
from sipm import composite as c
from sipm.cli import main, sweep_from_config
from sipm.geometry import build_rotor
from sipm.mapdl import (
    MAT_LAYER0,
    mech_model_from_rotor,
    mech_model_ring,
    parse_sleeve_nodes,
    sleeve_failure_outputs,
    sleeve_interface_outputs,
    write_deck,
)
from sipm.maxwell import MECH_ONLY, em_model_key
from sipm.params import (
    MAX_SLEEVE_LAYERS,
    coerce,
    default_params,
    layer_radii,
    parse_layer_prestress,
    parse_sleeve_layers,
)
from sipm.sweep import KEY_OMITTED_AT_DEFAULT, point_key
from sipm.validation import MULTILAYER_CASE, disc_body, layer_bodies, multilayer_ring_reference, winding_prestress

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

FIX = Path(__file__).parent / "fixtures"
SPEC = dict(sleeve_fibre="carbon_sm")
NEW_FIELDS = ("sleeve_layers", "sleeve_layer_prestress")
A, B = 71.1, 72.6
W = 20000 * 2 * math.pi / 60


# ---- parsing and validation -----------------------------------------------------------------
def test_parse_sleeve_layers():
    assert parse_sleeve_layers("") == () and parse_sleeve_layers("  ") == ()
    assert parse_sleeve_layers("0:2, 20:1 ,0:1") == ((0.0, 2.0), (20.0, 1.0), (0.0, 1.0))
    assert parse_sleeve_layers("90:0.5") == ((90.0, 0.5),)


@pytest.mark.parametrize("bad, match", [
    ("0", "angle:weight"), ("0:1:2", "angle:weight"), ("a:1", "not a number"), ("0:nan", "not finite"),
    ("-5:1", r"\[0, 90\]"), ("95:1", r"\[0, 90\]"), ("0:0", "> 0"), ("0:-1", "> 0"), ("0:1,", "angle:weight"),
    (",".join(["0:1"] * (MAX_SLEEVE_LAYERS + 1)), "at most"),
])
def test_parse_sleeve_layers_rejects(bad, match):
    with pytest.raises(ValueError, match=match):
        parse_sleeve_layers(bad)


def test_parse_layer_prestress():
    assert parse_layer_prestress("") == () and parse_layer_prestress("100") == (100.0,)
    assert parse_layer_prestress("100, 0,50.5") == (100.0, 0.0, 50.5)
    for bad in ("x", "100,", "-1", "inf"):
        with pytest.raises(ValueError):
            parse_layer_prestress(bad)


def test_coerce_arrays_and_numbers():
    assert coerce("sleeve_layers", ["0:2", "20:1"]) == "0:2,20:1"
    assert coerce("sleeve_layer_prestress", [100, 50.5]) == "100,50.5"
    assert coerce("sleeve_layer_prestress", 80) == "80"


def test_defaults_single_layer_unchanged():
    p = default_params("sleeved_bridged")
    assert p.sleeve_layers == "" and p.sleeve_layer_prestress == "" and not p.is_layered
    assert p.sleeve_layer_count == 1 and p.sleeve_layer_fractions() == (1.0,) and p.sleeve_layer_angles() == (0.0,)
    assert p.sleeve_layer_prestress_mpa() == (0.0,) and p.sleeve_layer_prestrain() == (0.0,)
    assert p.sleeve_layer_radii() == (p.r_lam, p.r_sleeve_out)


@pytest.mark.parametrize("bad, match", [
    (dict(sleeve_layers="0:1,20:1"), "needs a composite spec"),
    (dict(**SPEC, sleeve_layers="0:1,20:1", sleeve_wind_angle=10.0), "sleeve_wind_angle must be 0"),
    (dict(**SPEC, sleeve_layers="0:1,20:1,0:1", sleeve_layer_prestress="1,2"), "3 sleeve layer"),
    (dict(sleeve_layer_prestress="1,2"), "1 sleeve layer"),
    (dict(sleeve_layer_prestress="-5"), ">= 0"),
    (dict(**SPEC, sleeve_layers="0:30,20:1"), "thinnest layer"),          # 1.5 mm / 31 = 0.048 mm
    (dict(**SPEC, sleeve_layers="0:1,abc:1"), "not a number"),
    (dict(**SPEC, sleeve_layers=",".join(["0:1"] * 11)), "at most"),
])
def test_invalid_layer_specs_rejected(bad, match):
    with pytest.raises(ValueError, match=match):
        default_params("sleeved_bridged", **bad)


def test_thinnest_layer_limit_follows_thickness():
    default_params("sleeved_bridged", **SPEC, sleeve_layers="0:1,20:1,0:1,20:1,0:1,20:1,0:1,20:1,0:1,20:1")
    default_params("sleeved_bridgeless", **SPEC, sleeve_layers=",".join(["0:1"] * 10))   # 0.5 mm / 10 = 0.05 mm
    with pytest.raises(ValueError, match="thinnest"):
        default_params("sleeved_bridgeless", **SPEC, sleeve_layers="0:10,0:1")           # 0.5 mm / 11


def test_layer_helpers():
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1", sleeve_layer_prestress="100,40,100")
    assert p.is_layered and p.sleeve_layer_count == 3
    assert p.sleeve_layer_fractions() == (0.5, 0.25, 0.25) and p.sleeve_layer_angles() == (0.0, 20.0, 0.0)
    r = p.sleeve_layer_radii()
    assert r[0] == p.r_lam and r[-1] == p.r_sleeve_out and np.allclose(np.diff(r), [0.75, 0.375, 0.375])
    lams = p.sleeve_layer_laminates()
    assert [lam.angle_deg for lam in lams] == [0.0, 20.0, 0.0]
    eps = p.sleeve_layer_prestrain()
    expected = [s / (lam.E_hoop * 1e3) for s, lam in zip((100, 40, 100), lams)]
    assert eps == pytest.approx(expected, rel=1e-14)
    assert p.replace(sleeve_layer_prestress="60").sleeve_layer_prestress_mpa() == (60.0, 60.0, 60.0)
    # direct constants, one layer: prestrain from sleeve_E_hoop
    q = default_params("sleeved_bridged", sleeve_layer_prestress="125")
    assert q.sleeve_layer_prestrain() == (125 / 125e3,)
    assert layer_radii(10.0, 12.0, (0.5, 0.25, 0.25)) == (10.0, 11.0, 11.5, 12.0)


# ---- the thirteen constants of a layered sleeve ------------------------------------------------
def test_stack_constants_are_membrane_homogenisation():
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1")
    stack = c.layer_stack(p.sleeve_layer_laminates(), (2, 1, 1))
    for k, v in c.sleeve_params(stack).items():
        assert getattr(p, k) == pytest.approx(v, rel=1e-14, abs=1e-20), k
    assert isinstance(p.sleeve_laminate(), c.LayerStack)
    # independent check of the membrane rule: the in-plane stiffness (hoop, axial, hoop-axial) at zero
    # out-of-plane stress is the thickness average of the layers' plane-stress stiffnesses
    def q_inplane(S):
        idx = [0, 1, 5]
        return np.linalg.inv(S[np.ix_(idx, idx)])
    avg = sum(f * q_inplane(lam.S) for f, lam in zip(p.sleeve_layer_fractions(), p.sleeve_layer_laminates()))
    assert np.allclose(q_inplane(p.sleeve_compliance()), avg, rtol=1e-12)
    # hoop stiffness between the layers' and the density the common one
    lams = p.sleeve_layer_laminates()
    Eh = [lam.E_hoop for lam in lams]
    assert min(Eh) < p.sleeve_E_hoop < max(Eh) and p.sleeve_density == pytest.approx(lams[0].density)
    # free thermal strain: at the stack's in-plane free strain and zero out-of-plane stress the layers'
    # in-plane stresses (sig_I = Q_I (eps_I - alpha_I) at sig_O = 0) average to zero
    idx = [0, 1, 5]
    s_lay = [q_inplane(lam.S) @ (stack.alpha[idx] - lam.alpha[idx]) for lam in lams]
    s_in = sum(f * s for f, s in zip(p.sleeve_layer_fractions(), s_lay))
    scale = max(np.abs(s).max() for s in s_lay)
    assert scale > 0 and np.abs(s_in).max() < 1e-12 * scale


def test_one_layer_spec_equals_single_layer():
    a = default_params("sleeved_bridged", **SPEC, sleeve_layers="20:1")
    b = default_params("sleeved_bridged", **SPEC, sleeve_wind_angle=20.0)
    for k in c.SLEEVE_KEYS:
        assert getattr(a, k) == pytest.approx(getattr(b, k), rel=1e-12, abs=1e-20), k


def test_layer_stack_rejects():
    lam = c.angle_ply(c.lamina(c.FIBRES["carbon_sm"], c.MATRICES["epoxy"], 0.55), 0.0)
    with pytest.raises(ValueError):
        c.layer_stack([lam], (1, 2))
    with pytest.raises(ValueError):
        c.layer_stack([lam, lam], (1, 0))
    with pytest.raises(TypeError, match="LayerStack has no plies"):
        c.sleeve_failure(np.zeros(6), c.STRENGTHS["carbon_sm_epoxy"], c.layer_stack([lam], (1,)))


# ---- cache key, EM key --------------------------------------------------------------------------
def test_new_fields_key_and_mech_only():
    assert set(NEW_FIELDS) <= set(KEY_OMITTED_AT_DEFAULT) and set(NEW_FIELDS) <= MECH_ONLY
    p = default_params("sleeved_bridged")
    base = point_key(p, "mech")
    q = p.replace(**SPEC, sleeve_layers="0:2,20:1,0:1", sleeve_layer_prestress="50")
    assert em_model_key(q) == em_model_key(p)
    assert point_key(p.replace(sleeve_layer_prestress="50"), "mech") != base
    assert point_key(q, "mech") != point_key(q.replace(sleeve_layers="0:1,20:1,0:2"), "mech")
    assert point_key(p.replace(sleeve_layer_prestress=""), "mech") == base


# ---- closed form: N = 1 is the single-ring solution, bit for bit ---------------------------------
def _steel(alpha=12e-6):
    return an.isotropic(200e3, 0.3, 7850e-12, alpha=alpha)


def _cfrp(alpha_r=30e-6, alpha_t=0.0):
    return an.orthotropic(8.8e3, 125e3, 0.28, 1800e-12, alpha_r=alpha_r, alpha_t=alpha_t)


@pytest.mark.parametrize("w, delta, dT", [(W, 0.1, 0.0), (W, 0.1, 100.0), (0.0, 0.05, -40.0), (3 * W, 0.02, 0.0)])
def test_single_ring_bit_identical(w, delta, dT):
    D, R = _steel(), _cfrp()
    old = an.disc_ring(A, B, D, R, w, delta, dT)
    new = an.disc_rings([A, B], D, [R], w, delta, dT)
    assert np.array_equal(old["coeffs"], new.coeffs)
    assert old["p_interface"] == new.p_interface and old["in_contact"] == new.in_contact
    rr, rd = np.linspace(A, B, 9), np.linspace(1.0, A, 9)
    assert all(np.array_equal(x, y) for x, y in zip(old["ring"](rr), new.ring(0, rr)))
    assert all(np.array_equal(x, y) for x, y in zip(old["disc"](rd), new.disc(rd)))
    assert an.rings_coupled_pressure([A, B], D, [R], w, delta, dT) == an.coupled_pressure(A, B, D, R, w, delta, dT)
    assert an.rings_liftoff_speed([A, B], D, [R], delta, dT) == an.liftoff_speed(A, B, D, R, delta, dT)
    if w == 3 * W:
        assert not new.in_contact   # the separated branch is covered too


def test_ring_split_in_identical_layers_is_one_ring():
    """Three bonded rings of the same material are one ring: same fields to round-off."""
    D, R = _steel(), _cfrp()
    one = an.disc_rings([A, B], D, [R], W, 0.1, 50.0, [4e-4])
    three = an.disc_rings([A, 71.6, 72.0, B], D, [R, R, R], W, 0.1, 50.0, [4e-4] * 3)
    assert three.p_interface == pytest.approx(one.p_interface, rel=1e-9)
    for i, (r0, r1) in enumerate(zip(three.radii, three.radii[1:])):
        rr = np.linspace(r0, r1, 5)
        for x, y in zip(three.ring(i, rr), one.ring(0, rr)):
            assert np.allclose(x, y, rtol=1e-9, atol=1e-9)


# ---- closed form: ODE, continuity, boundary conditions -----------------------------------------
def _fd_residuals(body, f, r0, r1, w, fr, ft, n=41):
    """Equilibrium (central differences) and constitutive law with the stress-free strains fr, ft,
    relative to the largest stress on [r0, r1]."""
    h = 5e-3
    r = np.linspace(r0 + h, r1 - h, n)
    u, sr, st = f(r)
    scale = np.abs(np.r_[sr, st]).max()
    up = (f(r + h)[0] - f(r - h)[0]) / (2 * h)
    dsr = (f(r + h)[1] - f(r - h)[1]) / (2 * h)
    eq = dsr + (sr - st) / r + body.rho * w * w * r
    er = body.Q11 * (up - fr) + body.Q12 * (u / r - ft) - sr
    et = body.Q12 * (up - fr) + body.Q22 * (u / r - ft) - st
    return np.abs(eq).max() / scale, np.abs(er).max() / scale, np.abs(et).max() / scale


@pytest.fixture(scope="module")
def three_layer():
    p = default_params("sleeved_bridged", **{**MULTILAYER_CASE})
    radii = layer_radii(A, B, p.sleeve_layer_fractions())
    bodies = layer_bodies(p)
    dT = p.mech_delta_t
    eps = p.sleeve_layer_prestrain()
    sol = an.disc_rings(radii, disc_body(p), bodies, W, p.interference, dT, eps)
    return p, radii, bodies, dT, eps, sol


def test_multi_ring_satisfies_ode_continuity_and_bcs(three_layer):
    p, radii, bodies, dT, eps, sol = three_layer
    assert sol.in_contact and sol.n == 3
    D = disc_body(p)
    eq, er, et = _fd_residuals(D, sol.disc, 1.0, A, W, D.alpha_r * dT, D.alpha_t * dT)
    assert max(eq, er, et) < 1e-7, (eq, er, et)
    for i, body in enumerate(bodies):
        eq, er, et = _fd_residuals(body, lambda r, i=i: sol.ring(i, r), radii[i], radii[i + 1], W,
                                   body.alpha_r * dT, body.alpha_t * dT - eps[i])
        assert max(eq, er, et) < 1e-6, (i, eq, er, et)
    scale = 400.0
    ud, srd, _ = sol.disc(np.array([A]))
    u1, sr1, _ = sol.ring(0, np.array([A]))
    assert math.isclose(u1[0] - ud[0], p.interference, rel_tol=1e-9)
    assert abs(srd[0] - sr1[0]) < 1e-9 * scale and math.isclose(-sr1[0], sol.p_interface, rel_tol=1e-12)
    for i in range(2):
        r = np.array([radii[i + 1]])
        ui, si, ti = sol.ring(i, r)
        uo, so, to = sol.ring(i + 1, r)
        assert abs(ui[0] - uo[0]) < 1e-10 * abs(ui[0]) and abs(si[0] - so[0]) < 1e-9 * scale
        assert abs(ti[0] - to[0]) > 50.0       # the hoop stress jumps at a material interface
        assert sol.sr_interfaces[i] == pytest.approx(si[0], rel=1e-12)
    assert abs(sol.ring(2, np.array([B]))[1][0]) < 1e-9 * scale


def test_multi_ring_superposition_and_liftoff(three_layer):
    """Linear in the loads while in contact: the pressure is the sum of its parts, p(w) = p0 - c w^2."""
    p, radii, bodies, dT, eps, _ = three_layer
    D = disc_body(p)
    parts = [an.rings_coupled_pressure(radii, D, bodies, 0.0, p.interference, 0.0),
             an.rings_coupled_pressure(radii, D, bodies, 0.0, 0.0, dT),
             an.rings_coupled_pressure(radii, D, bodies, 0.0, 0.0, 0.0, eps),
             an.rings_coupled_pressure(radii, D, bodies, W, 0.0, 0.0)]
    total = an.rings_coupled_pressure(radii, D, bodies, W, p.interference, dT, eps)
    assert total == pytest.approx(sum(parts), rel=1e-10)
    assert all(x > 0 for x in parts[:3]) and parts[3] < 0
    w_lift = an.rings_liftoff_speed(radii, D, bodies, p.interference, dT, eps)
    assert abs(an.rings_coupled_pressure(radii, D, bodies, w_lift, p.interference, dT, eps)) < 1e-8
    sol = an.disc_rings(radii, D, bodies, 1.01 * w_lift, p.interference, dT, eps)
    assert not sol.in_contact and sol.p_interface == 0.0
    # prestress alone holds the sleeve on: lift-off without interference is finite and > 0
    assert 0 < an.rings_liftoff_speed(radii, D, bodies, 0.0, 0.0, eps) < w_lift


def test_free_ring_with_hoop_eigenstrain():
    """A free ring with a uniform hoop eigenstrain (stress-free hoop strain -eps0) is self-stressed;
    closed form against a numerical boundary-value solution with both faces traction-free, and the
    hoop force of a free ring vanishes."""
    from scipy.integrate import solve_bvp

    R, eps0 = _cfrp(), 8e-4
    sol = an.disc_rings([A, B], _steel(), [R], 0.0, -1.0, 0.0, [eps0])   # 1 mm clearance: free ring
    assert not sol.in_contact

    def f(r, y):
        u, sr = y
        up = (sr - R.Q12 * (u / r + eps0)) / R.Q11
        st = R.Q12 * up + R.Q22 * (u / r + eps0)
        return np.vstack([up, -(sr - st) / r])

    bvp = solve_bvp(f, lambda ya, yb: np.array([ya[1], yb[1]]), np.linspace(A, B, 50), np.zeros((2, 50)), tol=1e-12)
    r = np.linspace(A, B, 11)
    u, sr, st = sol.ring(0, r)
    assert u.max() < -0.05                                      # it shrinks by about eps0 r
    assert np.allclose(u, bvp.sol(r)[0], rtol=0, atol=1e-9)
    assert np.allclose(sr, bvp.sol(r)[1], rtol=0, atol=1e-6)
    assert abs(sol.mean_hoop(0)) < 1e-6 and abs(sol.mean_hoop(0, n=32)) < 1e-6
    assert np.abs(st).max() > 0.1                               # incompatible: not stress-free
    # with an isotropic ring and an equal radial eigenstrain it would be a free shrink, stress-free
    Ri = an.isotropic(70e3, 0.33, 2700e-12)
    u, sr, st = an.free_strain_terms(Ri, r, -eps0, -eps0)
    assert np.all(u == 0) and np.allclose(sr, st) and np.allclose(sr, sr[0])   # b_r = b_t: all in A r


def test_held_ring_carries_sigma0():
    """Definition of the prestress: held at zero hoop strain with no radial stress, the layer carries
    sigma0 = E_hoop eps0 (E_hoop = hoop_modulus)."""
    R = _cfrp()
    eps0 = 1e-3
    Et = an.hoop_modulus(R)
    assert Et == pytest.approx(125e3, rel=1e-14)
    # eps_t = 0, sigma_r = 0 -> eps_r from sigma_r = Q11 eps_r + Q12 (0 + eps0) = 0
    eps_r = -R.Q12 * eps0 / R.Q11
    st = R.Q12 * eps_r + R.Q22 * eps0
    assert st == pytest.approx(Et * eps0, rel=1e-14)


def test_rings_summary_and_reference(three_layer):
    p, radii, bodies, dT, eps, sol = three_layer
    s = multilayer_ring_reference(p, A, B - A, rpm=20000.0)
    assert s["radii"] == radii and s["p_interface"] == pytest.approx(sol.p_interface, rel=1e-12)
    assert len(s["layers"]) == 3 and s["sh_ring_max"] == max(d["sh_max"] for d in s["layers"])
    for i, d in enumerate(s["layers"]):
        assert d["sh_bore"] == pytest.approx(float(sol.ring(i, radii[i])[2]), rel=1e-12)
        assert d["sh_min"] <= d["sh_mean"] <= d["sh_max"]
    # the helical middle layer carries far less hoop stress than the hoop layers
    assert s["layers"][1]["sh_max"] < 0.5 * s["layers"][0]["sh_max"]


# ---- sequential winding ---------------------------------------------------------------------------
def test_sequential_winding_single_layer_is_exact():
    D, R = _steel(0.0), _cfrp(0.0, 0.0)
    out = an.sequential_winding([A, B], D, [R], 120.0)
    assert out["prestrain"][0] == pytest.approx(out["sequential"]["steps"][0], rel=1e-8)
    assert out["sequential"]["layers"][0]["sh_mean"] == pytest.approx(120.0, rel=1e-8)
    assert out["max_hoop_dev"] < 1e-6
    assert out["simultaneous"].mean_hoop(0) == pytest.approx(120.0, rel=1e-8)


def test_sequential_winding_three_layers():
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1")
    radii, D, bodies = p.sleeve_layer_radii(), disc_body(p), layer_bodies(p)
    out = an.sequential_winding(radii, D, bodies, 100.0)
    sim, seq = out["simultaneous"], out["sequential"]
    # the equivalent eigenstrains reproduce the radial stress at the disc and at every interface ...
    assert -sim.p_interface == pytest.approx(seq["sr"][0], rel=1e-8)
    for i in range(2):
        assert sim.sr_interfaces[i] == pytest.approx(seq["sr"][i + 1], rel=1e-8)
    # ... hence every layer's hoop force (mean hoop stress)
    for i in range(3):
        assert sim.mean_hoop(i) == pytest.approx(seq["layers"][i]["sh_mean"], rel=1e-8)
    # the last layer keeps its winding stress, the inner ones are relaxed by the layers above them
    assert seq["layers"][2]["sh_mean"] == pytest.approx(100.0, rel=1e-7)
    assert seq["layers"][0]["sh_mean"] < seq["layers"][1]["sh_mean"] < 100.0
    # within a layer the hoop gradient differs only slightly (thin layers)
    assert out["max_hoop_dev"] < 1e-2
    assert out["prestress"] == pytest.approx([e * an.hoop_modulus(b) for e, b in zip(out["prestrain"], bodies)])
    # linear in the winding stress, per layer
    two = an.sequential_winding(radii, D, bodies, [200.0, 200.0, 200.0])
    assert two["prestress"] == pytest.approx([2 * x for x in out["prestress"]], rel=1e-9)
    # sub-rings: still exact on the interface pressures, a coarser per-layer approximation of the gradient
    fine = an.sequential_winding(radii, D, bodies, 100.0, sublayers=4)
    assert -fine["simultaneous"].p_interface == pytest.approx(fine["sequential"]["sr"][0], rel=1e-8)
    assert fine["max_hoop_dev"] < 1.0
    with pytest.raises(ValueError):
        an.sequential_winding(radii, D, bodies, [1.0, 2.0])


def test_winding_prestress_from_params():
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:1,0:1,0:1")
    out = winding_prestress(p, 100.0)
    vals = [float(x) for x in out["spec"].split(",")]
    assert vals == pytest.approx(out["prestress"], rel=1e-5)
    q = p.replace(sleeve_layer_prestress=out["spec"])   # valid Params input
    assert q.sleeve_layer_prestress_mpa() == tuple(vals)


# ---- deck ---------------------------------------------------------------------------------------
def _deck(tmp_path, p, ring=False, **kw):
    model = mech_model_ring(A, B - A) if ring else mech_model_from_rotor(build_rotor(p))
    return write_deck(p, model, tmp_path / "d.inp", **kw).read_text(encoding="ascii")


def test_layered_deck_geometry_materials_mesh(tmp_path):
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1")
    lines = _deck(tmp_path, p).splitlines()
    consts = p.sleeve_layer_constants()
    # one area and material per layer, element axes from the cylindrical system 11
    for i in range(3):
        mat = MAT_LAYER0 + 1 + i
        k = lines.index(f"ASEL,S,AREA,,{301 + i}")
        assert lines[k + 1] == f"AATT,{mat},1,1,11"
        assert f"MP,EY,{mat},{consts[i]['sleeve_E_hoop'] * 1e3:.10g}" in lines
        assert f"MP,EX,{mat},{consts[i]['sleeve_E_radial'] * 1e3:.10g}" in lines
        nu_rh = consts[i]["sleeve_nu_hoop_radial"] * consts[i]["sleeve_E_radial"] / consts[i]["sleeve_E_hoop"]
        assert f"MP,PRXY,{mat},{nu_rh:.10g}" in lines
    assert not any(ln.startswith("MP,EX,3,") for ln in lines)   # the single-layer material is not written
    # shared keypoints and lines: 4 radii x 2 cut lines, 3 + 3 radial lines and 4 arcs, areas from them
    assert sum(ln.startswith("K,300") or ln.startswith("K,301") for ln in lines) == 8
    assert sum(ln.startswith("LARC,300") for ln in lines) == 4
    assert lines[lines.index("NUMSTR,AREA,302") + 1] == "AL,30001,30008,30004,30007"
    # at least two elements through each layer; the bore arc sized like the lamination
    # (radial lines 30000-30002 at theta 0, 30003-30005 at the sector angle; hf = 0.15 mm)
    for ln, nd in ((30000, 5), (30003, 5), (30001, 3), (30004, 3), (30002, 3), (30005, 3)):
        assert f"LESIZE,{ln},,,{nd},,1" in lines
    assert "LESIZE,30006,0.5,,,,1" in lines
    q = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1", mech_mesh_fine=0.5)
    assert "LESIZE,30000,,,2,,1" in _deck(tmp_path, q).splitlines()   # never fewer than two
    # contact target on the bore of the whole sleeve (all layer materials)
    assert lines.count(f"ESEL,S,MAT,,{MAT_LAYER0 + 1},{MAT_LAYER0 + 3}") == 1
    # per-layer post: sorts and the node export appended per layer
    for i in range(3):
        assert f"*GET,SH_SLV_L{i + 1},SORT,,MAX" in lines
        assert f"*VFILL,SLA(1),RAMP,{i + 1},0" in lines
    assert lines.count("*CFOPEN,sipm_sleeve_nodes,txt") == 1
    assert lines.count("*CFOPEN,sipm_sleeve_nodes,txt,,APPEND") == 2
    assert lines.count("*VWRITE,SRA(1),STA(1),SSX(1),SSY(1),SSZ(1),SSXY(1),SLA(1)") == 3
    assert not any("INISTATE" in ln for ln in lines) and not any(ln.startswith("MP,ALPX") for ln in lines)


def test_layered_deck_thermal_and_prestress(tmp_path):
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1", sleeve_layer_prestress="100,0,80",
                       mech_temperature=120.0)
    lines = _deck(tmp_path, p, ring=True).splitlines()
    consts = p.sleeve_layer_constants()
    for i in range(3):
        assert f"MP,ALPY,{MAT_LAYER0 + 1 + i},{consts[i]['sleeve_cte_hoop']:.10g}" in lines
    eps = p.sleeve_layer_prestrain()
    k = lines.index("INISTATE,SET,CSYS,-2")
    assert lines[k + 1] == "INISTATE,SET,DTYP,EPEL"
    define = "INISTATE,DEFINE,ALL,ALL,,,0,{:.10g},0,0,0,0"
    assert lines[k + 2: k + 7] == [f"ESEL,S,MAT,,{MAT_LAYER0 + 1}", define.format(eps[0]),
                                   f"ESEL,S,MAT,,{MAT_LAYER0 + 3}", define.format(eps[2]), "ALLSEL"]
    # in /SOLU, after the loads, before the solution controls
    assert lines.index("/SOLU") < lines.index("BFUNIF,TEMP,120") < k < lines.index("AUTOTS,ON") < lines.index("SOLVE")


@pytest.mark.parametrize("topo", ["sleeved_bridged", "ring"])
def test_single_layer_prestress_only_adds_inistate(tmp_path, topo):
    """A single-layer sleeve with prestress: the default deck plus the INISTATE block on material 3."""
    p = default_params("sleeved_bridged")
    plain = _deck(tmp_path, p, ring=topo == "ring").splitlines()
    lines = _deck(tmp_path, p.replace(sleeve_layer_prestress="125"), ring=topo == "ring").splitlines()
    k = lines.index("INISTATE,SET,CSYS,-2")
    block = lines[k - 1: k + 5]
    assert block == ["! winding prestress: hoop eigenstrain per sleeve layer as initial elastic strain (element axes)",
                     "INISTATE,SET,CSYS,-2", "INISTATE,SET,DTYP,EPEL", "ESEL,S,MAT,,3",
                     "INISTATE,DEFINE,ALL,ALL,,,0,0.001,0,0,0,0", "ALLSEL"]
    assert lines[:k - 1] + lines[k + 5:] == plain
    # zero prestress writes nothing: the default deck
    assert _deck(tmp_path, p.replace(sleeve_layer_prestress="0"), ring=topo == "ring").splitlines() == plain


def test_default_deck_hash_unchanged(tmp_path):
    """Same pin as tests/test_load_cases.py: the layered code does not touch the default deck."""
    text = _deck(tmp_path, default_params("sleeved_bridged"))
    assert hashlib.sha256(text.encode()).hexdigest() == \
        "0dbe06c6d82f923ce5b554ae71bc9c0f1be707e23066a045eb49c26e912ea9e5"


# ---- node export, per-layer failure and interface shear ------------------------------------------
LAYERED = dict(**MULTILAYER_CASE)   # the fixture is a 3-layer ring export of this case (theta = 0 nodes)


def test_parse_layered_node_fixture():
    s = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes_layered.txt")
    assert s["layer"].tolist() == [1, 1, 1, 2, 2, 2, 3, 3, 3] and s["layer"].dtype.kind == "i"
    # interface nodes appear once per layer with that layer's own hoop stress
    assert s["r"][2] == s["r"][3] == 71.85 and s["s_hoop"][2] > 2 * s["s_hoop"][3]
    assert s["s_radial"][2] == pytest.approx(s["s_radial"][3], rel=1e-3)   # traction continuous


def test_layered_failure_per_layer():
    p = default_params("sleeved_bridged", **LAYERED)
    nodes = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes_layered.txt")
    out = sleeve_failure_outputs(p, nodes)
    lams, st = p.sleeve_layer_laminates(), p.sleeve_strength()
    sig = np.zeros((9, 6))
    sig[:, 0], sig[:, 1], sig[:, 2], sig[:, 4] = nodes["s_hoop"], nodes["s_axial"], nodes["s_radial"], nodes["s_shear"]
    per = [c.sleeve_failure(sig[3 * k: 3 * k + 3], st, lams[k], p.mech_delta_t) for k in range(3)]
    for crit in c.CRITERIA[:1] + c.CRITERIA[2:]:
        assert out[f"sleeve_fi_{crit}"] == pytest.approx(max(f[crit].max() for f in per), rel=1e-14)
    for k in range(3):
        assert out[f"sleeve_sf_L{k + 1}"] == pytest.approx(per[k]["tsai_wu_sr"].min(), rel=1e-14)
    sf_all = min(out[f"sleeve_sf_L{k + 1}"] for k in range(3))
    assert out["sleeve_sf"] == sf_all and out["sleeve_fi_tsai_wu"] == pytest.approx(1 / sf_all, rel=1e-14)
    # the +/-30 deg layer governs (as in the MAPDL run), at the node where it does
    assert out["sleeve_sf_layer"] == 2 and out["sleeve_sf_r"] in nodes["r"][3:6]
    # evaluating the interface node with the wrong layer's laminate would give another value
    wrong = c.sleeve_failure(sig[3:4], st, lams[0], p.mech_delta_t)["tsai_wu_sr"][0]
    assert not math.isclose(wrong, per[1]["tsai_wu_sr"][0], rel_tol=1e-3)


def test_layered_failure_needs_layer_column():
    p = default_params("sleeved_bridged", **LAYERED)
    nodes = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes_layered.txt")
    nodes.pop("layer")
    with pytest.raises(ValueError, match="layer"):
        sleeve_failure_outputs(p, nodes)


def test_interface_shear():
    nodes = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes_layered.txt")
    t = nodes["s_shear"]
    expected = max(abs(0.5 * (t[2] + t[3])), abs(0.5 * (t[5] + t[6])))
    assert sleeve_interface_outputs(nodes) == {"tau_sleeve_interface_max": pytest.approx(expected, rel=1e-14)}
    assert sleeve_interface_outputs(parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes.txt")) == {}
    one = {k: v[:3] for k, v in nodes.items()}
    assert sleeve_interface_outputs(one) == {}


# ---- CLI and configuration -------------------------------------------------------------------------
def test_cli_composite_layers(capsys):
    assert main(["composite", "--fibre", "carbon_sm", "--layers", "0:2,20:1,0:1", "--toml"]) == 0
    out = capsys.readouterr().out
    rows = {ln.split()[0]: ln.split()[2:] for ln in out.splitlines() if ln.startswith("  sleeve_")}
    assert set(rows) == set(c.SLEEVE_KEYS) and all(len(v) == 4 for v in rows.values())
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:2,20:1,0:1")
    assert float(rows["sleeve_E_hoop"][1]) == pytest.approx(p.sleeve_layer_constants()[1]["sleeve_E_hoop"], rel=1e-5)
    assert float(rows["sleeve_E_hoop"][3]) == pytest.approx(p.sleeve_E_hoop, rel=1e-5)
    spec = tomllib.loads(out[out.index("# [base]"):])
    assert spec["sleeve_layers"] == "0:2,20:1,0:1"
    default_params("sleeved_bridged", **spec)
    with pytest.raises(SystemExit):
        main(["composite", "--fibre", "carbon_sm", "--layers", "0:1,20:1", "--angle", "10"])


def test_cli_winding(capsys):
    assert main(["winding", "--tension", "100", "--set", "sleeve_fibre=carbon_sm", "sleeve_layers=0:1,0:1,0:1"]) == 0
    out = capsys.readouterr().out
    spec = tomllib.loads(out[out.index("sleeve_layer_prestress"):])["sleeve_layer_prestress"]
    p = default_params("sleeved_bridged", **SPEC, sleeve_layers="0:1,0:1,0:1")
    assert [float(x) for x in spec.split(",")] == pytest.approx(winding_prestress(p, 100.0)["prestress"], rel=1e-5)
    with pytest.raises(SystemExit):
        main(["winding", "--tension", "1,2", "--set", "sleeve_fibre=carbon_sm", "sleeve_layers=0:1,0:1,0:1"])


def test_cli_set_keeps_string_fields():
    from sipm.cli import _parse_sets

    assert _parse_sets(["sleeve_layer_prestress=100", "interference=0.2"]) == dict(sleeve_layer_prestress="100",
                                                                                 interference=0.2)


def test_example_toml_resolves():
    cfg = tomllib.loads((Path(__file__).parents[1] / "examples/configs/mech_multilayer_sleeve.toml").read_text())
    sw = sweep_from_config(cfg)
    plan = sw.plan()
    assert all(err is None for *_, err in plan)
    ps = [p for _, _, p, _, _ in plan]
    assert {p.sleeve_layer_count for p in ps} >= {1, 3}
    assert any(p.sleeve_layer_prestress for p in ps) and any(p.is_layered and p.interference == 0 for p in ps)
    assert len({k for _, _, _, k, _ in plan}) == len(plan)


def test_labels_for_new_columns():
    from sipm.loadcases import ENVELOPE
    from sipm.post import LABELS

    assert ("tau_sleeve_interface_max", "max") in ENVELOPE
    for k in ("tau_sleeve_interface_max", "sleeve_sf_layer", "sh_sleeve_max_L1", f"sleeve_sf_L{MAX_SLEEVE_LAYERS}",
              "tau_sleeve_interface_max__worst", "sh_sleeve_max_L2__overspeed_hot"):
        assert k in LABELS, k
