"""Thermal load case (params, deck, closed form) and the sleeve stress export (deck, parser)."""

import math
import re
from pathlib import Path

import numpy as np
import pytest

from sipm.analytic import Body, disc_ring, isotropic, orthotropic, summary
from sipm.geometry import build_rotor
from sipm.mapdl import (
    SLEEVE_COLUMNS,
    SLEEVE_FILE,
    mech_model_from_rotor,
    mech_model_ring,
    parse_sleeve_nodes,
    write_deck,
)
from sipm.maxwell import MECH_ONLY, em_model_key
from sipm.params import CTE_FIELDS, default_params

FIX = Path(__file__).parent / "fixtures"
A, B = 71.1, 72.6
W = 20000 * 2 * math.pi / 60
NEW_FIELDS = (*CTE_FIELDS, "mech_temperature", "mech_ref_temperature")


# ---- params ---------------------------------------------------------------------------------
def test_thermal_defaults_isothermal():
    p = default_params("sleeved_bridged")
    assert p.mech_temperature == p.mech_ref_temperature and p.mech_delta_t == 0.0
    assert (p.lam_cte, p.pm_cte, p.sleeve_cte_hoop, p.sleeve_cte_radial, p.sleeve_cte_axial) == (
        12e-6, 2e-6, 0.0, 30e-6, 30e-6)
    assert math.isclose(p.replace(mech_temperature=120.0).mech_delta_t, 100.0)


@pytest.mark.parametrize("bad", [
    dict(lam_cte=12.0),             # 1e-6/K value without the factor
    dict(sleeve_cte_radial=30.0),
    dict(pm_cte=-1e-4),
    dict(sleeve_cte_hoop=float("nan")),
    dict(mech_temperature=-300.0),
    dict(mech_ref_temperature=-273.15),
    dict(mech_temperature=600.0),
])
def test_thermal_invalid_rejected(bad):
    with pytest.raises(ValueError):
        default_params("sleeved_bridged", **bad)


def test_thermal_valid_extremes_accepted():
    default_params("sleeved_bridged", sleeve_cte_hoop=-1e-6, mech_temperature=-40.0, mech_ref_temperature=150.0)


def test_new_fields_do_not_change_em_key():
    assert set(NEW_FIELDS) <= MECH_ONLY
    p = default_params("sleeved_bridged")
    q = p.replace(lam_cte=11e-6, pm_cte=5e-6, sleeve_cte_hoop=-0.5e-6, sleeve_cte_radial=25e-6,
                  sleeve_cte_axial=26e-6, mech_temperature=150.0, mech_ref_temperature=25.0)
    assert em_model_key(p) == em_model_key(q)


# ---- deck -----------------------------------------------------------------------------------
def _deck(tmp_path, topo="sleeved_bridged", **kw):
    p = default_params(topo, **kw)
    return p, write_deck(p, mech_model_from_rotor(build_rotor(p)), tmp_path / "x.inp").read_text()


def test_deck_isothermal_has_no_thermal_commands(tmp_path):
    _, d = _deck(tmp_path)
    assert "MP,ALP" not in d and "TREF" not in d and "BFUNIF" not in d


def test_deck_thermal_commands(tmp_path):
    p, d = _deck(tmp_path, mech_temperature=120.0, mech_ref_temperature=25.0, sleeve_cte_hoop=-0.5e-6)
    get = lambda k, m: float(re.search(rf"^MP,{k},{m},(\S+)$", d, re.M).group(1))  # noqa: E731
    for k in ("ALPX", "ALPY", "ALPZ"):
        assert math.isclose(get(k, 1), 12e-6) and math.isclose(get(k, 2), 2e-6)
    # sleeve in its element (cylindrical) system: x radial, y hoop, z axial
    assert math.isclose(get("ALPX", 3), 30e-6)
    assert math.isclose(get("ALPY", 3), -0.5e-6)
    assert math.isclose(get("ALPZ", 3), 30e-6)
    lines = d.splitlines()
    assert "TREF,25" in lines and "BFUNIF,TEMP,120" in lines
    # TREF in /PREP7, before the sleeve's element system is assigned; the temperature in /SOLU
    assert lines.index("/PREP7") < lines.index("TREF,25") < lines.index("AATT,3,1,1,11")
    assert lines.index("/SOLU") < lines.index("BFUNIF,TEMP,120") < lines.index("SOLVE")


def test_deck_ring_thermal(tmp_path):
    p = default_params("sleeved_bridged", mech_temperature=120.0)
    d = write_deck(p, mech_model_ring(A, 1.5), tmp_path / "r.inp").read_text()
    assert "TREF,20" in d and "BFUNIF,TEMP,120" in d and "MP,ALPY,3,0" in d


@pytest.mark.parametrize("topo", ["sleeved_bridged", "sleeved_bridgeless"])
def test_deck_sleeve_export(tmp_path, topo):
    _, d = _deck(tmp_path, topo)
    lines = d.splitlines()
    i = lines.index("*VWRITE,SRA(1),STA(1),SSX(1),SSY(1),SSZ(1),SSXY(1)")
    assert lines[i + 1] == "(6E21.12)" and lines[i + 2] == "*CFCLOS"
    stem, ext = SLEEVE_FILE.split(".")
    assert f"*CFOPEN,{stem},{ext}" in lines
    block = lines[lines.index("NSLE,S,CORNER") - 1:i]
    assert block[0] == "ESEL,S,MAT,,3"
    assert "RSYS,11" in block and "CSYS,11" in block   # stresses and node locations cylindrical
    for comp in ("X", "Y", "Z", "XY"):
        assert any(ln.startswith("*VGET,") and ln.endswith(f",NODE,1,S,{comp}") for ln in block), comp
    assert lines[i + 3] == "CSYS,0"
    # new scalars: max radial tension and both extremes of the r-theta shear, in RSYS 11
    for key, how in (("SRT_SLV", "MAX"), ("TAU_SLV_HI", "MAX"), ("TAU_SLV_LO", "MIN")):
        assert f"*GET,{key},SORT,,{how}" in d
        assert f"('{key.lower()} = ',E20.12)" in d
    for old in ("sh_slv", "sr_slv", "svm_slv", "sr_lam"):   # existing keys unchanged
        assert f"('{old} = ',E20.12)" in d


def test_parse_sleeve_nodes_fixture():
    s = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes.txt")
    assert tuple(s) == SLEEVE_COLUMNS
    n = len(s["r"])
    assert n == 6 and all(v.shape == (n,) for v in s.values())
    assert math.isclose(s["r"][0], 71.1) and math.isclose(s["theta"][-1], 22.5)
    assert math.isclose(s["s_hoop"].max(), 223.920281347, rel_tol=1e-12)
    assert math.isclose(s["s_radial"].min(), -4.63362799091, rel_tol=1e-12)
    assert math.isclose(s["s_shear"][1], -1.25e-2, rel_tol=1e-12)
    assert np.all(s["s_axial"] == 0.0)


def test_parse_sleeve_nodes_rejects_wrong_columns(tmp_path):
    f = tmp_path / "bad.txt"
    f.write_text("1.0 2.0 3.0\n")
    with pytest.raises(ValueError, match="columns"):
        parse_sleeve_nodes(f)


# ---- closed form ----------------------------------------------------------------------------
STEEL = dict(E=200e3, nu=0.3, rho=7850e-12)


def _cf_ring(alpha_r=30e-6, alpha_t=0.0):
    return orthotropic(8.8e3, 125e3, 0.28, 1800e-12, alpha_r=alpha_r, alpha_t=alpha_t)


def _residuals(body: Body, f, r0, r1, w, dT, n=41):
    """Equilibrium residual (central differences, h = 5e-3 mm) and the constitutive law with
    thermal strain, on [r0, r1], each divided by the largest stress there. The FD error is
    ~1e-9 (disc, truncation) and ~1e-8 (ring, round-off of the r^+-k basis) of that scale."""
    r = np.linspace(r0, r1, n)
    h = 5e-3
    u, sr, st = f(r)
    scale = np.abs(np.r_[sr, st]).max()
    up = (f(r + h)[0] - f(r - h)[0]) / (2 * h)
    dsr = (f(r + h)[1] - f(r - h)[1]) / (2 * h)
    eq = dsr + (sr - st) / r + body.rho * w * w * r
    er = (body.Q11 * (up - body.alpha_r * dT) + body.Q12 * (u / r - body.alpha_t * dT)) - sr
    et = (body.Q12 * (up - body.alpha_r * dT) + body.Q22 * (u / r - body.alpha_t * dT)) - st
    return np.abs(eq).max() / scale, np.abs(er).max() / scale, np.abs(et).max() / scale


@pytest.mark.parametrize("dT", [100.0, -60.0])
def test_thermal_satisfies_ode_bcs_and_interference(dT):
    D = isotropic(**STEEL, alpha=12e-6)
    R = _cf_ring()
    s = disc_ring(A, B, D, R, W, 0.1, dT)
    assert s["in_contact"]
    scale = abs(s["p_interface"])
    for body, f, r0, r1 in ((D, s["disc"], 1.0, A), (R, s["ring"], A, B)):
        eq, er, et = _residuals(body, f, r0, r1, W, dT)
        assert eq < 1e-7 and er < 1e-7 and et < 1e-7, (eq, er, et, body)
    ud, srd, _ = s["disc"](np.array([A]))
    ur, sra, _ = s["ring"](np.array([A]))
    assert math.isclose(ur[0] - ud[0], 0.1, rel_tol=1e-10)
    assert math.isclose(srd[0], sra[0], rel_tol=1e-10) and math.isclose(-sra[0], s["p_interface"], rel_tol=1e-10)
    assert abs(s["ring"](np.array([B]))[1][0]) < 1e-9 * scale


def test_thermal_reduces_to_isothermal():
    D0, R0 = isotropic(**STEEL), orthotropic(8.8e3, 125e3, 0.28, 1800e-12)
    D, R = isotropic(**STEEL, alpha=12e-6), _cf_ring()
    a0, a1 = summary(A, B, D0, R0, W, 0.1), summary(A, B, D, R, W, 0.1, dT=0.0)
    assert a0 == a1   # bit-identical: the thermal terms vanish exactly at dT = 0
    s = disc_ring(A, B, D, R, W, 0.1, 1e-9)
    assert math.isclose(s["p_interface"], a0["p_interface"], rel_tol=1e-9)


def test_isotropic_free_expansion_is_stress_free():
    alpha, dT = 12e-6, 150.0
    D = isotropic(**STEEL, alpha=alpha)
    R = isotropic(70e3, 0.33, 2700e-12, alpha=alpha)
    s = disc_ring(A, B, D, R, 0.0, 0.0, dT)
    for f, r in ((s["disc"], np.linspace(1.0, A, 20)), (s["ring"], np.linspace(A, B, 20))):
        u, sr, st = f(r)
        assert np.allclose(u, alpha * dT * r, rtol=1e-12, atol=0)
        assert np.abs(sr).max() < 1e-9 and np.abs(st).max() < 1e-9


def test_orthotropic_ring_free_expansion():
    """A free orthotropic ring heated with alpha_r != alpha_t is not stress-free (uniform
    thermal strains with alpha_r != alpha_t are incompatible in a ring); check the closed form
    against a direct numerical boundary-value solution with both faces traction-free."""
    from scipy.integrate import solve_bvp

    R = _cf_ring()
    dT = 100.0
    s = disc_ring(A, B, isotropic(**STEEL), R, 0.0, -1.0, dT)   # large clearance: parts separate
    assert not s["in_contact"]

    def f(r, y):
        u, sr = y
        up = (sr - R.Q12 * u / r + R.beta_r * dT) / R.Q11
        st = R.Q12 * up + R.Q22 * u / r - R.beta_t * dT
        return np.vstack([up, -(sr - st) / r])

    sol = solve_bvp(f, lambda ya, yb: np.array([ya[1], yb[1]]), np.linspace(A, B, 50), np.zeros((2, 50)), tol=1e-12)
    r = np.linspace(A, B, 11)
    u, sr, st = s["ring"](r)
    assert np.abs(u).max() > 1e-3                              # not trivially zero
    assert np.allclose(u, sol.sol(r)[0], rtol=0, atol=1e-9)   # mm
    assert np.allclose(sr, sol.sol(r)[1], rtol=0, atol=1e-6)  # MPa
    assert np.allclose(st, R.Q12 * (np.gradient(u, r) - R.alpha_r * dT) + R.Q22 * (u / r - R.alpha_t * dT),
                       rtol=0, atol=1e-2)


@pytest.mark.parametrize("ratio", [1.0, 1.0 + 1e-12, 1.0 + 1e-6, 1.0 - 1e-6])
def test_degenerate_k_equal_one(ratio):
    """Q22 = Q11 with alpha_r != alpha_t: resonant r ln r particular term, continuous across k = 1."""
    E = 50e3
    R = orthotropic(E, E * ratio, 0.3, 1800e-12, alpha_r=30e-6, alpha_t=0.0)
    D = isotropic(**STEEL, alpha=12e-6)
    s = disc_ring(A, B, D, R, W, 0.1, 100.0)
    eq, er, et = _residuals(R, s["ring"], A, B, W, 100.0)
    assert eq < 1e-7 and er < 1e-7 and et < 1e-7, (eq, er, et)
    ref = disc_ring(A, B, D, orthotropic(E, E, 0.3, 1800e-12, alpha_r=30e-6, alpha_t=0.0), W, 0.1, 100.0)
    assert math.isclose(s["p_interface"], ref["p_interface"], rel_tol=1e-5)


def test_heating_raises_sleeve_stress():
    """Steel grows, the fibre direction does not: heating adds interference."""
    D, R = isotropic(**STEEL, alpha=12e-6), _cf_ring()
    cold, hot = summary(A, B, D, R, W, 0.1), summary(A, B, D, R, W, 0.1, dT=100.0)
    assert hot["p_interface"] > cold["p_interface"] and hot["sh_ring_max"] > cold["sh_ring_max"]
    # thermal interference ~ (alpha_steel - alpha_hoop) dT a = 0.085 mm, comparable with 0.1 mm
    assert 1.5 < hot["sh_ring_max"] / cold["sh_ring_max"] < 2.5


def test_ring_reference_uses_params():
    from sipm.validation import ring_params, ring_reference

    p0, p1 = ring_params("orthotropic"), ring_params("thermal")
    assert p0.mech_delta_t == 0 and p1.mech_delta_t == 100.0
    r0, r1 = ring_reference(p0, A, 1.5), ring_reference(p1, A, 1.5)
    assert r1["sh_ring_max"] > r0["sh_ring_max"]
    with pytest.raises(ValueError):
        ring_params("hot")
