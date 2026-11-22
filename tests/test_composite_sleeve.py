"""Composite sleeve integration: Params spec and derivation, validation, sweep cache key,
failure outputs of the mech point, and the `sipm composite` command."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

from sipm import composite as c
from sipm.cli import main
from sipm.mapdl import parse_sleeve_nodes, sleeve_failure_outputs
from sipm.maxwell import MECH_ONLY, em_model_key
from sipm.params import Params, coerce, default_params, field_names
from sipm.sweep import KEY_OMITTED_AT_DEFAULT, point_key, resolve

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

FIX = Path(__file__).parent / "fixtures"
SPEC_FIELDS = ("sleeve_fibre", "sleeve_matrix", "sleeve_vf", "sleeve_micromechanics", "sleeve_wind_angle")
STRENGTH_FIELDS = ("sleeve_Xt", "sleeve_Xc", "sleeve_Yt", "sleeve_Yc", "sleeve_S")
THERMAL_FIELDS = ("lam_cte", "pm_cte", "sleeve_cte_hoop", "sleeve_cte_radial", "sleeve_cte_axial",
                  "mech_temperature", "mech_ref_temperature")
LOAD_CASE_FIELDS = ("mech_cases", "overspeed_factor", "mech_temperature_cold")   # tests/test_load_cases.py
NEW_FIELDS = THERMAL_FIELDS + SPEC_FIELDS + STRENGTH_FIELDS + LOAD_CASE_FIELDS


# ---- 1. cache key ------------------------------------------------------------------------------
# point keys computed once with the code of `main` before this change (git archive main src)
MAIN_KEYS = {
    ("sleeved_bridged", "em"): "02f290e812750495",
    ("sleeved_bridged", "mech"): "3f9e2cce34901fa1",
    ("sleeved_bridged", "acloss"): "85062abd023983df",
    ("sleeved_bridgeless", "em"): "65b6195107de3c12",
    ("sleeved_bridgeless", "mech"): "5ba426d0dd188b43",
    ("sleeved_bridgeless", "acloss"): "631621a12076c718",
}


@pytest.mark.parametrize("topo_phys", list(MAIN_KEYS))
def test_default_point_key_unchanged_from_main(topo_phys):
    topo, phys = topo_phys
    assert point_key(default_params(topo), phys) == MAIN_KEYS[topo_phys]


def test_variant_point_key_unchanged_from_main():
    p = default_params("sleeved_bridged", sleeve_thickness=2.0, interference=0.3, gamma=120.0)
    assert point_key(p, "mech") == "d06f95f8a7c8f1ac"


def test_new_fields_omitted_only_at_default():
    assert set(KEY_OMITTED_AT_DEFAULT) == set(NEW_FIELDS)
    p = default_params("sleeved_bridged")
    base = point_key(p, "mech")
    for q in (p.replace(mech_temperature=80.0), p.replace(sleeve_Xt=1500.0), p.replace(lam_cte=11e-6),
              p.replace(sleeve_fibre="carbon_sm"), p.replace(sleeve_matrix="peek"), p.replace(sleeve_vf=0.6)):
        assert point_key(q, "mech") != base
    # a composite spec changes the derived constants too, and different specs differ
    a = point_key(p.replace(sleeve_fibre="carbon_sm", sleeve_wind_angle=15.0), "mech")
    b = point_key(p.replace(sleeve_fibre="carbon_sm", sleeve_wind_angle=30.0), "mech")
    assert a != b


def test_new_fields_mech_only_and_em_key():
    assert set(NEW_FIELDS) <= MECH_ONLY
    p = default_params("sleeved_bridged")
    q = p.replace(sleeve_fibre="carbon_im", sleeve_matrix="peek", sleeve_vf=0.6, sleeve_micromechanics="rom",
                  sleeve_wind_angle=20.0, sleeve_Xt=1000.0, sleeve_S=50.0, mech_temperature=90.0)
    assert em_model_key(p) == em_model_key(q)


# ---- 2. composite spec in Params ---------------------------------------------------------------
def test_defaults_are_table_ii_constants():
    p = default_params("sleeved_bridged")
    assert p.sleeve_fibre == "" and p.sleeve_laminate() is None
    assert (p.sleeve_E_hoop, p.sleeve_E_radial, p.sleeve_nu_hoop_radial, p.sleeve_density) == (125.0, 8.8, 0.28, 1800.0)
    spec = (p.sleeve_matrix, p.sleeve_vf, p.sleeve_micromechanics, p.sleeve_wind_angle)
    assert spec == ("epoxy", 0.55, "chamis", 0.0)


@pytest.mark.parametrize("angle", [0.0, 15.0, 30.0, 45.0, 90.0])
def test_spec_derives_all_thirteen_constants(angle):
    p = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_vf=0.6, sleeve_wind_angle=angle)
    lam = c.angle_ply(c.lamina(c.FIBRES["carbon_sm"], c.MATRICES["epoxy"], 0.6, "chamis"), angle)
    for k, v in c.sleeve_params(lam).items():
        assert getattr(p, k) == pytest.approx(v, rel=1e-14, abs=1e-20), k
    assert isinstance(p.sleeve_laminate(), c.Laminate) and p.sleeve_laminate().angle_deg == angle


def test_spec_overwrites_passed_constants_and_replace_rederives():
    p = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_E_hoop=999.0, sleeve_cte_hoop=5e-5)
    assert p.sleeve_E_hoop == pytest.approx(128.07, abs=0.01) and abs(p.sleeve_cte_hoop) < 1e-6
    q = p.replace(sleeve_vf=0.6)
    assert q.sleeve_E_hoop > p.sleeve_E_hoop
    assert p.replace(sleeve_E_hoop=1.0).sleeve_E_hoop == p.sleeve_E_hoop   # still derived
    r = q.replace(sleeve_fibre="")      # back to direct constants: the last values stay as entered
    assert r.sleeve_laminate() is None and r.sleeve_E_hoop == q.sleeve_E_hoop


def test_vf_0536_against_table_ii():
    p = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_vf=0.536)
    assert p.sleeve_E_hoop == pytest.approx(125.0, rel=1e-3)
    assert p.sleeve_E_radial == pytest.approx(7.98, abs=0.01)
    assert p.sleeve_nu_hoop_radial == pytest.approx(0.2696, abs=1e-4)
    assert p.sleeve_density == pytest.approx(1521.6)


@pytest.mark.parametrize("bad", [
    dict(sleeve_wind_angle=15.0),                         # angle without a fibre
    dict(sleeve_fibre="kevlar"),
    dict(sleeve_fibre="carbon_sm", sleeve_matrix="pps"),
    dict(sleeve_fibre="carbon_sm", sleeve_micromechanics="mori_tanaka"),
    dict(sleeve_fibre="carbon_sm", sleeve_vf=0.0),
    dict(sleeve_fibre="carbon_sm", sleeve_vf=1.0),
    dict(sleeve_fibre="carbon_sm", sleeve_wind_angle=-10.0),
    dict(sleeve_fibre="carbon_sm", sleeve_wind_angle=95.0),
    dict(sleeve_Xt=0.0),
    dict(sleeve_S=-1.0),
])
def test_spec_invalid_rejected(bad):
    with pytest.raises(ValueError):
        default_params("sleeved_bridged", **bad)


def test_string_fields_coerce_and_resolve_from_table():
    assert coerce("sleeve_fibre", "carbon_hm") == "carbon_hm" and coerce("sleeve_vf", "0.6") == 0.6
    p = resolve({"sleeve_wind_angle": 30, "mech_temperature": 120},
                base={"sleeve_fibre": "carbon_sm", "sleeve_matrix": "epoxy", "sleeve_vf": 0.55})
    assert p.sleeve_wind_angle == 30.0 and p.sleeve_laminate() is not None
    assert p.sleeve_nu_hoop_axial > 1   # a 30 deg angle ply, accepted (see test below)


def test_example_toml_resolves():
    cfg = tomllib.loads((Path(__file__).parents[1] / "examples/configs/mech_composite_sleeve.toml").read_text())
    from sipm.cli import sweep_from_config

    sw = sweep_from_config(cfg)
    plan = sw.plan()
    assert len(plan) == 8 and all(err is None for *_, err in plan)
    assert {p.sleeve_wind_angle for _, _, p, _, _ in plan} == {0.0, 15.0, 30.0, 45.0}
    assert len({k for _, _, _, k, _ in plan}) == 8


# ---- 3. compliance-based validation ------------------------------------------------------------
def test_angle_ply_constants_outside_half_accepted_directly():
    """Constants of a +/-30 deg winding (nu_hoop_axial > 1, nu_hoop_radial < 0) entered directly
    pass, because the 6x6 compliance is positive definite."""
    q = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_wind_angle=30.0)
    direct = {k: getattr(q, k) for k in c.SLEEVE_KEYS}
    p = default_params("sleeved_bridged", **direct)
    assert p.sleeve_nu_hoop_axial > 1 and p.sleeve_nu_hoop_radial < 0
    assert np.linalg.eigvalsh(p.sleeve_compliance()).min() > 0
    assert np.allclose(p.sleeve_compliance(), q.sleeve_laminate().S, rtol=1e-10, atol=1e-14)


@pytest.mark.parametrize("bad", [
    dict(sleeve_nu_hoop_radial=4.0),                   # nu_hr^2 Er/Eh > 1
    dict(sleeve_nu_radial_axial=1.2),                  # transverse pair, E_r = E_a: |nu| < 1 needed
    dict(sleeve_nu_hoop_axial=2.5, sleeve_nu_hoop_radial=2.5),
])
def test_non_positive_definite_sleeve_rejected(bad):
    with pytest.raises(ValueError, match="positive definite"):
        default_params("sleeved_bridged", **bad)


def test_negative_poisson_allowed_when_spd():
    p = default_params("sleeved_bridged", sleeve_nu_radial_axial=-0.2)
    assert p.sleeve_nu_radial_axial == -0.2


# ---- 4. strengths ------------------------------------------------------------------------------
def test_strength_defaults_and_object():
    p = default_params("sleeved_bridged")
    assert p.sleeve_strength() == c.STRENGTHS["carbon_sm_epoxy"]
    s = p.replace(sleeve_Yt=40.0).sleeve_strength()
    assert s.Yt == 40.0 and s.St is None


# ---- 5. failure outputs ------------------------------------------------------------------------
def test_failure_outputs_direct_constants_hoop_ply():
    p = default_params("sleeved_bridged")
    nodes = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes.txt")
    out = sleeve_failure_outputs(p, nodes)
    st = p.sleeve_strength()
    # hoop ply: 11 = hoop, 22 = axial (0 here), 33 = radial, 13 = r-theta shear
    a = c.failure_indices(nodes["s_hoop"], nodes["s_axial"], 0.0, st)
    b = c.failure_indices(nodes["s_hoop"], nodes["s_radial"], nodes["s_shear"], st)
    for crit in ("max_stress", "hashin_fibre", "hashin_matrix"):
        assert out[f"sleeve_fi_{crit}"] == pytest.approx(max(a[crit].max(), b[crit].max()), rel=1e-14)
    sr = np.minimum(a["tsai_wu_sr"], b["tsai_wu_sr"])
    i = int(np.argmin(sr))
    assert out["sleeve_sf"] == pytest.approx(sr[i]) and out["sleeve_sf_r"] == nodes["r"][i]
    # Tsai-Wu column = 1 / strength ratio (max over nodes), not the raw quadratic (negative here)
    assert out["sleeve_fi_tsai_wu"] == pytest.approx(1.0 / out["sleeve_sf"], rel=1e-14)
    assert out["sleeve_fi_tsai_wu"] == pytest.approx(np.max(1.0 / sr), rel=1e-14)
    assert max(a["tsai_wu"].max(), b["tsai_wu"].max()) < 0 < out["sleeve_fi_tsai_wu"]
    assert out["sleeve_sf_theta"] == nodes["theta"][i]
    # sanity: 224 MPa hoop against Xt = 2000 MPa and a few MPa radial compression
    assert 0.1 < out["sleeve_fi_max_stress"] < 0.2 and out["sleeve_sf"] > 1


def test_failure_outputs_spec_at_zero_angle_equals_hoop_ply():
    nodes = parse_sleeve_nodes(FIX / "mapdl_sleeve_nodes.txt")
    spec = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_vf=0.536, mech_temperature=120.0)
    direct = default_params("sleeved_bridged", **{k: getattr(spec, k) for k in c.SLEEVE_KEYS},
                            mech_temperature=120.0)
    a, b = sleeve_failure_outputs(spec, nodes), sleeve_failure_outputs(direct, nodes)
    for k in a:
        assert a[k] == pytest.approx(b[k], rel=1e-9, abs=1e-12), k


def test_failure_outputs_angle_ply_thermal_residual():
    """+/-30 deg plies heated with no mechanical stress: the ply mismatch alone loads them."""
    zero = {k: np.zeros(3) for k in ("s_radial", "s_hoop", "s_axial", "s_shear")}
    zero.update(r=np.array([71.1, 71.8, 72.6]), theta=np.zeros(3))
    cold = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_wind_angle=30.0)
    hot = cold.replace(mech_temperature=120.0)
    zc = sleeve_failure_outputs(cold, zero)
    assert zc["sleeve_fi_max_stress"] == 0.0 and zc["sleeve_fi_tsai_wu"] == 0.0 and math.isinf(zc["sleeve_sf"])
    out = sleeve_failure_outputs(hot, zero)
    assert out["sleeve_fi_max_stress"] > 0.05 and math.isfinite(out["sleeve_sf"])
    # the same heating in a hoop winding leaves the plies unstressed
    assert sleeve_failure_outputs(hot.replace(sleeve_wind_angle=0.0), zero)["sleeve_fi_max_stress"] < 1e-12


# ---- CLI ---------------------------------------------------------------------------------------
def test_cli_composite_prints_constants_and_toml(capsys):
    assert main(["composite", "--fibre", "carbon_sm", "--vf", "0.536", "--toml"]) == 0
    out = capsys.readouterr().out
    vals = {ln.split()[0]: float(ln.split()[2]) for ln in out.splitlines() if ln.startswith("  sleeve_")}
    assert set(vals) == set(c.SLEEVE_KEYS)
    assert vals["sleeve_E_hoop"] == pytest.approx(124.904, abs=1e-3)
    assert vals["sleeve_density"] == pytest.approx(1521.6)
    toml_part = out[out.index("# [base]"):]
    spec = tomllib.loads(toml_part)
    assert spec == dict(sleeve_fibre="carbon_sm", sleeve_matrix="epoxy", sleeve_vf=0.536,
                        sleeve_micromechanics="chamis", sleeve_wind_angle=0.0)
    p = default_params("sleeved_bridged", **spec)
    assert p.sleeve_E_hoop == pytest.approx(vals["sleeve_E_hoop"], rel=1e-5)
    # the commented direct constants are complete and valid Params input
    direct = tomllib.loads("\n".join(ln[2:] for ln in toml_part.splitlines() if ln.startswith("# sleeve_")))
    assert set(direct) == set(c.SLEEVE_KEYS)
    default_params("sleeved_bridged", **direct)


def test_cli_composite_angle_and_errors(capsys):
    assert main(["composite", "--fibre", "carbon_sm", "--angle", "30", "--model", "halpin_tsai"]) == 0
    assert "+/-30 deg" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(["composite", "--fibre", "basalt"])
    with pytest.raises(SystemExit):
        main(["composite", "--fibre", "carbon_sm", "--vf", "1.5"])


def test_params_json_carries_string_fields(capsys):
    """MATLAB's sipm.defaultParams reads `sipm params --json`; string fields round-trip via --set."""
    import json

    assert main(["params", "--json", "--set", "sleeve_fibre=carbon_sm", 'sleeve_matrix="peek"',
                 "sleeve_wind_angle=15"]) == 0
    d = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert d["sleeve_fibre"] == "carbon_sm" and d["sleeve_matrix"] == "peek" and d["sleeve_wind_angle"] == 15
    assert set(d) == set(field_names()) and isinstance(Params(**d), Params)
