import math

import pytest

from sipm.params import TOPOLOGIES, Params, coerce, default_params, field_names, parameter_markdown, parameter_table


@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_defaults_valid(topo):
    p = default_params(topo)
    assert p.topology == topo
    p.validate()


def test_airgap_fixed_sleeve_reduces_lamination():
    for t in (0.5, 1.0, 1.5, 2.0):
        p = default_params("sleeved_bridged", sleeve_thickness=t)
        assert math.isclose(p.r_sleeve_out, 147.84 / 2 - 1.32)
        assert math.isclose(p.r_lam, 147.84 / 2 - 1.32 - t)


def test_paper_values_present():
    p = default_params("sleeved_bridged")
    assert (p.stator_id, p.slots, p.poles, p.slot_opening, p.airgap, p.hairpin_layers) == (147.84, 48, 8, 1.9, 1.32, 7)
    assert (p.lam_density, p.lam_E, p.lam_nu, p.pm_density, p.pm_E, p.pm_nu) == (7850, 200, 0.3, 7400, 160, 0.24)
    assert (p.sleeve_density, p.sleeve_E_hoop, p.sleeve_E_radial) == (1800, 125, 8.8)
    assert p.mech_rpm == 20000 and p.current_rms == 420


def test_maxwell_values_present():
    """Defaults are the saved states of the Maxwell models (V angles 2 x RD_thetapm)."""
    p = default_params("sleeved_bridged")
    assert (p.stator_od, p.stack_length, p.shaft_diameter) == (210.0, 130.0, 60.0)
    assert math.isclose(p.pm1_angle / 2, 48.5039792373095, rel_tol=1e-15)
    assert math.isclose(p.pm2_angle / 2, 47.5144523955196, rel_tol=1e-15)
    assert (p.pm1_d_crossing, p.pm2_d_crossing) == (56.7305951647269, 45.0864009239181)
    assert (round(p.pm1_width, 3), round(p.pm2_width, 3)) == (3.649, 4.306)  # paper Table I: 3.65 / 4.31
    t = default_params("sleeved_bridgeless")
    assert (t.sleeve_thickness, t.pm1_rib_half_width, t.pm2_web_width, t.lip_depth, t.pm_end_clearance) == (
        0.5, 0.1, 0.1, 1.0, 0.0)


def test_minor_poisson_from_symmetry():
    p = default_params("sleeved_bridged")
    assert math.isclose(p.sleeve_nu_radial_hoop, 0.28 * 8.8 / 125)


@pytest.mark.parametrize("bad", [
    dict(topology="conventional"),
    dict(topology="sleeved_bridged", sleeve_thickness=0.0),
    dict(topology="sleeved_bridged", interference=-0.1),
    dict(topology="sleeved_bridged", slots=47),
    dict(topology="sleeved_bridged", pm1_angle=190.0),
    dict(topology="sleeved_bridged", plane="axisym"),
    dict(topology="sleeved_bridged", lam_nu=0.6),
    dict(topology="sleeved_bridged", sleeve_thickness=50.0),
    dict(topology="sleeved_bridged", pm_side_clearance=5.0),
    dict(topology="sleeved_bridgeless", lip_depth=0.0),
    dict(topology="sleeved_bridged", pm2_inner_tilt=50.0),
    dict(topology="sleeved_bridged", magnet_contact="glued"),
    dict(topology="sleeved_bridged", pm1_length=13.0),
    dict(topology="banana"),
])
def test_invalid_rejected(bad):
    topo = bad.pop("topology")
    with pytest.raises((ValueError, KeyError)):
        if topo in TOPOLOGIES:
            default_params(topo, **bad)
        else:
            Params(topology=topo)


def test_unknown_override():
    with pytest.raises(KeyError):
        default_params("sleeved_bridged", sleeve_thicknes=1.0)
    with pytest.raises(KeyError):
        default_params("sleeved_bridged").replace(nope=1)


def test_replace_validates():
    p = default_params("sleeved_bridged")
    with pytest.raises(ValueError):
        p.replace(airgap=-1.0)
    q = p.replace(sleeve_thickness=2.0)
    assert q.sleeve_thickness == 2.0 and p.sleeve_thickness == 1.5


def test_coerce_types():
    assert coerce("slots", "48") == 48 and isinstance(coerce("slots", 48.0), int)
    with pytest.raises(ValueError):
        coerce("slots", 47.5)
    assert coerce("gamma", "130") == 130.0
    assert coerce("topology", "sleeved_bridged") == "sleeved_bridged"


def test_every_field_documented():
    rows = parameter_table()
    assert {r["name"] for r in rows} == set(field_names())
    for r in rows:
        assert r["source"] in ("maxwell", "paper", "illustrative"), r
        assert r["doc"], r
    md = parameter_markdown()
    for n in field_names():
        assert f"`{n}`" in md


def test_child_environment_strips_matlab_settings():
    import os

    from sipm.config import child_environment

    base = {"MW_RESOURCE_PATH": r"C:\M\R2024b", "KMP_STACKSIZE": "512k", "OMP_NUM_THREADS": "4",
            "PATH": os.pathsep.join([r"C:\M\R2024b\extern\bin\win64", r"C:\Tools"])}
    env = child_environment(base)
    assert "KMP_STACKSIZE" not in env and env["OMP_NUM_THREADS"] == "4"
    assert env["PATH"] == r"C:\Tools" or env["PATH"].endswith("Tools")
    viasystem = {"MATLAB_ARCH": "win64", "KMP_STACKSIZE": "512k", "PATH": "x"}
    assert "KMP_STACKSIZE" not in child_environment(viasystem)
    plain = {"KMP_STACKSIZE": "1g", "PATH": "x"}
    assert child_environment(plain) == plain  # untouched outside MATLAB
