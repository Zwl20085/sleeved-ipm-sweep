"""Integration tests (need licensed Ansys). Run with:  pytest -m ansys -s

Set SIPM_WORKDIR to choose where AEDT/MAPDL files go (default <tempdir>/sipm_runs).
"""

import math

import pytest

from sipm.config import Config, ansys_available
from sipm.params import TOPOLOGIES, default_params

pytestmark = pytest.mark.ansys

need_mapdl = pytest.mark.skipif(not ansys_available("mapdl"), reason="MAPDL not found")
need_aedt = pytest.mark.skipif(not ansys_available("aedt"), reason="AEDT not found")


@pytest.fixture(scope="module")
def cfg():
    c = Config()
    return Config(workdir=c.workdir / "integration_tests")


@need_mapdl
@pytest.mark.parametrize("case", ["isotropic", "orthotropic", "thermal"])
def test_ring_against_closed_form(cfg, case):
    import numpy as np

    from sipm.validation import ring_check

    r = ring_check(case, cfg=cfg, workdir=cfg.workdir / f"ring_{case}")
    for row in r["rows"]:
        print(f"[ring {case}, dT = {r['delta_t']:g} K] {row['quantity']}: FE {row['fe']:.4f} closed form "
              f"{row['closed_form']:.4f} error {100 * row['rel_error']:+.3f} %")
        assert abs(row["rel_error"]) < 5e-3
    # sleeve node export: cylindrical locations, extremes equal to the scalar results
    fe, s = r["fe"], r["fe"]["_sleeve_nodes"]
    assert np.allclose([s["r"].min(), s["r"].max()], [r["radius"], r["radius"] + r["t_sleeve"]], atol=1e-6)
    assert np.allclose([s["theta"].min(), s["theta"].max()], [0.0, 22.5], atol=1e-6)
    assert np.isclose(s["s_hoop"].max(), fe["sh_sleeve_max"], rtol=1e-9)
    assert np.isclose(s["s_radial"].min(), fe["sr_sleeve_min"], rtol=1e-9)
    assert np.isclose(s["s_radial"].max(), fe["sr_sleeve_max"], rtol=1e-9)
    assert np.isclose(np.abs(s["s_shear"]).max(), fe["tau_sleeve_max"], rtol=1e-9)
    # axisymmetric ring: no radial tension and no shear beyond discretisation noise
    assert fe["sr_sleeve_max"] < 1e-3 * fe["sh_sleeve_max"]
    assert fe["tau_sleeve_max"] < 5e-3 * fe["sh_sleeve_max"]


@need_mapdl
@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_mech_point(cfg, topo):
    from sipm.mapdl import solve_mech

    r = solve_mech(default_params(topo), cfg, cfg.workdir / f"mech_{topo}")
    print(f"[mech {topo}] svm_rotor_max={r['svm_rotor_max']:.1f} MPa (lam {r['svm_lam_max']:.1f} at r="
          f"{r['svm_lam_r']:.2f} mm, theta={r['svm_lam_theta']:.2f} deg; PM {r['svm_pm_max']:.1f}) "
          + (f"sleeve hoop max {r['sh_sleeve_max']:.1f} MPa, open contact fraction {r['sleeve_open_fraction']:.3f} "
             if "sh_sleeve_max" in r else "")
          + f"stabilisation/strain energy {r.get('stab_energy_ratio', float('nan')):.2e} "
          f"nodes {r['n_nodes']} wall {r['mech_wall_s']:.1f} s")
    assert r["svm_rotor_max"] > 0
    if "stab_energy_ratio" in r:
        assert r["stab_energy_ratio"] < 1e-3
    s = r["_sleeve_nodes"]
    print(f"[mech {topo}] sleeve: {len(s['r'])} corner nodes, sr max {r['sr_sleeve_max']:.3f} MPa, "
          f"|tau_rt| max {r['tau_sleeve_max']:.3f} MPa")
    assert abs(s["s_hoop"].max() - r["sh_sleeve_max"]) <= 1e-9 * r["sh_sleeve_max"]
    assert r["sr_sleeve_min"] <= r["sr_sleeve_max"] and r["tau_sleeve_max"] >= 0
    print(f"[mech {topo}] sleeve failure: max-stress {r['sleeve_fi_max_stress']:.3f}, "
          f"Tsai-Wu {r['sleeve_fi_tsai_wu']:.3f}, Hashin fibre {r['sleeve_fi_hashin_fibre']:.3f}, "
          f"matrix {r['sleeve_fi_hashin_matrix']:.3f}, "
          f"SF {r['sleeve_sf']:.2f} at r={r['sleeve_sf_r']:.2f} mm, theta={r['sleeve_sf_theta']:.2f} deg")
    assert r["sleeve_sf"] > 1 and 0 < r["sleeve_fi_max_stress"] < 1


@need_mapdl
def test_mech_point_composite_angle_ply(cfg):
    """A +/-30 deg winding (nu_hoop_axial > 1, nu_hoop_radial < 0) is accepted by MAPDL and solves."""
    from sipm.mapdl import solve_mech

    p = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_wind_angle=30.0)
    assert p.sleeve_nu_hoop_axial > 1 and p.sleeve_nu_hoop_radial < 0
    r = solve_mech(p, cfg, cfg.workdir / "mech_composite_30")
    print(f"[mech composite +/-30] sleeve hoop {r['sh_sleeve_max']:.1f} MPa, Tsai-Wu SF {r['sleeve_sf']:.2f}, "
          f"svm_lam_max {r['svm_lam_max']:.1f} MPa")
    assert r["svm_rotor_max"] > 0 and r["sleeve_sf"] > 0


@need_mapdl
def test_mech_point_thermal(cfg):
    """Heating the default bridge-retained rotor raises the sleeve hoop stress (steel grows, the
    fibre direction does not); the isothermal default is unchanged by the thermal option."""
    from sipm.mapdl import solve_mech

    p = default_params("sleeved_bridged")
    cold = solve_mech(p, cfg, cfg.workdir / "mech_thermal_cold")
    hot = solve_mech(p.replace(mech_temperature=p.mech_ref_temperature + 100.0), cfg, cfg.workdir / "mech_thermal_hot")
    print(f"[mech thermal] sleeve hoop {cold['sh_sleeve_max']:.1f} -> {hot['sh_sleeve_max']:.1f} MPa, "
          f"svm_lam_max {cold['svm_lam_max']:.1f} -> {hot['svm_lam_max']:.1f} MPa, contact p max "
          f"{cold['contact_p_max']:.2f} -> {hot['contact_p_max']:.2f} MPa")
    assert math.isclose(cold["sh_sleeve_max"], 223.920281347, rel_tol=1e-6)  # docs/validation.md section 2
    assert hot["sh_sleeve_max"] > 1.3 * cold["sh_sleeve_max"]
    assert hot["contact_p_min"] > cold["contact_p_min"]


@pytest.fixture(scope="module")
def session(cfg):
    from sipm.maxwell import MaxwellSession

    with MaxwellSession(cfg) as s:
        yield s


@need_aedt
@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_em_point(session, topo):
    r = session.solve(default_params(topo))
    print(f"[em {topo}] gamma=130 deg, 420 Arms: torque_avg={r['torque_avg']:.2f} N m, ripple pk-pk/avg="
          f"{r['torque_ripple']:.3f}, wall {r['em_wall_s']:.1f} s (incl. model build)")
    assert r["torque_avg"] > 0  # motoring: torque along the direction of rotation


@need_aedt
def test_em_open_circuit(session):
    from sipm.validation import open_circuit_check

    p = default_params("sleeved_bridged")
    chk = open_circuit_check(session.open_circuit(p), p)
    for k, v in chk["phases"].items():
        print(f"[open circuit] phase {k}: psi phase {v['psi_phase']:.4f} deg (predicted {v['axis_predicted']:.1f}), "
              f"EMF amplitude {v['emf_amp']:.3f} V, half-wave asymmetry {v['half_wave_asym']:.2e}")
    assert all(abs(v["phase_error"]) < 0.5 for v in chk["phases"].values())
    assert all(v["half_wave_asym"] < 1e-3 for v in chk["phases"].values())
    assert chk["emf_amp_spread"] < 1e-3
    assert math.isclose(chk["phase_shift_AB"], 120, abs_tol=1.0)
