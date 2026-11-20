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
@pytest.mark.parametrize("case", ["isotropic", "orthotropic"])
def test_ring_against_closed_form(cfg, case):
    from sipm.validation import ring_check

    r = ring_check(case, cfg=cfg, workdir=cfg.workdir / f"ring_{case}")
    for row in r["rows"]:
        print(f"[ring {case}] {row['quantity']}: FE {row['fe']:.4f} closed form {row['closed_form']:.4f} "
              f"error {100 * row['rel_error']:+.3f} %")
        assert abs(row["rel_error"]) < 5e-3


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
