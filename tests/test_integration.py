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


@need_mapdl
def test_ring_liftoff_against_closed_form(cfg):
    """FE lift-off speed of the disc + orthotropic sleeve, extrapolated element by element from
    standstill and 20000 rpm, against the closed form, isothermal and 100 K hot."""
    from sipm.validation import ring_liftoff_check

    r = ring_liftoff_check(cfg=cfg, workdir=cfg.workdir / "ring_liftoff")
    for row in r["rows"]:
        print(f"[ring lift-off {row['label']}, {row['temperature']:g} degC] FE {row['fe']:.1f} rpm, closed form "
              f"{row['closed_form']:.1f} rpm, error {100 * row['rel_error']:+.4f} % ({row['basis']})")
        assert abs(row["rel_error"]) < 1e-3
    fe = r["fe"]
    assert fe["mech_case"] == "rated_hot" and fe["sh_sleeve_max"] == fe["sh_sleeve_max__rated_hot"]
    assert math.isclose(fe["sh_sleeve_max__rated_cold"], 192.4089, rel_tol=1e-5)   # = the single-case run


@need_mapdl
def test_mech_load_cases(cfg):
    """The six cases of the default bridge-retained point, 120 degC hot: the rated cold case is the
    single-case default point, the envelope picks the expected cases."""
    from sipm.mapdl import solve_mech
    from sipm.params import MECH_CASES

    p = default_params("sleeved_bridged", mech_temperature=120.0, mech_cases=",".join(MECH_CASES))
    r = solve_mech(p, cfg, cfg.workdir / "mech_load_cases")
    for c in MECH_CASES:
        print(f"[mech cases] {c:<15s} svm_lam_max {r[f'svm_lam_max__{c}']:7.1f}  sh_sleeve_max "
              f"{r[f'sh_sleeve_max__{c}']:7.1f}  contact p min {r[f'contact_p_min__{c}']:6.3f}  "
              f"open {r[f'sleeve_open_fraction__{c}']:.3f}  SF {r[f'sleeve_sf__{c}']:.2f}")
    print(f"[mech cases] lift-off cold {r['liftoff_rpm_cold']:.0f} rpm ({r['liftoff_basis_cold']}), "
          f"hot {r['liftoff_rpm_hot']:.0f} rpm ({r['liftoff_basis_hot']}), wall {r['mech_wall_s']:.1f} s")
    assert math.isclose(r["sh_sleeve_max__rated_cold"], 223.920281347, rel_tol=1e-6)   # single default point
    assert r["mech_case"] == "rated_hot" and r["sh_sleeve_max"] == r["sh_sleeve_max__rated_hot"]
    col = lambda k, c: r[f"{k}__{c}"]  # noqa: E731
    for t in ("cold", "hot"):
        # centrifugal load: sleeve hoop stress rises and the interface pressure falls with speed
        sh = [col("sh_sleeve_max", f"{s}_{t}") for s in ("standstill", "rated", "overspeed")]
        pm = [col("contact_p_min", f"{s}_{t}") for s in ("standstill", "rated", "overspeed")]
        assert sh[0] < sh[1] < sh[2] and pm[0] > pm[1] > pm[2]
        assert col("svm_lam_max", f"rated_{t}") < col("svm_lam_max", f"overspeed_{t}")
    for s in ("standstill", "rated", "overspeed"):
        # heating adds interference (steel grows, the fibre direction does not)
        assert col("sh_sleeve_max", f"{s}_hot") > col("sh_sleeve_max", f"{s}_cold")
        assert col("contact_p_min", f"{s}_hot") > col("contact_p_min", f"{s}_cold")
    # the envelope is the extreme over the case columns
    for k, how in (("svm_lam_max", max), ("sh_sleeve_max", max), ("sleeve_fi_tsai_wu", max),
                   ("sleeve_sf", min), ("contact_p_min", min)):
        vals = {c: col(k, c) for c in MECH_CASES}
        assert r[f"{k}__worst"] == how(vals.values()) and vals[r[f"{k}__worst_case"]] == r[f"{k}__worst"]
    # names only where the margin is large: hoop stress hot vs cold +50 % (overspeed vs rated +5 %, but
    # monotonic in speed, asserted above); contact pressure cold vs hot -50 %, overspeed vs rated -13 %
    assert r["sh_sleeve_max__worst_case"] == "overspeed_hot"
    assert r["contact_p_min__worst_case"] == "overspeed_cold"
    assert r["liftoff_rpm_hot"] > r["liftoff_rpm_cold"] > p.overspeed_factor * p.mech_rpm   # 51.6k, 39.4k, 24k


@need_mapdl
def test_mech_liftoff_direct(cfg):
    """A direct run at the extrapolated lift-off speed of the default bridge-retained rotor opens the
    first contact element (min pressure ~0), a run at 0.95 x that speed has none open."""
    from sipm.validation import liftoff_check

    r = liftoff_check(default_params("sleeved_bridged"), "cold", cfg=cfg, workdir=cfg.workdir / "mech_liftoff")
    at, below = r["direct"]["at"], r["direct"]["below"]
    print(f"[mech lift-off] {r['liftoff_rpm']:.1f} rpm ({r['basis']}); p_min at standstill "
          f"{r['p_min_standstill']:.4f} MPa, at lift-off {at['contact_p_min']:.5f} MPa "
          f"(open {at['sleeve_open_fraction']:.4f}), at 0.95x {below['contact_p_min']:.4f} MPa "
          f"(open {below['sleeve_open_fraction']:.4f})")
    # docs/validation.md section 11: 0.026 MPa (1.0 % of the standstill minimum), one element open
    assert at["sleeve_open_fraction"] > 0 and abs(at["contact_p_min"]) < 2e-2 * r["p_min_standstill"]
    assert below["sleeve_open_fraction"] == 0 and below["contact_p_min"] > 0


@need_mapdl
def test_ring_multilayer_against_closed_form(cfg):
    """3-layer sleeve (hoop / +/-30 deg / hoop) with interference, winding prestress, 20000 rpm and
    100 K: every layer's hoop stress, the interface stresses and the contact pressure against the
    disc + N rings closed form, and the lift-off speed cold and hot (docs/validation.md section 12)."""
    import numpy as np

    from sipm.validation import multilayer_ring_check

    r = multilayer_ring_check(cfg=cfg, workdir=cfg.workdir / "ring_multilayer")
    for row in r["rows"]:
        print(f"[ring 3 layers] {row['quantity']}: FE {row['fe']:.4f} closed form {row['closed_form']:.4f} "
              f"error {100 * row['rel_error']:+.3f} %")
        assert abs(row["rel_error"]) < 5e-3
    for row in r["liftoff"]:
        print(f"[ring 3 layers lift-off {row['label']}] FE {row['fe']:.1f} rpm, closed form {row['closed_form']:.1f} "
              f"rpm, error {100 * row['rel_error']:+.4f} % ({row['basis']})")
        assert abs(row["rel_error"]) < 1e-3
    fe, s = r["fe"], r["fe"]["_sleeve_nodes"]
    assert sorted(set(s["layer"].tolist())) == [1, 2, 3]
    assert fe["sh_sleeve_max"] == max(fe[f"sh_sleeve_max_L{k}"] for k in (1, 2, 3))
    assert fe["sleeve_sf"] == min(fe[f"sleeve_sf_L{k}"] for k in (1, 2, 3))
    assert fe["sleeve_sf_layer"] == 2          # the helical layer governs
    for k in (1, 2, 3):
        assert np.isclose(s["s_hoop"][s["layer"] == k].max(), fe[f"sh_sleeve_max_L{k}"], rtol=1e-9)
    # axisymmetric: the interlaminar shear is discretisation noise
    assert fe["tau_sleeve_interface_max"] < 1e-3 * fe["sh_sleeve_max"]


@need_mapdl
def test_ring_prestress_single_layer(cfg):
    """Winding prestress on the single-layer sleeve (INISTATE on material 3), no interference, hot,
    against the N = 1 closed form with hoop eigenstrain."""
    from sipm import analytic
    from sipm.mapdl import mech_model_ring, solve_mech
    from sipm.validation import disc_body, layer_bodies, ring_params

    p = ring_params("thermal", interference=0.0, sleeve_layer_prestress="150", mech_mesh=0.5, mech_mesh_fine=0.5)
    fe = solve_mech(p, cfg, cfg.workdir / "ring_prestress", model=mech_model_ring(71.1, 1.5))
    ref = analytic.rings_summary([71.1, 72.6], disc_body(p), layer_bodies(p), p.mech_rpm * 2 * math.pi / 60, 0.0,
                                 p.mech_delta_t, p.sleeve_layer_prestrain())
    for name, f, c in (("sleeve max hoop", fe["sh_sleeve_max"], ref["sh_ring_max"]),
                       ("contact p min", fe["contact_p_min"], ref["p_interface"]),
                       ("disc max von Mises", fe["svm_lam_max"], ref["svm_disc_max"])):
        print(f"[ring prestress] {name}: FE {f:.4f} closed form {c:.4f} error {100 * (f - c) / c:+.3f} %")
        assert abs(f - c) < 5e-3 * abs(c)
    assert "sleeve_sf_layer" not in fe and "tau_sleeve_interface_max" not in fe   # single-layer columns only


@need_mapdl
def test_mech_point_multilayer(cfg):
    """The default bridge-retained rotor with a hoop / helical / hoop sleeve and winding prestress
    solves; the envelope columns are the extremes of the per-layer ones."""
    from sipm.mapdl import solve_mech

    p = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_layers="0:2,20:1,0:1",
                       sleeve_layer_prestress="50")
    r = solve_mech(p, cfg, cfg.workdir / "mech_multilayer")
    hoop = [r[f"sh_sleeve_max_L{k}"] for k in (1, 2, 3)]
    print(f"[mech 3 layers] svm_lam_max {r['svm_lam_max']:.1f} MPa, sleeve hoop per layer "
          f"{', '.join(f'{h:.1f}' for h in hoop)} MPa, contact p min {r['contact_p_min']:.3f} MPa, "
          f"SF {r['sleeve_sf']:.2f} (layer {r['sleeve_sf_layer']}), "
          f"interlaminar shear {r['tau_sleeve_interface_max']:.3f} MPa, nodes {r['n_nodes']}, "
          f"wall {r['mech_wall_s']:.1f} s")
    assert r["sh_sleeve_max"] == max(hoop) and hoop[1] < hoop[0]
    assert r["sleeve_sf"] == min(r[f"sleeve_sf_L{k}"] for k in (1, 2, 3)) and r["sleeve_sf"] > 1
    assert r["sleeve_sf_layer"] in (1, 2, 3) and r["sleeve_open_fraction"] == 0
    assert 0 < r["tau_sleeve_interface_max"] < r["tau_sleeve_max"] + 1e-9
    s = r["_sleeve_nodes"]
    assert sorted(set(s["layer"].tolist())) == [1, 2, 3]


@need_mapdl
def test_mech_point_multilayer_bridgeless(cfg):
    """A hoop / helical / hoop sleeve with winding prestress on the bridgeless rotor (0.5 mm sleeve, layers
    0.25 / 0.125 / 0.125 mm) solves with the per-layer outputs (docs/validation.md section 12)."""
    from sipm.mapdl import solve_mech

    p = default_params("sleeved_bridgeless", sleeve_fibre="carbon_sm", sleeve_layers="0:2,20:1,0:1",
                       sleeve_layer_prestress="50")
    r = solve_mech(p, cfg, cfg.workdir / "mech_multilayer_bridgeless")
    hoop = [r[f"sh_sleeve_max_L{k}"] for k in (1, 2, 3)]
    print(f"[mech 3 layers bridgeless] svm_lam_max {r['svm_lam_max']:.1f} MPa ({r['svm_lam_where']}), sleeve hoop "
          f"per layer {', '.join(f'{h:.1f}' for h in hoop)} MPa, contact p min {r['contact_p_min']:.3f} MPa, "
          f"SF {r['sleeve_sf']:.2f} (layer {r['sleeve_sf_layer']}), interlaminar shear "
          f"{r['tau_sleeve_interface_max']:.3f} MPa, nodes {r['n_nodes']}")
    assert r["sh_sleeve_max"] == max(hoop) and hoop[1] < min(hoop[0], hoop[2])
    assert r["sleeve_sf"] == min(r[f"sleeve_sf_L{k}"] for k in (1, 2, 3)) and r["sleeve_sf"] > 1
    assert r["sleeve_open_fraction"] == 0 and r["contact_p_min"] > 0
    assert sorted(set(r["_sleeve_nodes"]["layer"].tolist())) == [1, 2, 3]
    assert 0 < r["tau_sleeve_interface_max"] <= r["tau_sleeve_max"] + 1e-9


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
def test_em_point_conducting_sleeve(session):
    """2D sleeve eddy loss at 1e4 S/m axial (docs/validation.md section 13): converged period average,
    loss of the validated size, torque lowered by about the drag torque loss / omega."""
    p0 = default_params("sleeved_bridged")
    r0 = session.solve(p0)
    r = session.solve(p0.replace(sleeve_sigma_transverse=1e4))
    drag = r["sleeve_loss_w"] / (2 * math.pi * p0.em_rpm / 60)
    print(f"[em conducting sleeve] loss {r['sleeve_loss_w']:.4f} W (change over last period "
          f"{r['sleeve_loss_rel_change']:.2e}), torque {r0['torque_avg']:.4f} -> {r['torque_avg']:.4f} N m, "
          f"drag torque {drag:.4f} N m, wall {r['em_wall_s']:.1f} s")
    assert r["sleeve_loss_rel_change"] < 1e-3
    assert 4.8 < r["sleeve_loss_w"] < 5.3
    assert r["sleeve_sigma_axial"] == 1e4 and r["sleeve_skin_ratio"] > 10
    assert 0.5 * drag < r0["torque_avg"] - r["torque_avg"] < 1.5 * drag


@need_aedt
def test_em_point_conducting_sleeve_bridgeless(session):
    """The conducting-sleeve path on the bridgeless rotor (0.5 mm sleeve, 1e4 S/m axial): converged loss,
    torque lowered by about the drag torque (docs/validation.md section 13)."""
    p0 = default_params("sleeved_bridgeless")
    r0 = session.solve(p0)
    r = session.solve(p0.replace(sleeve_sigma_transverse=1e4))
    drag = r["sleeve_loss_w"] / (2 * math.pi * p0.em_rpm / 60)
    print(f"[em conducting sleeve bridgeless] loss {r['sleeve_loss_w']:.4f} W "
          f"(change {r['sleeve_loss_rel_change']:.2e}), torque {r0['torque_avg']:.4f} -> {r['torque_avg']:.4f} N m, "
          f"drag torque {drag:.4f} N m")
    assert r["sleeve_loss_rel_change"] < 1e-3
    assert 3.0 < r["sleeve_loss_w"] < 3.7          # 3.352 W measured
    assert 0.5 * drag < r0["torque_avg"] - r["torque_avg"] < 1.5 * drag


@need_aedt
def test_sleeve_loss_check(session):
    """`sipm validate sleeve-loss`: FE loss against the thin-shell estimate at 1e4 S/m, ~0.993 at load
    (docs/validation.md section 13 (c))."""
    from sipm.validation import sleeve_loss_check

    r = sleeve_loss_check(session=session)
    print(f"[sleeve loss check] FE {r['fe_loss']:.4f} W, estimate {r['estimate']:.4f} W, FE/estimate {r['ratio']:.4f}")
    assert abs(r["ratio"] - 0.993) < 5e-3


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
