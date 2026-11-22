"""Conducting sleeve (2D axial-current eddy loss): conductivity fields and the axial-conductivity
helper, cache keys, the averaging and thin-shell estimate functions, the conducting EM run wiring
(fake AEDT session), resolution checks, line-report parsing and the CLI."""

import math
from pathlib import Path

import numpy as np
import pytest

from sipm import composite as c
from sipm import eddy
from sipm import maxwell as mx
from sipm.cli import load_config, main, sweep_from_config
from sipm.maxwell import MECH_ONLY, em_model_key
from sipm.params import default_params
from sipm.sweep import KEY_OMITTED_AT_DEFAULT, point_key

ROOT = Path(__file__).resolve().parents[1]
SIGMA = ("sleeve_sigma_fibre", "sleeve_sigma_transverse")


# ---- 1. Params ---------------------------------------------------------------------------------
def test_defaults_non_conducting():
    p = default_params("sleeved_bridged")
    assert (p.sleeve_sigma_fibre, p.sleeve_sigma_transverse) == (0.0, 0.0)
    assert p.sleeve_conductivity_axial == 0.0


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf"), 2e8])
@pytest.mark.parametrize("name", SIGMA)
def test_sigma_validation(name, bad):
    with pytest.raises(ValueError):
        default_params("sleeved_bridged", **{name: bad})


def test_axial_conductivity_through_params():
    # hoop winding without a composite spec: axial = across-fibre value
    p = default_params("sleeved_bridged", sleeve_sigma_fibre=3e4, sleeve_sigma_transverse=50.0)
    assert p.sleeve_conductivity_axial == pytest.approx(50.0)
    q = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_wind_angle=90.0,
                       sleeve_sigma_fibre=3e4, sleeve_sigma_transverse=50.0)
    assert q.sleeve_conductivity_axial == pytest.approx(3e4)
    r = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_layers="0:3,90:1",
                       sleeve_sigma_fibre=3e4, sleeve_sigma_transverse=50.0)
    assert r.sleeve_conductivity_axial == pytest.approx(0.75 * 50.0 + 0.25 * 3e4)


# ---- 2. axial-conductivity helper ---------------------------------------------------------------
def test_ply_axial_conductivity_limits_and_tensor():
    sl, st = 4e4, 30.0
    assert c.ply_axial_conductivity(sl, st, 0.0) == pytest.approx(st)
    assert c.ply_axial_conductivity(sl, st, 90.0) == pytest.approx(sl)
    assert c.ply_axial_conductivity(sl, st, 45.0) == pytest.approx(0.5 * (sl + st))
    # independent check: rotate the ply tensor diag(sl, st) in the hoop-axial plane
    for ang in (0.0, 10.0, 33.0, 60.0, 90.0):
        t = math.radians(ang)
        R = np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])   # columns = ply axes
        S = R @ np.diag([sl, st]) @ R.T
        assert c.ply_axial_conductivity(sl, st, ang) == pytest.approx(S[1, 1], rel=1e-12)
        # the +theta / -theta off-diagonal terms cancel
        Rm = np.array([[math.cos(-t), -math.sin(-t)], [math.sin(-t), math.cos(-t)]])
        Sm = Rm @ np.diag([sl, st]) @ Rm.T
        assert 0.5 * (S[0, 1] + Sm[0, 1]) == pytest.approx(0.0, abs=1e-9 * sl)


def test_axial_conductivity_layers():
    sl, st = 2e4, 100.0
    assert c.axial_conductivity(sl, st, [0.0]) == pytest.approx(st)
    assert c.axial_conductivity(sl, st, [90.0]) == pytest.approx(sl)
    # thickness-weighted (layers in parallel), weights normalised
    v = c.axial_conductivity(sl, st, [0.0, 30.0, 90.0], [2.0, 1.0, 1.0])
    ref = 0.5 * st + 0.25 * (sl * 0.25 + st * 0.75) + 0.25 * sl
    assert v == pytest.approx(ref)
    assert c.axial_conductivity(sl, st, [0.0, 90.0]) == pytest.approx(0.5 * (sl + st))   # equal shares
    assert c.axial_conductivity(0.0, 0.0, [0.0, 45.0], [1, 1]) == 0.0


@pytest.mark.parametrize("args", [(-1.0, 0.0, [0.0]), (1.0, 0.0, [91.0]), (1.0, 0.0, []),
                                  (1.0, 0.0, [0.0, 10.0], [1.0]), (1.0, 0.0, [0.0], [0.0])])
def test_axial_conductivity_rejects(args):
    with pytest.raises(ValueError):
        c.axial_conductivity(*args)


# ---- 3. cache keys -----------------------------------------------------------------------------
# em_model_key of these points computed with the code of HEAD before this change (git archive)
OLD_EM_KEYS = {
    "bridged": "853349de788e", "bridgeless": "8a34fa043398", "bridged_t2_g120": "46e4a985684e",
    "bridged_cf30": "853349de788e", "bridged_mode": "fc9cbcffef2b", "bridgeless_mode": "97ad8771d471",
}


def test_em_model_key_unchanged_for_non_conducting_sleeves():
    d = default_params
    assert em_model_key(d("sleeved_bridged")) == OLD_EM_KEYS["bridged"]
    assert em_model_key(d("sleeved_bridgeless")) == OLD_EM_KEYS["bridgeless"]
    assert em_model_key(d("sleeved_bridged", sleeve_thickness=2.0, gamma=120.0)) == OLD_EM_KEYS["bridged_t2_g120"]
    assert em_model_key(d("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_wind_angle=30.0)) == \
        OLD_EM_KEYS["bridged_cf30"]
    assert em_model_key(d("sleeved_bridged"), "torque|span=60.0") == OLD_EM_KEYS["bridged_mode"]
    assert em_model_key(d("sleeved_bridgeless"), "torque|span=60.0") == OLD_EM_KEYS["bridgeless_mode"]
    # ply conductivities that give no axial conductivity (hoop winding, sigma_T = 0): same model
    assert em_model_key(d("sleeved_bridged", sleeve_sigma_fibre=3e4)) == OLD_EM_KEYS["bridged"]


def test_sigma_fields_keys():
    assert set(SIGMA) <= set(KEY_OMITTED_AT_DEFAULT)
    assert not set(SIGMA) & MECH_ONLY            # they are EM fields
    p = default_params("sleeved_bridged")
    q = p.replace(sleeve_sigma_transverse=1e4)
    assert em_model_key(q) != em_model_key(p)
    assert point_key(q, "em") != point_key(p, "em")
    assert point_key(p.replace(sleeve_sigma_fibre=3e4), "em") != point_key(p, "em")   # a different point...
    # ... and back at the defaults the old key returns
    assert point_key(q.replace(sleeve_sigma_transverse=0.0), "em") == point_key(p, "em")
    # the EM key sees the axial conductivity only: equal axial conductivity -> same model
    a = default_params("sleeved_bridged", sleeve_sigma_transverse=1e4)
    b = default_params("sleeved_bridged", sleeve_fibre="carbon_sm", sleeve_wind_angle=90.0, sleeve_sigma_fibre=1e4)
    assert em_model_key(a) == em_model_key(b)
    assert em_model_key(a) != em_model_key(a.replace(sleeve_sigma_transverse=2e4))
    # winding angle changes the EM key of a conducting sleeve (through the axial conductivity) ...
    e1 = b.replace(sleeve_sigma_transverse=10.0)
    assert em_model_key(e1) != em_model_key(e1.replace(sleeve_wind_angle=45.0))
    # ... but not of a non-conducting one
    nc = default_params("sleeved_bridged", sleeve_fibre="carbon_sm")
    assert em_model_key(nc) == em_model_key(nc.replace(sleeve_wind_angle=45.0))


# ---- 4. averaging --------------------------------------------------------------------------------
def test_period_average_synthetic():
    T = 1e-3
    t = np.linspace(0.0, 3 * T, 3 * 90 + 1)
    y = 5.0 + 2.0 * np.sin(2 * np.pi * t / T) + 0.7 * np.cos(4 * np.pi * t / T)
    r = eddy.period_average(t, y, T)
    assert r["mean"] == pytest.approx(5.0, rel=1e-12)
    assert r["previous"] == pytest.approx(5.0, rel=1e-12)
    assert r["rel_change"] < 1e-12
    assert r["window"] == pytest.approx((2 * T, 3 * T))
    # a decaying start-up transient: the last period is converged, the check sees the decay
    y2 = y + 3.0 * np.exp(-t / (0.2 * T))
    r2 = eddy.period_average(t, y2, T)
    exact = 5.0 + 3.0 * 0.2 * (np.exp(-10) - np.exp(-15))
    assert r2["mean"] == pytest.approx(exact, rel=1e-4)
    assert r2["rel_change"] > 1e-5
    # two periods per block; window edges between samples are interpolated
    r3 = eddy.period_average(t[:-7], y[:-7], T, n_periods=1)
    assert r3["mean"] == pytest.approx(5.0, rel=1e-3)
    assert math.isnan(eddy.period_average(t[: 90 + 50], y[: 90 + 50], T)["previous"])


def test_period_average_rejects_short_record():
    t = np.linspace(0, 1, 11)
    with pytest.raises(ValueError):
        eddy.period_average(t, t, 2.0)
    with pytest.raises(ValueError):
        eddy.period_average(t[::-1], t, 0.5)


# ---- 5. thin-shell estimate --------------------------------------------------------------------
def _wave_grid(n=7, f=500.0, nt=60, nth=720, periods=1):
    w = 2 * np.pi * f
    t = np.linspace(0.0, periods / f, periods * nt + 1)
    th = np.arange(nth) * 2 * np.pi / nth
    return w, t, th


def test_shell_loss_potential_travelling_wave_known_answer():
    n, A0, sigma, L = 7, 2e-3, 1e4, 0.13
    w, t, th = _wave_grid(n)
    radii = np.linspace(0.071, 0.0725, 4)
    a = np.array([[A0 * np.sin(n * th - w * tk) for tk in t] for _ in radii])
    area = np.pi * (radii[-1] ** 2 - radii[0] ** 2)
    exact = sigma * L * area * (w * A0) ** 2 / 2
    # backward difference (as Maxwell's backward Euler): amplitude factor 2 sin(w dt / 2) / (w dt), exact
    dt = t[1] - t[0]
    fac = (2 * np.sin(w * dt / 2) / (w * dt)) ** 2
    rb = eddy.shell_loss_from_potential(t, th, a, sigma, radii, L)
    assert rb["loss"] == pytest.approx(exact * fac, rel=1e-9)
    # central differences converge to the exact value
    rc = eddy.shell_loss_from_potential(t, th, a, sigma, radii, L, method="central")
    assert rc["loss"] == pytest.approx(exact, rel=5e-3)
    # linear in sigma, quadratic in the frequency
    assert eddy.shell_loss_from_potential(t, th, a, 2 * sigma, radii, L)["loss"] == pytest.approx(2 * rb["loss"])


def test_shell_loss_net_current_options():
    n, A0, U, sigma, L = 5, 1e-3, 4e-4, 1e4, 0.1
    w, t, th = _wave_grid(n)
    radii = np.array([0.07, 0.0715])
    dt = t[1] - t[0]
    fac = (2 * np.sin(w * dt / 2) / (w * dt)) ** 2
    area = np.pi * (radii[-1] ** 2 - radii[0] ** 2)
    # travelling wave plus a potential uniform over the ring (a net current if not removed)
    a = np.array([[A0 * np.sin(n * th - w * tk) + U * np.sin(w * tk) for tk in t] for _ in radii])
    ring = eddy.shell_loss_from_potential(t, th, a, sigma, radii, L, net_current="ring")["loss"]
    none = eddy.shell_loss_from_potential(t, th, a, sigma, radii, L, net_current="none")["loss"]
    assert ring == pytest.approx(sigma * L * area * (w * A0) ** 2 / 2 * fac, rel=1e-9)
    assert none == pytest.approx(sigma * L * area * (w ** 2) * (A0 ** 2 + U ** 2) / 2 * fac, rel=1e-9)
    # anti-periodic field over a pole pitch (odd n): sector means do not vanish; 'sector' with 2n
    # sectors removes them, 'ring' does not need to
    sec = eddy.shell_loss_from_potential(t, th, a, sigma, radii, L, net_current="sector", sectors=2 * n)["loss"]
    # the sector mean of E = w A0 cos(n th - w t) over half a wavelength is (2/pi) w A0 |sin(w t)|:
    # its time-averaged square (2/pi^2) (w A0)^2 is a fraction 4/pi^2 of the power (the U term goes too)
    assert sec == pytest.approx(ring * (1 - 4 / np.pi ** 2), rel=1e-3)   # discrete sector means
    with pytest.raises(ValueError):
        eddy.shell_loss_from_potential(t, th, a, sigma, radii, L, net_current="sector", sectors=7)


def test_shell_loss_from_br_matches_potential():
    n, B0, sigma, L = 11, 0.05, 2e4, 0.13
    r_in, r_out = 0.0711, 0.0726
    rm = 0.5 * (r_in + r_out)
    w, t, th = _wave_grid(n, nt=64, nth=512)
    tt = t[:-1]                                      # endpoint excluded
    br = np.array([B0 * np.cos(n * th - w * tk) for tk in tt])
    r = eddy.shell_loss_from_br(tt, th, br, sigma, r_in, r_out, L)
    exact = sigma * L * (r_out - r_in) * 2 * np.pi * rm * (w * rm * B0 / n) ** 2 / 2
    assert r["loss"] == pytest.approx(exact, rel=1e-9)
    assert r["spectrum"][0][0] == n and r["spectrum"][0][1] == pytest.approx(500.0)
    # the wave moves forward at w / n relative to the shell: it drags the shell forward, T = P n / w
    assert r["torque"] == pytest.approx(exact * n / w, rel=1e-9)
    # a stator-fixed harmonic seen from a rotor turning at om moves at -om: drag -P / om
    om = w / n
    br_s = np.array([B0 * np.cos(n * (th + om * tk)) for tk in tt])
    rs = eddy.shell_loss_from_br(tt, th, br_s, sigma, r_in, r_out, L)
    assert rs["torque"] == pytest.approx(-rs["loss"] / om, rel=1e-9)
    # same field through the potential (B_r = (1/r) dA/dtheta -> A = r B0 sin(n th - w t) / n)
    radii = np.array([r_in, r_out])
    a = np.array([[rm * B0 / n * np.sin(n * th - w * tk) for tk in t] for _ in radii])
    rc = eddy.shell_loss_from_potential(t, th, a, sigma, radii, L, method="central")
    assert rc["loss"] == pytest.approx(exact * (2 * np.pi * rm * (r_out - r_in))
                                       / (np.pi * (r_out ** 2 - r_in ** 2)), rel=5e-3)


def test_rotor_frame_removes_rotation():
    # field fixed to the rotor, sampled at stator angles while the rotor turns: static after resampling
    om = 2 * np.pi * 50.0
    phi = np.arange(360) * 2 * np.pi / 360
    t = np.linspace(0, 0.004, 21)
    f = np.array([np.cos(4 * (phi - om * tk)) for tk in t])
    out = eddy.rotor_frame(phi, t, f, om, phi)
    assert np.max(np.abs(out - np.cos(4 * phi))) < 2e-3


def test_sleeve_shell_estimate_synthetic_export():
    """A rotor-frame travelling wave sampled like MaxwellSession.sleeve_lines (stationary arcs over one
    pole, odd periodicity, whole run) gives the known loss through the full pipeline."""
    from sipm.validation import sleeve_shell_estimate

    p = default_params("sleeved_bridged")
    n, A0, sigma = 12, 1e-3, 1e4                      # order 3p: anti-periodic over a pole
    fe = p.em_rpm * p.pole_pairs / 60
    w = 2 * np.pi * 6 * fe                            # rotor-frame frequency (6 f_e)
    om = 2 * np.pi * p.em_rpm / 60
    t = np.linspace(0, 180 / 360 / fe, 91)
    phi = np.linspace(0, 2 * np.pi / p.poles, 721)
    radii = np.linspace(p.r_lam, p.r_sleeve_out, 5) * 1e-3
    a = np.array([[A0 * np.sin(n * (phi - om * tk) - w * tk) for tk in t] for _ in radii])
    r = sleeve_shell_estimate(dict(t=t, phi=phi, radii=radii, Flux_Lines=a), p, sigma)
    dt = t[1] - t[0]
    exact = sigma * p.stack_length * 1e-3 * np.pi * (radii[-1] ** 2 - radii[0] ** 2) * (w * A0) ** 2 / 2
    fac = (2 * np.sin(w * dt / 2) / (w * dt)) ** 2
    assert r["seam_mismatch"] < 1e-12
    assert r["loss_ring"] == pytest.approx(exact * fac, rel=2e-3)
    assert r["loss_none"] == pytest.approx(r["loss_ring"], rel=1e-9)
    assert r["loss_central"] == pytest.approx(exact, rel=2e-2)
    with pytest.raises(ValueError):
        sleeve_shell_estimate(dict(t=t, phi=phi, radii=radii, Flux_Lines=a), p)   # p does not conduct


def test_periods_and_skin_depth():
    assert eddy.loss_period_deg(48, 8) == 60.0
    assert eddy.slot_period_deg(48, 8) == 30.0
    assert eddy.slot_frequency(48, 3000.0) == 2400.0
    with pytest.raises(ValueError):
        eddy.loss_period_deg(36, 8)
    # copper at 50 Hz: about 9.3 mm
    assert eddy.skin_depth(5.8e7, 50.0) == pytest.approx(9.35e-3, rel=1e-2)
    assert eddy.skin_depth(0.0, 50.0) == math.inf


# ---- 6. run length, resolution and the conducting run (fake AEDT) ------------------------------
def test_em_run_span_and_resolution():
    p = default_params("sleeved_bridged")
    assert mx.em_run_span(p) == p.em_span                       # unchanged for a non-conducting sleeve
    q = p.replace(sleeve_sigma_transverse=1e4)
    assert mx.em_run_span(q) == 180.0
    assert mx.em_run_span(q.replace(em_span=150.0)) == 240.0
    assert mx.sleeve_mesh_length(p) == p.mesh_rotor
    assert mx.sleeve_mesh_length(q) == q.mesh_rotor          # same mesh as the non-conducting model
    r = mx.sleeve_resolution(q)
    assert r["sleeve_steps_per_slot"] == pytest.approx(30.0 / q.em_step)
    assert r["sleeve_skin_depth_mm"] == pytest.approx(eddy.skin_depth(1e4, 2400.0) * 1e3)
    mx.check_sleeve_resolution(q)
    with pytest.raises(ValueError, match="em_step"):
        mx.check_sleeve_resolution(q.replace(em_step=3.0))
    with pytest.raises(ValueError, match="skin depth"):
        mx.check_sleeve_resolution(q.replace(sleeve_sigma_transverse=9e7, em_rpm=20000.0))


class _FakeSession(mx.MaxwellSession):
    """MaxwellSession with the AEDT calls replaced by a synthetic run."""

    def __init__(self, loss_mean=5.0, torque=330.0):
        super().__init__()
        self.loss_mean, self.torque, self.built = loss_mean, torque, []

    def _new_project(self, p, mode, sleeve_mesh=None):
        self.built.append(mode)
        self.p = p
        self.key = mx.em_model_key(p, mode)

    def _set_operating_point(self, p, rpm=None):
        pass

    def _analyze(self):
        pass

    def report(self, expressions, name="sipm_q"):
        p = self.p
        fe = p.em_rpm * p.pole_pairs / 60
        span = mx.em_run_span(p)
        t = np.linspace(0, span / 360 / fe, round(span / p.em_step) + 1)
        x = 2 * np.pi * 6 * fe * t
        d = {"Time": t, "Moving1.Torque": self.torque + 10 * np.cos(x)}
        if "SolidLoss" in expressions:
            d["SolidLoss"] = self.loss_mean * (1 + 0.2 * np.sin(x))
            d["SolidLoss"][0] = 0.0          # Maxwell's t = 0 sample: static solution, no eddy current
        return d


def test_conducting_run_outputs():
    s = _FakeSession(loss_mean=5.0, torque=330.0)
    p = default_params("sleeved_bridged", sleeve_sigma_transverse=1e4)
    r = s.solve(p)
    assert s.built == ["torque|conducting"]
    assert r["sleeve_loss_w"] == pytest.approx(5.0, rel=1e-9)
    assert r["sleeve_loss_rel_change"] < 1e-9           # the t = 0 sample lies in the lead-in
    assert r["torque_avg"] == pytest.approx(330.0, rel=1e-9)
    assert r["sleeve_loss_frac"] == pytest.approx(5.0 / (330.0 * 2 * np.pi * 3000 / 60))
    assert r["sleeve_sigma_axial"] == 1e4 and r["sleeve_mesh_mm"] == 1.2
    assert len(r["_wave"]["Time"]) == 91 and "SolidLoss" in r["_wave"]
    with pytest.raises(ValueError):
        s.solve(p, span=60.0)
    with pytest.raises(ValueError, match="em_step"):
        s.solve(p.replace(em_step=4.0))


def test_non_conducting_run_unchanged_and_zero_axial():
    s = _FakeSession()
    p = default_params("sleeved_bridged")
    r = s.solve(p)
    assert s.built == ["torque|span=60.0"]
    assert set(r) == {"torque_avg", "torque_max", "torque_min", "torque_ripple", "em_wall_s", "_wave"}
    assert len(r["_wave"]["Time"]) == 31
    # ply conductivities with zero axial conductivity: same model, loss reported as exactly 0
    q = p.replace(sleeve_sigma_fibre=3e4)
    s2 = _FakeSession()
    r2 = s2.solve(q)
    assert s2.built == ["torque|span=60.0"]
    assert r2["sleeve_loss_w"] == 0.0 and r2["sleeve_sigma_axial"] == 0.0
    assert r2["torque_avg"] == pytest.approx(r["torque_avg"])


# ---- 7. line report parsing ---------------------------------------------------------------------
def test_parse_line_report_csv(tmp_path):
    f = tmp_path / "line.csv"
    f.write_text('"Distance [cm]","Flux_Lines [] - Time=\'27777.7777778ns\'","Flux_Lines [] - Time=\'0ns\'"\n'
                 "0,-0.003,-0.004\n0.5,0.001,0.002\n", encoding="utf-8")
    dist, t, v = mx.parse_line_report_csv(f)
    assert dist == pytest.approx([0.0, 0.005])
    assert t == pytest.approx([0.0, 27777.7777778e-9])          # sorted by time
    assert v[0] == pytest.approx([-0.004, 0.002]) and v[1] == pytest.approx([-0.003, 0.001])
    f.write_text('"Angle [deg]","x [] - Time=\'0ns\'"\n0,1\n', encoding="utf-8")
    with pytest.raises(ValueError):
        mx.parse_line_report_csv(f)


# ---- 8. CLI and configuration ------------------------------------------------------------------
def test_cli_composite_prints_axial_conductivity(capsys):
    assert main(["composite", "--fibre", "carbon_sm", "--angle", "30", "--sigma-fibre", "2e4",
                 "--sigma-transverse", "100", "--toml"]) == 0
    out = capsys.readouterr().out
    ref = 2e4 * 0.25 + 100 * 0.75
    assert f"sleeve_conductivity_axial = {ref:.6g} S/m" in out
    assert "sleeve_sigma_fibre = 20000.0" in out
    assert main(["composite", "--fibre", "carbon_sm", "--layers", "0:1,90:1", "--sigma-fibre", "2e4"]) == 0
    out = capsys.readouterr().out
    assert "axial conductivity of each layer (S/m): 0, 20000" in out
    assert f"= {1e4:.6g} S/m" in out
    assert main(["composite", "--fibre", "carbon_sm"]) == 0
    assert "conductivity" not in capsys.readouterr().out


def test_example_config_resolves():
    cfg = load_config(ROOT / "examples" / "configs" / "em_sleeve_loss.toml")
    sw = sweep_from_config(cfg)
    plan = sw.plan()
    assert len(plan) == 6 and all(err is None for *_, err in plan)
    sig = sorted({round(p.sleeve_conductivity_axial, 6) for _, _, p, _, _ in plan})
    assert sig == [100.0, round(2e4 * 0.25 + 100 * 0.75, 6), round(2e4 * 0.75 + 100 * 0.25, 6)]


def test_post_labels():
    from sipm.post import LABELS

    for k in ("sleeve_loss_w", "sleeve_loss_frac", "sleeve_sigma_axial", "sleeve_sigma_fibre",
              "sleeve_sigma_transverse"):
        assert k in LABELS
