"""Validation runs: MAPDL path versus closed form, and a mesh-refinement study."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from . import analytic
from .config import Config
from .mapdl import mech_model_ring, solve_mech
from .params import Params, default_params

RING_CASES = ("isotropic", "orthotropic", "thermal")
THERMAL_DELTA_T = 100.0   # K, rotor temperature above the fitting temperature in the 'thermal' case


def ring_params(case: str, **kw) -> Params:
    """Material set for the disc + sleeve check. ``case`` = 'isotropic', 'orthotropic' or
    'thermal' (the orthotropic set, with the whole rotor ``THERMAL_DELTA_T`` above the
    temperature at which the interference is defined, default CTEs: steel 12e-6 /K, sleeve
    0 /K along the fibre and 30e-6 /K across it)."""
    if case not in RING_CASES:
        raise ValueError(f"case must be one of {RING_CASES}")
    if case == "thermal":
        ref = kw.get("mech_ref_temperature", default_params("sleeved_bridged").mech_ref_temperature)
        kw.setdefault("mech_temperature", ref + THERMAL_DELTA_T)
    p = default_params("sleeved_bridged", **kw)
    if case == "isotropic":
        # an isotropic aluminium-like ring on a steel disc (two different materials)
        E, nu = 70.0, 0.33
        G = E / (2 * (1 + nu))
        p = p.replace(sleeve_E_hoop=E, sleeve_E_radial=E, sleeve_E_axial=E, sleeve_nu_hoop_radial=nu,
                      sleeve_nu_hoop_axial=nu, sleeve_nu_radial_axial=nu, sleeve_G_hoop_radial=G,
                      sleeve_G_hoop_axial=G, sleeve_G_radial_axial=G, sleeve_density=2700.0)
    return p


def _ring_bodies(p: Params):
    disc = analytic.isotropic(p.lam_E * 1e3, p.lam_nu, p.lam_density * 1e-12, alpha=p.lam_cte)
    ring = analytic.orthotropic(p.sleeve_E_radial * 1e3, p.sleeve_E_hoop * 1e3, p.sleeve_nu_hoop_radial,
                                p.sleeve_density * 1e-12, alpha_r=p.sleeve_cte_radial, alpha_t=p.sleeve_cte_hoop)
    return disc, ring


def ring_reference(p: Params, a: float, t: float) -> dict:
    """Closed form for the disc + ring with the materials, speed, interference and uniform
    temperature change of ``p`` (plane stress)."""
    w = p.mech_rpm * 2 * math.pi / 60
    disc, ring = _ring_bodies(p)
    return analytic.summary(a, a + t, disc, ring, w, p.interference, dT=p.mech_delta_t)


def ring_liftoff_reference(p: Params, a: float, t: float, temperature: float | None = None) -> float:
    """Closed-form lift-off speed (rpm) of the disc + ring with the materials and interference of
    ``p`` at the uniform ``temperature`` (degC; default ``p.mech_temperature``) against
    ``p.mech_ref_temperature``."""
    temp = p.mech_temperature if temperature is None else temperature
    disc, ring = _ring_bodies(p)
    w = analytic.liftoff_speed(a, a + t, disc, ring, p.interference, dT=temp - p.mech_ref_temperature)
    return w * 60 / (2 * math.pi)


def ring_liftoff_check(radius: float = 71.1, t_sleeve: float = 1.5, interference: float = 0.1,
                       rpm: float = 20000.0, hot: float = 120.0, mesh: float = 0.5, cfg: Config | None = None,
                       workdir: str | Path | None = None) -> dict:
    """FE lift-off speed of the disc + orthotropic sleeve (load cases standstill and rated, cold =
    the fitting temperature and ``hot``) against the closed form. Returns per temperature label
    the FE speed, its basis, the closed-form speed and the relative error, and the FE result."""
    cfg = cfg or Config()
    p = ring_params("orthotropic", interference=interference, mech_rpm=rpm, mech_mesh=mesh, mech_mesh_fine=mesh,
                    mech_temperature=hot, mech_cases="standstill_cold,rated_cold,standstill_hot,rated_hot")
    p = p.replace(mech_temperature_cold=p.mech_ref_temperature)
    wd = Path(workdir) if workdir else cfg.workdir / "validation_ring_liftoff"
    fe = solve_mech(p, cfg, wd, model=mech_model_ring(radius, t_sleeve))
    rows = []
    for label, temp in (("cold", p.mech_temperature_cold), ("hot", p.mech_temperature)):
        cf = ring_liftoff_reference(p, radius, t_sleeve, temp)
        f = fe[f"liftoff_rpm_{label}"]
        rows.append(dict(label=label, temperature=temp, fe=f, basis=fe[f"liftoff_basis_{label}"], closed_form=cf,
                         rel_error=(f - cf) / cf))
    return dict(rows=rows, fe=fe, params=p)


def liftoff_check(p: Params, label: str = "cold", factor: float = 0.95, cfg: Config | None = None,
                  workdir: str | Path | None = None) -> dict:
    """Lift-off speed of a rotor from its standstill, rated and overspeed cases at one temperature
    (``label`` 'cold' or 'hot'), then two direct single-case runs: at that speed (the minimum
    contact pressure should be ~0) and at ``factor`` x that speed (no element should be open)."""
    if label not in ("cold", "hot"):
        raise ValueError("label must be 'cold' or 'hot'")
    cfg = cfg or Config()
    base = Path(workdir) if workdir else cfg.workdir / f"validation_liftoff_{p.topology}_{label}"
    q = p.replace(mech_cases=f"standstill_{label},rated_{label},overspeed_{label}")
    cases = solve_mech(q, cfg, base / "cases")
    rpm = cases[f"liftoff_rpm_{label}"]
    temp = q.case_speed_temperature(f"standstill_{label}")[1]
    out = dict(topology=p.topology, label=label, temperature=temp, liftoff_rpm=rpm,
               basis=cases[f"liftoff_basis_{label}"], cases=cases, direct={})
    if not (math.isfinite(rpm) and rpm > 0):
        return out
    for tag, s in (("at", rpm), ("below", factor * rpm)):
        r = solve_mech(p.replace(mech_cases="", mech_rpm=s, mech_temperature=temp), cfg, base / f"direct_{tag}")
        out["direct"][tag] = dict(rpm=s, contact_p_min=r["contact_p_min"], contact_p_max=r["contact_p_max"],
                                  sleeve_open_fraction=r["sleeve_open_fraction"])
    out["p_min_standstill"] = cases[f"contact_p_min__standstill_{label}"]
    return out


def ring_check(case: str = "orthotropic", radius: float = 71.1, t_sleeve: float = 1.5,
               interference: float = 0.1, rpm: float = 20000.0, mesh: float = 0.5,
               cfg: Config | None = None, workdir: str | Path | None = None) -> dict:
    """Run the disc + sleeve case through the MAPDL driver and compare with the closed form."""
    cfg = cfg or Config()
    p = ring_params(case, interference=interference, mech_rpm=rpm, mech_mesh=mesh,
                    mech_mesh_fine=mesh)
    model = mech_model_ring(radius, t_sleeve)
    wd = Path(workdir) if workdir else cfg.workdir / f"validation_ring_{case}"
    fe = solve_mech(p, cfg, wd, model=model)
    ref = ring_reference(p, radius, t_sleeve)
    pairs = [("disc max von Mises", fe["svm_lam_max"], ref["svm_disc_max"]),
             ("sleeve max hoop stress", fe["sh_sleeve_max"], ref["sh_ring_max"]),
             ("interface pressure (max over contact)", fe["contact_p_max"], ref["p_interface"]),
             ("interface pressure (min over contact)", fe["contact_p_min"], ref["p_interface"])]
    rows = [dict(quantity=q, fe=f, closed_form=c, rel_error=(f - c) / c if c else float("nan"))
            for q, f, c in pairs]
    return dict(case=case, radius=radius, t_sleeve=t_sleeve, interference=interference, rpm=rpm, mesh=mesh,
                delta_t=p.mech_delta_t, rows=rows, n_nodes=fe["n_nodes"], wall_s=fe["mech_wall_s"], reference=ref,
                fe=fe)


# ---------------------------------------------------------------------------------------------
# layered sleeve: disc + N bonded rings
# ---------------------------------------------------------------------------------------------
# The 3-layer check: hoop / +/-30 deg helical / hoop (2:1:1, bore outward) carbon_sm/epoxy, winding
# prestress per layer, 0.1 mm interference, 20000 rpm and 100 K above the fitting temperature.
MULTILAYER_CASE = dict(sleeve_fibre="carbon_sm", sleeve_layers="0:2,30:1,0:1", sleeve_layer_prestress="100,40,100",
                       interference=0.1, mech_rpm=20000.0, mech_temperature=120.0)


def disc_body(p: Params) -> analytic.Body:
    """The lamination as a solid isotropic disc (closed form)."""
    return analytic.isotropic(p.lam_E * 1e3, p.lam_nu, p.lam_density * 1e-12, alpha=p.lam_cte)


def layer_bodies(p: Params) -> list[analytic.Body]:
    """One closed-form ring body per sleeve layer (bore outward), from each layer's constants."""
    return [analytic.orthotropic(c["sleeve_E_radial"] * 1e3, c["sleeve_E_hoop"] * 1e3, c["sleeve_nu_hoop_radial"],
                                 c["sleeve_density"] * 1e-12, alpha_r=c["sleeve_cte_radial"],
                                 alpha_t=c["sleeve_cte_hoop"]) for c in p.sleeve_layer_constants()]


def multilayer_ring_params(**kw) -> Params:
    """Parameters of the 3-layer disc + ring check (:data:`MULTILAYER_CASE`, overridable)."""
    return default_params("sleeved_bridged", **{**MULTILAYER_CASE, **kw})


def multilayer_ring_reference(p: Params, a: float, t: float, temperature: float | None = None,
                              rpm: float | None = None) -> dict:
    """Closed form (:func:`sipm.analytic.rings_summary`) of the disc (radius ``a``) + the sleeve layers of
    ``p`` (total thickness ``t``) with the interference, prestress, speed and temperature of ``p``
    (or the given ones); adds ``radii``."""
    temp = p.mech_temperature if temperature is None else temperature
    w = (p.mech_rpm if rpm is None else rpm) * 2 * math.pi / 60
    from .params import layer_radii

    radii = layer_radii(a, a + t, p.sleeve_layer_fractions())
    out = analytic.rings_summary(radii, disc_body(p), layer_bodies(p), w, p.interference,
                                 dT=temp - p.mech_ref_temperature, prestrain=p.sleeve_layer_prestrain())
    out["radii"] = radii
    return out


def _layer_node_mean(nodes: dict, layer: int, radius: float, comp: str) -> float:
    m = (nodes["layer"] == layer) & (np.abs(nodes["r"] - radius) < 1e-4)
    return float(np.mean(nodes[comp][m])) if m.any() else float("nan")


def multilayer_ring_check(radius: float = 71.1, t_sleeve: float = 1.5, mesh: float = 0.5, cfg: Config | None = None,
                          workdir: str | Path | None = None, liftoff: bool = True, **kw) -> dict:
    """The layered disc + sleeve (:data:`MULTILAYER_CASE`, overridable by ``kw``) through the MAPDL driver
    against the closed form: disc peak von Mises, contact pressure, and per layer the max hoop stress,
    the hoop stress at the layer's bore and OD (node export of that layer) and the radial stress at the
    layer interfaces. With ``liftoff`` also the lift-off speed from the load cases standstill and rated,
    cold (= the fitting temperature) and hot, against :func:`sipm.analytic.rings_liftoff_speed`."""
    cfg = cfg or Config()
    p = multilayer_ring_params(mech_mesh=mesh, mech_mesh_fine=mesh, sleeve_thickness=t_sleeve, **kw)
    if liftoff:
        p = p.replace(mech_cases="standstill_cold,rated_cold,standstill_hot,rated_hot",
                      mech_temperature_cold=p.mech_ref_temperature)
    wd = Path(workdir) if workdir else cfg.workdir / "validation_ring_multilayer"
    fe = solve_mech(p, cfg, wd, model=mech_model_ring(radius, t_sleeve))
    ref = multilayer_ring_reference(p, radius, t_sleeve)
    nodes, radii = fe["_sleeve_nodes"], ref["radii"]
    pairs = [("disc max von Mises", fe["svm_lam_max"], ref["svm_disc_max"]),
             ("interface pressure (max over contact)", fe["contact_p_max"], ref["p_interface"]),
             ("interface pressure (min over contact)", fe["contact_p_min"], ref["p_interface"])]
    for i, lay in enumerate(ref["layers"]):
        k = i + 1
        bore, od = _layer_node_mean(nodes, k, radii[i], "s_hoop"), _layer_node_mean(nodes, k, radii[i + 1], "s_hoop")
        pairs += [(f"layer {k} max hoop stress", fe[f"sh_sleeve_max_L{k}"], lay["sh_max"]),
                  (f"layer {k} hoop stress at its bore", bore, lay["sh_bore"]),
                  (f"layer {k} hoop stress at its OD", od, lay["sh_od"])]
    for i, s in enumerate(ref["sr_interfaces"]):
        fe_sr = 0.5 * (_layer_node_mean(nodes, i + 1, radii[i + 1], "s_radial")
                       + _layer_node_mean(nodes, i + 2, radii[i + 1], "s_radial"))
        pairs.append((f"radial stress at interface {i + 1}/{i + 2}", fe_sr, s))
    rows = [dict(quantity=q, fe=f, closed_form=c, rel_error=(f - c) / c if c else float("nan")) for q, f, c in pairs]
    lift = []
    if liftoff:
        disc, bodies = disc_body(p), layer_bodies(p)
        for label, temp in (("cold", p.mech_temperature_cold), ("hot", p.mech_temperature)):
            w = analytic.rings_liftoff_speed(radii, disc, bodies, p.interference, temp - p.mech_ref_temperature,
                                             p.sleeve_layer_prestrain())
            cf = w * 60 / (2 * math.pi)
            f = fe[f"liftoff_rpm_{label}"]
            lift.append(dict(label=label, temperature=temp, fe=f, basis=fe[f"liftoff_basis_{label}"], closed_form=cf,
                             rel_error=(f - cf) / cf))
    return dict(rows=rows, liftoff=lift, reference=ref, fe=fe, params=p, n_nodes=fe["n_nodes"],
                wall_s=fe["mech_wall_s"])


def winding_prestress(p: Params, tension, sublayers: int = 1) -> dict:
    """Equivalent ``sleeve_layer_prestress`` of a sleeve wound layer by layer onto the rotor with the
    winding stress ``tension`` (MPa, one value or one per layer): :func:`sipm.analytic.sequential_winding`
    on the closed-form model of ``p`` (the lamination as a solid isotropic disc of radius ``r_lam``, the
    sleeve layers of ``p`` as bonded rings). Adds ``spec`` (the value string for
    ``sleeve_layer_prestress``) and ``radii``."""
    radii = p.sleeve_layer_radii()
    out = analytic.sequential_winding(radii, disc_body(p), layer_bodies(p), tension, sublayers)
    out["spec"] = ",".join(f"{v:.6g}" for v in out["prestress"])
    out["radii"] = radii
    return out


def harmonic(theta_deg, y, order=1):
    """Amplitude and phase (deg) of y ~ A cos(order*theta - phase) over a whole period."""
    th = np.radians(np.asarray(theta_deg, float))
    y = np.asarray(y, float)
    c = 2 * np.mean(y * np.cos(order * th))
    s = 2 * np.mean(y * np.sin(order * th))
    return float(np.hypot(c, s)), float(np.degrees(np.arctan2(s, c)))


def open_circuit_check(data: dict, p: Params) -> dict:
    """Checks on a 360-electrical-degree no-load run (``MaxwellSession.open_circuit``).

    * fundamental flux-linkage amplitudes of A, B, C equal (balanced);
    * phase of each flux linkage equals the axis computed from the winding layout
      (psi_k = Psi cos(theta_d - phi_k) with theta_d = 0 at t = 0), which confirms the
      magnet polarity, the coil polarity convention and the current phasing;
    * half-wave symmetry psi(theta + 180) = -psi(theta) (master/slave odd periodicity);
    * no-load torque average ~ 0 (cogging only).
    """
    from .geometry import build_stator, winding_layout

    t = np.asarray(data["Time"])
    fe = p.em_rpm * p.pole_pairs / 60.0
    theta = 360.0 * fe * t
    n = len(theta) - 1
    _, axes = winding_layout(p, build_stator(p).slot_angles)
    out = dict(phases={})
    for k in "ABC":
        psi = np.asarray(data[f"FluxLinkage(Phase{k})"])[:n]
        emf = np.asarray(data[f"InducedVoltage(Phase{k})"])[:n]
        amp, ph = harmonic(theta[:n], psi)
        eamp, eph = harmonic(theta[:n], emf)
        half = n // 2
        sym = float(np.max(np.abs(psi[:half] + psi[half:2 * half])) / amp) if n % 2 == 0 else float("nan")
        thd = math.sqrt(sum(harmonic(theta[:n], emf, h)[0] ** 2 for h in range(2, 40))) / eamp
        dphi = (ph - axes[k] + 180) % 360 - 180
        out["phases"][k] = dict(psi_amp=amp, psi_phase=ph, axis_predicted=axes[k], phase_error=dphi,
                                emf_amp=eamp, emf_phase=eph, emf_thd=thd, half_wave_asym=sym)
    amps = [out["phases"][k]["emf_amp"] for k in "ABC"]
    out["emf_amp_spread"] = (max(amps) - min(amps)) / np.mean(amps)
    out["emf_rms_line_fund"] = float(np.mean(amps) / math.sqrt(2) * math.sqrt(3))
    tq = np.asarray(data["Moving1.Torque"])[:n]
    out["cogging_pk_pk"] = float(tq.max() - tq.min())
    out["cogging_mean"] = float(tq.mean())
    pa, pb, pc = (out["phases"][k]["emf_phase"] for k in "ABC")
    out["phase_shift_AB"] = (pb - pa) % 360
    out["phase_shift_AC"] = (pc - pa) % 360
    return out


def mesh_study(p: Params, sizes=((0.8, 0.3), (0.5, 0.15), (0.3, 0.08), (0.2, 0.04)), cfg: Config | None = None,
               workdir: str | Path | None = None) -> list[dict]:
    """Solve one design with successively finer meshes; returns one row per mesh."""
    cfg = cfg or Config()
    base = Path(workdir) if workdir else cfg.workdir / "mesh_study"
    rows = []
    for h, hf in sizes:
        q = p.replace(mech_mesh=h, mech_mesh_fine=hf)
        r = solve_mech(q, cfg, base / f"h{h}_hf{hf}")
        sections = {k: v for k, v in r.items() if k.startswith("svm_") and k not in (
            "svm_lam_max", "svm_lam_max_all", "svm_pm_max", "svm_lam_r", "svm_lam_theta", "svm_lam_all_r",
            "svm_lam_all_theta", "svm_rotor_max", "svm_sleeve_max")}
        rows.append(dict(mech_mesh=h, mech_mesh_fine=hf, n_nodes=r["n_nodes"],
                         svm_lam_max=r["svm_lam_max"], svm_lam_max_all=r.get("svm_lam_max_all"),
                         svm_pm_max=r.get("svm_pm_max"),
                         sh_sleeve_max=r.get("sh_sleeve_max"), svm_lam_r=r["svm_lam_r"],
                         svm_lam_theta=r["svm_lam_theta"], wall_s=r["mech_wall_s"], **sections))
    return rows
