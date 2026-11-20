"""Validation runs: MAPDL path versus closed form, and a mesh-refinement study."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from . import analytic
from .config import Config
from .mapdl import mech_model_ring, solve_mech
from .params import Params, default_params


def ring_params(case: str, **kw) -> Params:
    """Material set for the disc + sleeve check. ``case`` = 'isotropic' or 'orthotropic'."""
    p = default_params("sleeved_bridged", **kw)
    if case == "isotropic":
        # an isotropic aluminium-like ring on a steel disc (two different materials)
        E, nu = 70.0, 0.33
        G = E / (2 * (1 + nu))
        p = p.replace(sleeve_E_hoop=E, sleeve_E_radial=E, sleeve_E_axial=E, sleeve_nu_hoop_radial=nu,
                      sleeve_nu_hoop_axial=nu, sleeve_nu_radial_axial=nu, sleeve_G_hoop_radial=G,
                      sleeve_G_hoop_axial=G, sleeve_G_radial_axial=G, sleeve_density=2700.0)
    elif case != "orthotropic":
        raise ValueError("case must be 'isotropic' or 'orthotropic'")
    return p


def ring_reference(p: Params, a: float, t: float) -> dict:
    w = p.mech_rpm * 2 * math.pi / 60
    disc = analytic.isotropic(p.lam_E * 1e3, p.lam_nu, p.lam_density * 1e-12)
    ring = analytic.orthotropic(p.sleeve_E_radial * 1e3, p.sleeve_E_hoop * 1e3, p.sleeve_nu_hoop_radial,
                                p.sleeve_density * 1e-12)
    return analytic.summary(a, a + t, disc, ring, w, p.interference)


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
                rows=rows, n_nodes=fe["n_nodes"], wall_s=fe["mech_wall_s"], reference=ref)


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
