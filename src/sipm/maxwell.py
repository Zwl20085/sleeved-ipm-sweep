"""Maxwell 2D transient driver (PyAEDT).

One pole is modelled (sector 0 .. 360/poles deg) with odd periodicity
(master/slave, B_slave = -B_master), model depth = stack length and symmetry
multiplier = poles, so torque, flux linkage and loss are reported for the whole
machine. The rotor objects sit inside a rotating band; the band boundary is in
the middle of the physical air gap (sleeve OD to stator bore).

Excitation: three current-driven windings, i_k(t) = Imax cos(w t + theta_d0 + gamma - phi_k),
where theta_d0 is the electrical position of the rotor d-axis at t = 0, phi_k the
magnetic axis of phase k computed from the winding layout and gamma the current
advance angle from the +d axis (90 deg = pure q-axis current). The model is
rotated by 360/(2*poles) deg so the d-axis of the modelled pole lies on the
sector bisector; theta_d0 and phi_k are evaluated in that same frame.

Only ``gamma`` and ``Imax`` are design variables; any other parameter change
creates a new project (geometry, mesh and setup depend on it).

The sleeve is non-conducting unless its axial conductivity (``Params.sleeve_conductivity_axial``)
is > 0; then it is a passive solid conductor with eddy effects and :meth:`MaxwellSession.solve`
also reports the 2D sleeve eddy-current loss (docs/methods.md, 2.1).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import time
from pathlib import Path

import numpy as np

from .config import Config, sanitize_process_environment
from .eddy import loss_period_deg, period_average, skin_depth, slot_frequency, slot_period_deg
from .geometry import Loop, build_rotor, build_stator, winding_layout
from .geometry.stator import slot_profile
from .materials import load_bh, magnet_coercivity
from .params import Params

# Parameters that do not affect the EM model (mechanical-only) or are design variables.
MECH_ONLY = {
    "lam_density", "lam_E", "lam_nu", "pm_density", "pm_E", "pm_nu", "sleeve_density", "sleeve_E_hoop",
    "sleeve_E_radial", "sleeve_E_axial", "sleeve_nu_hoop_radial", "sleeve_nu_hoop_axial",
    "sleeve_nu_radial_axial", "sleeve_G_hoop_radial", "sleeve_G_hoop_axial", "sleeve_G_radial_axial",
    "mech_rpm", "interference", "mu_sleeve", "mu_magnet", "plane", "mech_mesh", "mech_mesh_fine",
    "magnet_contact", "stress_corner_margin",
    "lam_cte", "pm_cte", "sleeve_cte_hoop", "sleeve_cte_radial", "sleeve_cte_axial", "mech_temperature",
    "mech_ref_temperature",
    "sleeve_fibre", "sleeve_matrix", "sleeve_vf", "sleeve_micromechanics", "sleeve_wind_angle",
    "sleeve_Xt", "sleeve_Xc", "sleeve_Yt", "sleeve_Yc", "sleeve_S",
    "mech_cases", "overspeed_factor", "mech_temperature_cold",
    "sleeve_layers", "sleeve_layer_prestress",
}
EM_VARIABLES = {"gamma", "current_rms"}
# Sleeve ply conductivities: the EM model sees them (and the winding angle and layer spec, which are
# in MECH_ONLY) only through the derived axial conductivity, which enters the key instead - and only
# when it is > 0, so the key of every non-conducting sleeve is the key of earlier versions.
SLEEVE_SIGMA_FIELDS = {"sleeve_sigma_fibre", "sleeve_sigma_transverse"}


def em_model_key(p: Params, extra: str = "") -> str:
    """Hash of everything that defines the Maxwell model (not the design variables gamma and
    current_rms). A conducting sleeve adds its axial conductivity (12 significant digits); a
    non-conducting sleeve (axial conductivity 0, whatever the ply values) hashes as before."""
    d = {k: v for k, v in p.to_dict().items() if k not in MECH_ONLY | EM_VARIABLES | SLEEVE_SIGMA_FIELDS}
    sig = p.sleeve_conductivity_axial
    if sig > 0:
        d["sleeve_conductivity_axial"] = float(f"{sig:.12g}")
    s = json.dumps(d, sort_keys=True) + extra
    return hashlib.sha1(s.encode()).hexdigest()[:12]


# ======================================================================================
# report CSV parsing (licence free, unit tested)
# ======================================================================================
_SCALE = {
    "s": 1.0, "ms": 1e-3, "us": 1e-6, "ns": 1e-9,
    "newtonmeter": 1.0, "mnewtonmeter": 1e-3, "knewtonmeter": 1e3, "n-m": 1.0, "nm": 1.0,
    "wb": 1.0, "mwb": 1e-3, "uwb": 1e-6,
    "v": 1.0, "mv": 1e-3, "kv": 1e3,
    "a": 1.0, "ma": 1e-3, "ka": 1e3,
    "w": 1.0, "mw": 1e-3, "kw": 1e3,
    "deg": 1.0, "": 1.0,
}


def parse_report_csv(path: str | Path) -> dict[str, np.ndarray]:
    """Read an AEDT report CSV; columns are returned in SI units keyed by quantity name."""
    import csv

    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], [r for r in rows[1:] if r]
    out = {}
    for j, h in enumerate(header):
        h = h.strip()
        name, unit = h, ""
        if h.endswith("]") and "[" in h:
            name = h[:h.rindex("[")].strip()
            unit = h[h.rindex("[") + 1:-1].strip()
        # drop the variation suffix AEDT appends, e.g. "Moving1.Torque - gamma='130deg'"
        name = name.split(" - ")[0].strip()
        key = unit.lower().replace(" ", "")
        if key not in _SCALE:
            raise ValueError(f"unknown unit '{unit}' in column '{h}' of {path}")
        out[name] = np.array([float(r[j]) for r in body]) * _SCALE[key]
    return out


_LENGTH = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "meter": 1.0}
_FIELD_UNIT = {"": 1.0, "wb/m": 1.0, "tesla": 1.0, "t": 1.0, "mtesla": 1e-3, "a_per_m2": 1.0, "a/m^2": 1.0,
               "a_per_meter2": 1.0, "a_per_meter_sq": 1.0, "w_per_m3": 1.0, "w/m^3": 1.0, "w_per_meter3": 1.0}


def parse_line_report_csv(path: str | Path):
    """Read an AEDT Fields line report exported with Distance as the x axis and one column per saved
    time ("<qty> [unit] - Time='<value><unit>'"). Returns (distance in m, times in s, values of shape
    (n_t, n_distance) in SI units). ``Flux_Lines`` of Maxwell 2D (A_z, Wb/m) is exported without a
    unit ('[]') and is taken as Wb/m."""
    import csv
    import re

    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.reader(fh) if r]
    header, body = rows[0], np.array([[float(x) for x in r] for r in rows[1:]])
    m = re.fullmatch(r"\s*Distance\s*\[(\w*)\]\s*", header[0])
    if not m or m.group(1).lower() not in _LENGTH:
        raise ValueError(f"{path}: first column must be Distance with a length unit, got {header[0]!r}")
    dist = body[:, 0] * _LENGTH[m.group(1).lower()]
    times, cols = [], []
    for j, h in enumerate(header[1:], start=1):
        mm = re.fullmatch(r"\s*(.+?)\s*\[([^\]]*)\]\s*-\s*Time='([-+0-9.eE]+)\s*([a-z]*)'\s*", h)
        if not mm:
            raise ValueError(f"{path}: cannot parse column header {h!r}")
        unit, tval, tunit = mm.group(2).strip().lower(), float(mm.group(3)), mm.group(4).lower()
        if unit not in _FIELD_UNIT or tunit not in _SCALE:
            raise ValueError(f"{path}: unknown unit in column header {h!r}")
        times.append(tval * _SCALE[tunit])
        cols.append(body[:, j] * _FIELD_UNIT[unit])
    order = np.argsort(times)
    return dist, np.array(times)[order], np.array(cols)[order]


# ======================================================================================
# geometry helpers
# ======================================================================================
def _poly_args(loop: Loop):
    """Points and segment types for PyAEDT create_polyline (closed, 3-point arcs)."""
    pts = [[float(loop.pts[0][0]), float(loop.pts[0][1]), 0.0]]
    seg = []
    n = len(loop)
    for i in range(n):
        q = loop.pts[(i + 1) % n]
        if loop.arcs[i] is None:
            seg.append("Line")
        else:
            m = loop.arc_mid(i)
            pts.append([float(m[0]), float(m[1]), 0.0])
            seg.append("Arc")
        pts.append([float(q[0]), float(q[1]), 0.0])
    return pts, seg


def _sector(r_in, r_out, a0, a1) -> Loop:
    def P(r, a):
        return np.array([r * math.cos(a), r * math.sin(a)])

    if r_in <= 0:
        return Loop(np.array([[0.0, 0.0], P(r_out, a0), P(r_out, a1)]), [None, (0.0, 0.0), None])
    return Loop(np.array([P(r_in, a0), P(r_out, a0), P(r_out, a1), P(r_in, a1)]),
                [None, (0.0, 0.0), None, (0.0, 0.0)])


# ======================================================================================
# session
# ======================================================================================
class MaxwellSession:
    """Owns one AEDT process for a sequence of points; always use as a context manager.

    >>> with MaxwellSession(cfg) as s:
    ...     res = s.solve(params)
    """

    def __init__(self, cfg: Config | None = None):
        self.cfg = cfg or Config()
        self.desktop = None
        self.pid = None
        self.app = None
        self.key = None
        self.project_path = None

    # -- process lifetime ------------------------------------------------------------------
    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def open(self):
        from ansys.aedt.core import Desktop

        sanitize_process_environment()  # PyAEDT launches AEDT with this process's environment
        ver = self.cfg.resolved_aedt_version()
        # AEDT writes batch.log etc. into its working directory: start it inside the workdir
        self.cfg.workdir.mkdir(parents=True, exist_ok=True)
        here = os.getcwd()
        os.chdir(self.cfg.workdir)
        try:
            self.desktop = Desktop(version=ver, non_graphical=self.cfg.aedt_non_graphical, new_desktop=True,
                                   close_on_exit=True)
        finally:
            os.chdir(here)
        self.pid = self.desktop.aedt_process_id

    def close(self):
        try:
            self._close_project()
        finally:
            if self.desktop is not None:
                try:
                    self.desktop.release_desktop(close_projects=True, close_on_exit=True)
                except Exception:
                    pass
                self.desktop = None
            _ensure_dead(self.pid)
            self.pid = None

    def _close_project(self):
        if self.app is not None:
            name = self.app.project_name
            try:
                self.app.close_project(name=name, save=False)
            except Exception:
                pass
            self.app = None
            if self.project_path and not self.cfg.keep_aedt_projects:
                for q in (self.project_path, self.project_path.with_suffix(".aedtresults"),
                          self.project_path.with_suffix(".pyaedt"), Path(str(self.project_path) + ".lock")):
                    if q.is_dir():
                        shutil.rmtree(q, ignore_errors=True)
                    elif q.exists():
                        try:
                            q.unlink()
                        except OSError:
                            pass
        self.key = None

    # -- model -----------------------------------------------------------------------------
    def _new_project(self, p: Params, mode: str, sleeve_mesh: float | None = None):
        from ansys.aedt.core import Maxwell2d

        self._close_project()
        key = em_model_key(p, mode)
        d = self.cfg.workdir / "aedt"
        d.mkdir(parents=True, exist_ok=True)
        stamp = time.strftime("%H%M%S")
        self.project_path = d / f"sipm_{key}_{stamp}.aedt"
        self.app = Maxwell2d(project=str(self.project_path), design="em", solution_type="TransientXY",
                             version=self.cfg.resolved_aedt_version(), new_desktop=False,
                             non_graphical=self.cfg.aedt_non_graphical)
        build_model(self.app, p, mode, sleeve_mesh=sleeve_mesh)
        self.key = key if sleeve_mesh is None else None

    def _set_operating_point(self, p: Params, rpm: float | None = None):
        imax = math.sqrt(2) * p.current_rms
        self.app["Imax"] = f"{imax:.10g}A"
        self.app["gamma"] = f"{p.gamma:.10g}deg"

    def messages(self, severity: int = 2) -> list[str]:
        """AEDT message-manager entries (2 = errors, 1 = warnings and above)."""
        try:
            msgs = self.app.odesktop.GetMessages(self.app.project_name, self.app.design_name, severity)
            return [str(x).strip() for x in msgs if str(x).strip()]
        except Exception as e:  # pragma: no cover - depends on AEDT
            return [f"(could not read AEDT messages: {e})"]

    def _analyze(self):
        ok = self.app.analyze_setup("Setup1", cores=self.cfg.aedt_cores)
        if not ok:
            raise RuntimeError("Maxwell analysis failed: " + " | ".join(self.messages(1)[-8:]))

    # -- public ----------------------------------------------------------------------------
    def solve(self, p: Params, span: float | None = None, quantities=("Moving1.Torque",)) -> dict:
        """Torque run at (gamma, current_rms). Returns waveform and average torque.

        With a conducting sleeve (``p.sleeve_conductivity_axial > 0``) the run also gives the 2D
        sleeve eddy-current loss (:meth:`_solve_conducting`)."""
        if p.sleeve_conductivity_axial > 0:
            if span is not None:
                raise ValueError("solve(span=...) is not supported with a conducting sleeve (the run length "
                                 "follows from em_span and the loss period)")
            return self._solve_conducting(p, quantities)
        t0 = time.time()
        mode = f"torque|span={span or p.em_span}"
        if self.key != em_model_key(p, mode):
            self._new_project(p, mode)
            if span is not None:
                _setup_times(self.app, p, span=span)
        self._set_operating_point(p)
        self._analyze()
        data = self.report(list(quantities))
        res = _torque_summary(data, time.time() - t0)
        if p.sleeve_sigma_fibre or p.sleeve_sigma_transverse:
            # ply conductivities given, but the sleeve does not conduct axially (e.g. hoop winding with
            # sigma_transverse = 0): the 2D model has no sleeve current, the 2D loss is exactly 0
            res.update(sleeve_resolution(p))
            res.update(sleeve_loss_w=0.0, sleeve_loss_rel_change=0.0, sleeve_loss_frac=0.0)
        return res

    def _solve_conducting(self, p: Params, quantities=("Moving1.Torque",), sleeve_mesh: float | None = None) -> dict:
        """Torque and 2D sleeve eddy-current loss with a conducting sleeve (docs/methods.md, 2.1).

        The run covers a lead-in of one loss period (start-up transient) and then
        ``max(2, ceil(em_span / P))`` loss periods P (:func:`sipm.eddy.loss_period_deg`). Torque is
        averaged over the last ``em_span`` degrees as in a torque run; the loss (Maxwell ``SolidLoss``,
        only the sleeve is a solid conductor) over the last loss period, with its relative change from
        the period before as the convergence check."""
        check_sleeve_resolution(p)
        t0 = time.time()
        mode = "torque|conducting"
        if sleeve_mesh is not None or self.key != em_model_key(p, mode):
            self._new_project(p, mode, sleeve_mesh=sleeve_mesh)
        self._set_operating_point(p)
        self._analyze()
        qs = list(dict.fromkeys(list(quantities) + ["SolidLoss"]))
        data = self.report(qs)
        t = data["Time"]
        m = max(1, round(p.em_span / p.em_step))
        tail = {k: v[-m - 1:-1] for k, v in data.items()}   # last em_span degrees, end point excluded
        res = _torque_summary({k: np.append(v, v[:1]) for k, v in tail.items()}, 0.0)
        fe = p.em_rpm * p.pole_pairs / 60.0
        per = loss_period_deg(p.slots, p.poles) / 360.0 / fe
        avg = period_average(t, data["SolidLoss"], per)
        loss = float(avg["mean"])
        shaft = res["torque_avg"] * 2 * math.pi * p.em_rpm / 60.0
        res.update(sleeve_resolution(p))
        if sleeve_mesh is not None:
            res["sleeve_mesh_mm"] = sleeve_mesh
        res.update(sleeve_loss_w=loss, sleeve_loss_rel_change=float(avg["rel_change"]),
                   sleeve_loss_frac=loss / shaft if shaft > 0 else float("nan"),
                   em_wall_s=time.time() - t0)
        res["_wave"] = {k: v.tolist() for k, v in data.items()}
        return res

    def report(self, expressions: list[str], name: str = "sipm_q") -> dict[str, np.ndarray]:
        rs = self.app.post.oreportsetup
        try:
            rs.DeleteReports([name])
        except Exception:
            pass
        rs.CreateReport(name, "Transient", "Rectangular Plot", "Setup1 : Transient", ["Domain:=", "Sweep"],
                        ["Time:=", ["All"]], ["X Component:=", "Time", "Y Component:=", expressions])
        out = self.cfg.workdir / "aedt" / f"{name}_{int(time.time() * 1000)}.csv"
        rs.ExportToFile(name, str(out))
        data = parse_report_csv(out)
        out.unlink()
        return data

    def sleeve_lines(self, p: Params, radii=None, n_points: int = 361, quantities=("Flux_Lines",),
                     sleeve_mesh: float | None = None) -> dict:
        """Field quantities on stationary arcs through the sleeve over a whole run (validation of the
        sleeve loss, docs/validation.md, 13).

        One arc per radius (mm; default five radii from the sleeve bore to its OD) spans the modelled
        pole 0 .. 360/poles deg in the model frame; the run is that of :meth:`solve` (``em_run_span``),
        with fields saved at every step. Returns ``t`` (s), ``phi`` (rad, model frame, n_points values
        including both ends), ``radii`` (m) and per quantity an array (n_r, n_t, n_points) in SI units
        (``Flux_Lines`` = A_z in Wb/m). ``sleeve_mesh`` overrides the sleeve mesh."""
        q = p
        r_in, r_out = q.r_lam, q.r_sleeve_out
        radii = list(np.linspace(r_in, r_out, 5) if radii is None else radii)
        self._close_project()
        self._new_project(q, f"lines|{time.time()}", sleeve_mesh=sleeve_mesh)
        pitch = 2 * math.pi / q.poles
        names = []
        for i, r in enumerate(radii):
            pts = [[r * math.cos(a), r * math.sin(a), 0.0] for a in (0.0, pitch / 2, pitch)]
            pl = self.app.modeler.create_polyline(pts, segment_type="Arc", name=f"SleeveArc{i}", non_model=True)
            if not pl:
                raise RuntimeError(f"could not create arc at r = {r} mm")
            names.append(pl.name)
        s = self.app.setups[0]
        s.props["SaveFieldsType"] = "Every N Steps"
        s.props["N Steps"] = "1"
        s.props["Steps From"] = "0s"
        s.props["Steps To"] = s.props["StopTime"]
        s.update()
        self._set_operating_point(q)
        self._analyze()
        rs = self.app.post.oreportsetup
        out = {qq: [] for qq in quantities}
        t_ref = None
        for nm, r in zip(names, radii):
            for qq in quantities:
                rep = f"line_{nm}_{qq}"
                rs.CreateReport(rep, "Fields", "Rectangular Plot", "Setup1 : Transient",
                                ["Context:=", nm, "PointCount:=", int(n_points)],
                                ["Distance:=", ["All"], "Time:=", ["All"]],
                                ["X Component:=", "Distance", "Y Component:=", [qq]])
                f = self.cfg.workdir / "aedt" / f"{rep}_{int(time.time() * 1000)}.csv"
                rs.ExportToFile(rep, str(f))
                dist, tt, vals = parse_line_report_csv(f)
                f.unlink()
                rs.DeleteReports([rep])
                if t_ref is None:
                    t_ref = tt
                elif len(tt) != len(t_ref) or np.max(np.abs(tt - t_ref)) > 1e-12:
                    raise RuntimeError("line reports returned different time sets")
                out[qq].append(vals)
            phi = dist / (r * 1e-3)
        res = dict(t=t_ref, phi=phi, radii=np.array(radii) * 1e-3)
        res.update({qq: np.array(v) for qq, v in out.items()})
        return res

    def open_circuit(self, p: Params, span: float = 360.0) -> dict:
        """No-load run over ``span`` electrical degrees: flux linkages and back-EMF."""
        q = p.replace(current_rms=0.0)
        mode = f"open|span={span}"
        if self.key != em_model_key(q, mode):
            self._new_project(q, mode)
            _setup_times(self.app, q, span=span)
        self._set_operating_point(q)
        self._analyze()
        ex = (["Moving1.Torque"] + [f"FluxLinkage(Phase{k})" for k in "ABC"]
              + [f"InducedVoltage(Phase{k})" for k in "ABC"])
        return self.report(ex)

    def no_load_flux_density(self, p: Params) -> dict:
        """|B| at no load, rotor d-axis on the middle tooth (t = 0): mid-tooth on the d-axis,
        mid-yoke near the q-axis, and the maximum over the stator and the rotor lamination."""
        q = p.replace(current_rms=0.0)
        self._new_project(q, "fields")
        sp = build_stator(q)
        hp = math.pi / q.poles
        x_b = sp.metrics["slot_bottom_radius"]
        x_w = slot_profile(q)[2][0]  # end of the tooth-tip taper
        r_tooth = 0.5 * (x_w + x_b)
        r_yoke = 0.5 * (x_b + q.stator_od / 2)
        # model frame = pole frame rotated by +hp; d-axis tooth at hp, q-axis at 0 (sector edge)
        pts = {"P_tooth": (r_tooth, hp), "P_yoke": (r_yoke, math.radians(1.0))}
        for nm, (r, th) in pts.items():
            self.app.modeler.create_point([r * math.cos(th), r * math.sin(th), 0.0], name=nm)
        _setup_times(self.app, q, span=q.em_step, step=q.em_step)
        s = self.app.setups[0]
        s.props["SaveFieldsType"] = "Every N Steps"
        s.props["N Steps"] = "1"
        s.props["Steps From"] = "0s"
        s.props["Steps To"] = s.props["StopTime"]
        s.update()
        self._set_operating_point(q)
        self._analyze()
        post = self.app.post
        sol = "Setup1 : Transient"
        out = {}
        for nm in pts:
            out[nm] = float(post.get_scalar_field_value("Mag_B", scalar_function="Value", solution=sol,
                                                        intrinsics={"Time": "0s"}, object_name=nm,
                                                        object_type="point"))
        for obj in ("Stator", "Rotor"):
            out[f"{obj}_max"] = float(post.get_scalar_field_value("Mag_B", scalar_function="Maximum", solution=sol,
                                                                  intrinsics={"Time": "0s"}, object_name=obj,
                                                                  object_type="surface"))
        return dict(b_tooth_mid=out["P_tooth"], b_yoke_mid=out["P_yoke"], b_stator_max=out["Stator_max"],
                    b_rotor_max=out["Rotor_max"], r_tooth=r_tooth, r_yoke=r_yoke)

    def ac_loss(self, p: Params, rpm: float, periods: float = 1.5, steps_per_period: int = 90) -> dict:
        """Solid-conductor run at ``rpm``: AC copper loss factor of the active part (see methods notes)."""
        q = p.replace(em_rpm=rpm)
        mode = f"acloss|{periods}|{steps_per_period}"
        t0 = time.time()
        if self.key != em_model_key(q, mode):
            self._new_project(q, mode)
            _setup_times(self.app, q, span=360.0 * periods, step=360.0 / steps_per_period)
        self._set_operating_point(q)
        self._analyze()
        d = self.report(["SolidLoss", "Moving1.Torque", "Current(PhaseA)"])
        t, loss = d["Time"], d["SolidLoss"]
        fe = rpm * q.pole_pairs / 60.0
        last = t >= t[-1] - 1.0 / fe - 1e-12
        tl, ll = t[last], loss[last]
        p_ac = float(np.trapezoid(ll, tl) / (tl[-1] - tl[0]))
        p_dc = dc_copper_loss(q)
        return dict(rpm=rpm, p_ac=p_ac, p_dc=p_dc, kac=p_ac / p_dc,
                    torque_avg=float(np.mean(d["Moving1.Torque"][last])),
                    em_wall_s=time.time() - t0)


def _torque_summary(data: dict, wall_s: float) -> dict:
    """Average, extremes and ripple of ``Moving1.Torque`` over a record whose last point duplicates
    the first of the next period (that point is left out of the average)."""
    tq = data["Moving1.Torque"]
    n = len(data["Time"]) - 1
    avg = float(np.mean(tq[:n]))
    ripple = float((np.max(tq) - np.min(tq)) / abs(avg)) if avg else float("nan")
    res = dict(torque_avg=avg, torque_max=float(np.max(tq)), torque_min=float(np.min(tq)),
               torque_ripple=ripple, em_wall_s=wall_s)
    res["_wave"] = {k: v.tolist() for k, v in data.items()}
    return res


# Conducting sleeve: resolution rules (docs/methods.md, 2.1 and docs/validation.md, 13)
# Time steps per slot-passing period below which a conducting run is refused. The default em_step of
# 2 deg gives 15 and a loss about 1.3 % below the step-converged value (1 deg: 0.4 %, about second
# order in the step); use em_step = 1 for sub-percent loss (docs/validation.md, 13).
MIN_STEPS_PER_SLOT = 15


def sleeve_mesh_length(p: Params) -> float:
    """Max element length in the sleeve: ``mesh_rotor``, conducting or not, so a conducting model differs
    from the non-conducting one only by the conductivity (a 0.3 mm sleeve mesh changed the loss by
    0.03 % but the torque by 0.07 % at the check point, docs/validation.md, 13)."""
    return p.mesh_rotor


def em_run_span(p: Params) -> float:
    """Electrical degrees simulated by an EM run: ``em_span`` for a non-conducting sleeve; with a
    conducting one a lead-in of one loss period plus ``max(2, ceil(em_span / P))`` loss periods."""
    if p.sleeve_conductivity_axial <= 0:
        return p.em_span
    per = loss_period_deg(p.slots, p.poles)
    return per * (1 + max(2, math.ceil(p.em_span / per - 1e-9)))


def sleeve_resolution(p: Params) -> dict:
    """Skin depth, mesh and time-step figures of a conducting sleeve at the run speed ``em_rpm``.

    The slot-passing frequency seen from the rotor is f_s = slots x rpm / 60 (the lowest frequency of
    the slot harmonics; the winding harmonics in the rotor frame are at 6 f_e = f_s / q). Columns:
    ``sleeve_sigma_axial`` (S/m), ``sleeve_skin_depth_mm`` at f_s (inf if not conducting),
    ``sleeve_skin_ratio`` = skin depth / sleeve thickness (>> 1: field uniform through the sleeve),
    ``sleeve_mesh_mm`` (max element length in the sleeve) and ``sleeve_steps_per_slot`` (time steps per
    slot-passing period = slot pitch in electrical degrees / ``em_step``, independent of speed)."""
    sig = p.sleeve_conductivity_axial
    delta = skin_depth(sig, slot_frequency(p.slots, p.em_rpm)) * 1e3
    return dict(sleeve_sigma_axial=sig, sleeve_skin_depth_mm=delta, sleeve_skin_ratio=delta / p.sleeve_thickness,
                sleeve_mesh_mm=sleeve_mesh_length(p),
                sleeve_steps_per_slot=slot_period_deg(p.slots, p.poles) / p.em_step)


def check_sleeve_resolution(p: Params) -> None:
    """Refuse a conducting-sleeve run whose time step or mesh cannot resolve the sleeve currents."""
    r = sleeve_resolution(p)
    err = []
    if r["sleeve_steps_per_slot"] < MIN_STEPS_PER_SLOT - 1e-9:
        err.append(f"em_step = {p.em_step:g} deg gives {r['sleeve_steps_per_slot']:.3g} time steps per slot-passing "
                   f"period, fewer than {MIN_STEPS_PER_SLOT}: use em_step <= "
                   f"{slot_period_deg(p.slots, p.poles) / MIN_STEPS_PER_SLOT:g} deg for the sleeve loss")
    if r["sleeve_mesh_mm"] > r["sleeve_skin_depth_mm"] / 2:
        err.append(f"sleeve mesh {r['sleeve_mesh_mm']:g} mm is coarser than half the skin depth "
                   f"({r['sleeve_skin_depth_mm']:.3g} mm at the slot-passing frequency): reduce mesh_rotor")
    if err:
        raise ValueError("conducting sleeve: " + "; ".join(err))


def dc_copper_loss(p: Params) -> float:
    """DC loss of the active length for the whole machine at current_rms (W)."""
    a_c = p.cond_width * p.cond_height * 1e-6
    r_c = p.stack_length * 1e-3 / (p.cu_conductivity * a_c)
    i_c = p.current_rms / p.parallel_paths
    return p.slots * p.hairpin_layers * i_c * i_c * r_c


def _ensure_dead(pid, timeout=30.0):
    if not pid:
        return
    try:
        import psutil
    except ImportError:
        return
    t0 = time.time()
    while time.time() - t0 < timeout:
        if not psutil.pid_exists(pid):
            return
        time.sleep(1.0)
    try:
        proc = psutil.Process(pid)
        if "ansysedt" in proc.name().lower():  # only the process this session started
            for c in proc.children(recursive=True):
                c.kill()
            proc.kill()
    except Exception:
        pass


def _setup_times(app, p: Params, span: float | None = None, step: float | None = None):
    fe = p.em_rpm * p.pole_pairs / 60.0
    span = p.em_span if span is None else span
    step = p.em_step if step is None else step
    nsteps = max(1, round(span / step))
    t_stop = span / 360.0 / fe
    dt = t_stop / nsteps
    s = app.setups[0]
    s.props["StopTime"] = f"{t_stop:.12g}s"
    s.props["TimeStep"] = f"{dt:.12g}s"
    s.update()


def build_model(app, p: Params, mode: str = "torque", sleeve_mesh: float | None = None):
    """Create geometry, materials, excitations, boundaries, mesh and setup in ``app``.

    ``sleeve_mesh`` overrides the max element length in the sleeve (validation only; default
    :func:`sleeve_mesh_length`)."""
    solid = mode.startswith("acloss")
    if solid:
        # The Maxwell model's coil region is one stranded region per slot; the strand geometry of
        # its AC-loss designs is not determinable from the projects (docs/methods.md).
        raise NotImplementedError("AC-loss (solid conductor) runs need the strand geometry, which is not "
                                  "defined for the rebased stator yet")
    rp, sp = build_rotor(p), build_stator(p)
    layout, axes = winding_layout(p, sp.slot_angles)
    hp_deg = 180.0 / p.poles
    rot = hp_deg  # pole frame -> sector 0 .. 2*hp
    pitch = math.radians(2 * hp_deg)
    m = app.modeler
    m.model_units = "mm"

    # ---- materials -----------------------------------------------------------------------
    mats = app.materials
    bh = load_bh(p.lam_bh)
    if bh is None:
        lam_name = p.lam_bh.split(":", 1)[1]
    else:
        lam_name = "sipm_lamination"
        lam = mats.add_material(lam_name)
        lam.permeability = [[float(h), float(b)] for h, b in bh]
        lam.conductivity = 0
    pm = mats.add_material("sipm_magnet")
    pm.permeability = p.pm_mur
    pm.conductivity = 0
    pm.set_magnetic_coercivity(value=-magnet_coercivity(p.pm_br, p.pm_mur), x=1, y=0, z=0)
    slv = mats.add_material("sipm_sleeve")
    slv.permeability = 1.0
    sigma_slv = p.sleeve_conductivity_axial   # 0 = non-conducting sleeve (the model of earlier versions)
    slv.conductivity = sigma_slv if sigma_slv > 0 else 0
    cu = mats.add_material("sipm_copper")
    cu.permeability = 1.0
    cu.conductivity = p.cu_conductivity

    def poly(loop: Loop, name, material):
        pts, seg = _poly_args(loop.transformed(rot_deg=rot))
        o = m.create_polyline(pts, segment_type=seg, cover_surface=True, close_surface=True, name=name,
                              material=material)
        if not o:
            raise RuntimeError(f"AEDT could not create object {name}")
        o.solve_inside = True  # covered polylines are created with SolveInside off
        return o

    # ---- rotor ---------------------------------------------------------------------------
    # lamination pieces (outer outline minus its holes); several pieces are united into "Rotor"
    pieces = []
    for j, face in enumerate(rp.lamination):
        outer, holes = face.solver_loops()
        name = "Rotor" if j == 0 else f"Rotor_{j + 1}"
        poly(outer, name, lam_name)
        if holes:
            tools = [poly(h, f"Pocket{j + 1}_{i + 1}", "vacuum").name for i, h in enumerate(holes)]
            m.subtract([name], tools, keep_originals=False)
        pieces.append(name)
    if len(pieces) > 1:
        m.unite(pieces)
    for mg in rp.magnets:
        name = f"PM{mg.layer}{'u' if mg.side > 0 else 'l'}"
        poly(mg.loop, name, "sipm_magnet")
    poly(rp.sleeve, "Sleeve", "sipm_sleeve")
    r_band = 0.5 * (rp.r_sleeve_out + p.r_stator_in)
    band = _sector(0.0, r_band, -math.radians(hp_deg), math.radians(hp_deg))
    poly(band, "Band", "vacuum")
    poly(_sector(rp.r_sleeve_out, r_band, -math.radians(hp_deg), math.radians(hp_deg)), "GapRotor", "vacuum")
    poly(_sector(r_band, p.r_stator_in, -math.radians(hp_deg), math.radians(hp_deg)), "GapStator", "vacuum")
    # magnet orientation: a coordinate system per magnet with x along the magnetisation
    for mg in rp.magnets:
        name = f"PM{mg.layer}{'u' if mg.side > 0 else 'l'}"
        c = mg.centre
        th = math.radians(rot)
        R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        cc, mm = R @ c, R @ mg.magnetisation
        m.create_coordinate_system(origin=[float(cc[0]), float(cc[1]), 0.0], name=f"CS_{name}", mode="axis",
                                   x_pointing=[float(mm[0]), float(mm[1]), 0.0],
                                   y_pointing=[float(-mm[1]), float(mm[0]), 0.0])
        m.set_working_coordinate_system("Global")
        m[name].part_coordinate_system = f"CS_{name}"

    # ---- stator ----------------------------------------------------------------------------
    poly(sp.loop, "Stator", lam_name)
    cond_names = []
    for c in sp.conductors:  # one coil region per slot carrying c.turns conductors (Maxwell model)
        nm = f"C{c.slot + 1}"
        poly(c.loop, nm, "sipm_copper")
        cond_names.append((nm, c.slot, c.turns))
    poly(_sector(0.0, p.stator_od / 2, -math.radians(hp_deg), math.radians(hp_deg)), "Region", "vacuum")

    # ---- design settings -------------------------------------------------------------------
    app.model_depth = f"{p.stack_length}mm"
    app.change_symmetry_multiplier(p.poles)

    # ---- boundaries ------------------------------------------------------------------------
    rmid = 0.5 * (p.r_stator_in + p.stator_od / 2)
    e_master = m.get_edgeid_from_position([rmid, 0.0, 0.0], "Region")
    e_slave = m.get_edgeid_from_position([rmid * math.cos(pitch), rmid * math.sin(pitch), 0.0], "Region")
    th_mid = pitch / 2
    e_outer = m.get_edgeid_from_position([p.stator_od / 2 * math.cos(th_mid), p.stator_od / 2 * math.sin(th_mid), 0.0],
                                         "Region")
    if not (e_master and e_slave and e_outer):
        raise RuntimeError("could not find Region edges for boundaries")
    app.assign_master_slave(independent=e_master, dependent=e_slave, reverse_master=False, reverse_slave=True,
                            same_as_master=False, boundary="Periodic")
    app.assign_vector_potential(e_outer, vector_value=0, boundary="ZeroA")

    # ---- windings --------------------------------------------------------------------------
    # In the pole frame the rotor d-axis is at 0 (theta_d0 = 0) and the phase axes are axes[k];
    # rotating the whole model does not change their difference, so theta_d0 - phi_k = -axes[k].
    fe = p.em_rpm * p.pole_pairs / 60.0
    app["Imax"] = f"{math.sqrt(2) * p.current_rms:.10g}A"
    app["gamma"] = f"{p.gamma:.10g}deg"
    app["fe"] = f"{fe:.10g}Hz"
    for ph in "ABC":
        expr = f"Imax*cos(2*pi*fe*Time + gamma + ({-axes[ph]:.10g}deg))"
        app.assign_winding(winding_type="Current", is_solid=solid, current=expr,
                           parallel_branches=p.parallel_paths, name=f"Phase{ph}")
    coils = {"A": [], "B": [], "C": []}
    for nm, slot, turns in cond_names:
        ph, s = layout[slot]
        cname = f"Coil_{nm}"
        app.assign_coil(assignment=[nm], conductors_number=turns, polarity="Positive" if s > 0 else "Negative",
                        name=cname)
        coils[ph].append(cname)
    for ph, cl in coils.items():
        app.add_winding_coils(f"Phase{ph}", coils=cl)
    if solid:
        app.eddy_effects_on([nm for nm, _, _ in cond_names], enable_eddy_effects=True,
                            enable_displacement_current=False)
        others = [o for o in m.object_names if o not in {nm for nm, _, _ in cond_names}]
        app.eddy_effects_on(others, enable_eddy_effects=False, enable_displacement_current=False)
    if sigma_slv > 0:
        _sleeve_eddy(app, p)

    # ---- motion ----------------------------------------------------------------------------
    app.assign_rotate_motion("Band", coordinate_system="Global", axis="Z", positive_movement=True,
                             start_position=0, has_rotation_limits=False, angular_velocity=f"{p.em_rpm}rpm")

    # ---- mesh ------------------------------------------------------------------------------
    rotor_objs = ["Rotor"] + [f"PM{mg.layer}{'u' if mg.side > 0 else 'l'}" for mg in rp.magnets]
    if sleeve_mesh is not None:   # validation: own sleeve mesh
        app.mesh.assign_length_mesh(rotor_objs, inside_selection=True, maximum_length=p.mesh_rotor,
                                    maximum_elements=None, name="L_rotor")
        app.mesh.assign_length_mesh(["Sleeve"], inside_selection=True, maximum_length=sleeve_mesh,
                                    maximum_elements=None, name="L_sleeve")
    else:   # the mesh of earlier versions (sleeve_mesh_length = mesh_rotor), conducting sleeve or not
        app.mesh.assign_length_mesh(rotor_objs + ["Sleeve"], inside_selection=True, maximum_length=p.mesh_rotor,
                                    maximum_elements=None, name="L_rotor")
    app.mesh.assign_length_mesh(["Stator"], inside_selection=True, maximum_length=p.mesh_stator,
                                maximum_elements=None, name="L_stator")
    cmax = min(p.mesh_stator, 0.3) if solid else p.mesh_stator
    app.mesh.assign_length_mesh([nm for nm, _, _ in cond_names], inside_selection=True, maximum_length=cmax,
                                maximum_elements=None, name="L_cond")
    app.mesh.assign_length_mesh(["GapRotor", "GapStator"], inside_selection=True, maximum_length=p.mesh_gap,
                                maximum_elements=None, name="L_gap")

    # ---- setup -----------------------------------------------------------------------------
    setup = app.create_setup("Setup1")
    setup.props["NonlinearSolverResidual"] = "1e-05"
    setup.props["SaveFieldsType"] = "None"
    setup.update()
    _setup_times(app, p, span=em_run_span(p))
    return dict(axes=axes, layout=layout, r_band=r_band)


# Validation only: force zero net axial current in the modelled sleeve sector (a solid current
# winding of 0 A). Physically wrong for a ring sleeve under odd periodicity, see _sleeve_eddy.
SLEEVE_SECTOR_ZERO_CURRENT = False


def _sleeve_eddy(app, p: Params):
    """Eddy effects on the conducting sleeve (a passive solid conductor).

    Net current: the sleeve is a full ring whose axial current returns at its ends; with an ideal
    (resistance-free) end return the only constraint is zero total current through the ring. The
    model holds one pole with odd periodicity, so the sleeve current density is anti-periodic over a
    pole pitch and the ring total is zero by symmetry, whatever the sector carries; a sector of one
    pole carries a net current of alternating sign from pole to pole. Maxwell 2D leaves the net
    current of an eddy-enabled object that is not part of a winding free (E_z = -dA/dt, no added
    uniform field; verified: identical loss to a 0 V solid voltage winding, docs/validation.md, 13),
    which is the physical case. Forcing zero net current in the sector (``SLEEVE_SECTOR_ZERO_CURRENT``,
    a 0 A solid current winding) would add a constraint the ring does not have; it is kept for the
    validation only."""
    # not app.eddy_effects_on: PyAEDT only lists objects whose material conducts >= 1e5 S/m
    # (Material.is_conductor) and silently skips a carbon sleeve below that
    data = ["NAME:Data", "Object Name:=", "Sleeve", "Eddy Effect:=", True]
    app.oboundary.SetEddyEffect(["NAME:Eddy Effect Setting", ["NAME:EddyEffectVector", data]])
    if not app.oboundary.GetEddyEffect("Sleeve"):
        raise RuntimeError("Maxwell did not accept eddy effects on the Sleeve")
    if SLEEVE_SECTOR_ZERO_CURRENT:
        app.assign_winding(assignment=["Sleeve"], winding_type="Current", is_solid=True, current=0, name="SleeveSector")


def phase_axes(p: Params) -> dict:
    sp = build_stator(p)
    return winding_layout(p, sp.slot_angles)[1]
