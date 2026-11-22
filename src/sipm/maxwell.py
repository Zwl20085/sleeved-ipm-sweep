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


def em_model_key(p: Params, extra: str = "") -> str:
    d = {k: v for k, v in p.to_dict().items() if k not in MECH_ONLY | EM_VARIABLES}
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
    def _new_project(self, p: Params, mode: str):
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
        build_model(self.app, p, mode)
        self.key = key

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
        """Torque run at (gamma, current_rms). Returns waveform and average torque."""
        t0 = time.time()
        mode = f"torque|span={span or p.em_span}"
        if self.key != em_model_key(p, mode):
            self._new_project(p, mode)
            if span is not None:
                _setup_times(self.app, p, span=span)
        self._set_operating_point(p)
        self._analyze()
        data = self.report(list(quantities))
        tq = data["Moving1.Torque"]
        t = data["Time"]
        n = len(t) - 1  # last point duplicates the first of the next period
        avg = float(np.mean(tq[:n]))
        ripple = float((np.max(tq) - np.min(tq)) / abs(avg)) if avg else float("nan")
        res = dict(torque_avg=avg, torque_max=float(np.max(tq)), torque_min=float(np.min(tq)),
                   torque_ripple=ripple, em_wall_s=time.time() - t0)
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


def build_model(app, p: Params, mode: str = "torque"):
    """Create geometry, materials, excitations, boundaries, mesh and setup in ``app``."""
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
    slv.conductivity = 0
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

    # ---- motion ----------------------------------------------------------------------------
    app.assign_rotate_motion("Band", coordinate_system="Global", axis="Z", positive_movement=True,
                             start_position=0, has_rotation_limits=False, angular_velocity=f"{p.em_rpm}rpm")

    # ---- mesh ------------------------------------------------------------------------------
    rotor_objs = ["Rotor"] + [f"PM{mg.layer}{'u' if mg.side > 0 else 'l'}" for mg in rp.magnets] + ["Sleeve"]
    app.mesh.assign_length_mesh(rotor_objs, inside_selection=True, maximum_length=p.mesh_rotor,
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
    _setup_times(app, p)
    return dict(axes=axes, layout=layout, r_band=r_band)


def phase_axes(p: Params) -> dict:
    sp = build_stator(p)
    return winding_layout(p, sp.slot_angles)[1]
