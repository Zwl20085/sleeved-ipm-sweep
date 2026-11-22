"""Mechanical (centrifugal + interference) analysis with Ansys Mechanical APDL.

The driver writes a self-contained APDL input deck, runs the MAPDL executable in
batch mode and reads a plain ``key = value`` text file that the deck writes in
/POST1. No PyMAPDL dependency: the deck can be inspected and re-run by hand.

Model (see docs/methods.md for the full list of assumptions):

* half a pole (0 <= theta <= pi/poles) in the pole frame, symmetry (zero tangential
  displacement in nodal cylindrical axes) on both radial cut lines; this is exact for
  centrifugal load and interference because the rotor is symmetric about both the
  d- and q-axes;
* 2D plane stress (default) or plane strain, PLANE183 with free triangular mesh;
* lamination isotropic, magnets isotropic, sleeve cylindrically orthotropic
  (element axes x = radial, y = hoop through ESYS of a local cylindrical system); a layered
  sleeve (``sleeve_layers``) is one area and material per layer on shared lines (perfect bond),
  post-processed layer by layer (nodal averages within a layer, not across its interfaces);
* optional winding prestress (``sleeve_layer_prestress``): a hoop eigenstrain per layer as an
  INISTATE initial elastic strain in the element axes, applied with the other loads;
* sleeve/lamination: surface-to-surface contact, interference applied as a contact
  offset (CNOF) ramped over the load step, friction ``mu_sleeve``;
* magnet/pocket (``magnet_contact``): 'bonded' (default, as the original Workbench model):
  bonded MPC contact that also ties the faces across the designed clearance (the gap is
  kept in the geometry and bridged by the bond); 'frictional': surface-to-surface contact
  with friction ``mu_magnet``, faces that touch start closed, faces with the designed
  clearance start open, contact stabilisation damping keeps the loose magnets from drifting
  (not supported for the bridgeless rotor);
* lamination bore traction-free (no shaft modelled); units mm, t, s, N, MPa;
* optional uniform temperature ``mech_temperature`` against the fitted reference
  ``mech_ref_temperature`` (thermal expansion of all three materials, the sleeve's in its
  element axes); the thermal commands are written only when the two differ, so an
  isothermal deck solves exactly as without them;
* headline lamination stress: the peak von Mises stress over named load-carrying sections
  (:mod:`sipm.geometry.locations`), read from a nodal dump of the lamination; the global
  lamination maximum is kept as a secondary, mesh-dependent quantity.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import Config, child_environment
from .geometry import Loop, RotorPole, build_rotor
from .params import Params, layer_radii

RESULTS_FILE = "sipm_results.txt"
NODES_FILE = "sipm_lam_nodes.txt"   # lamination nodes: x, y, von Mises (named sections)
SLEEVE_FILE = "sipm_sleeve_nodes.txt"   # sleeve corner nodes: r, theta, stresses in RSYS 11
# columns of SLEEVE_FILE and keys of parse_sleeve_nodes(): r (mm), theta (deg, pole frame) and
# the nodal-averaged stresses (MPa) in the cylindrical system: SX radial, SY hoop, SZ axial,
# SXY radial-hoop shear
SLEEVE_COLUMNS = ("r", "theta", "s_radial", "s_hoop", "s_axial", "s_shear")
# per-element sleeve contact export (decks written with contact_elements=True): element number,
# centroid x, y (mm, global), contact pressure CONT,PRES (MPa, > 0 compressive; the quantity
# behind contact_p_min) and contact status CONT,STAT (0/1 open, 2 sliding, 3 sticking)
CONTACT_FILE = "sipm_contact_elems.txt"
CONTACT_COLUMNS = ("elem", "x", "y", "pres", "status")
UNSUPPORTED_FRICTIONAL = (
    "magnet_contact='frictional' is not supported for sleeved_bridgeless: the loose magnets in "
    "open-clearance contact do not converge reliably in that rotor (docs/validation.md); use "
    "magnet_contact='bonded'")

MAT_LAM, MAT_PM, MAT_SLV = 1, 2, 3
# layered sleeve (sleeve_layers): layer i (1 = bore) has material MAT_LAYER0 + i and area 300 + i
MAT_LAYER0 = 40
# extra column of SLEEVE_FILE written by a layered deck: the layer (1 = bore) whose elements the
# node's stresses are averaged over; a node on a layer interface appears once per layer
SLEEVE_LAYER_COLUMN = "layer"
INTERFACE_TOL = 1e-4   # mm, nodes within this of a layer interface radius lie on it


@dataclass
class MechModel:
    """Geometry and data the deck writer needs (independent of how it was built)."""

    lamination: list            # [(outer Loop, [hole Loop, ...]), ...] one entry per lamination piece
    magnets: list
    sleeve: Loop | None
    r_lam: float
    sector_deg: float
    origin_node: bool = False   # solid disc: sector has a vertex on the axis
    fine_lines: float = 3.0     # lines shorter than this get the fine element size
    sections: list = field(default_factory=list)   # named stress sections (geometry.locations.Section)
    magnet_clearance: float = 0.0   # mm, designed gap between magnet faces and pocket walls (0 = none)
    meta: dict = field(default_factory=dict)


def _pol(r, th):
    return np.array([r * math.cos(th), r * math.sin(th)])


def mech_model_from_rotor(rp: RotorPole) -> MechModel:
    """Half pole (0 <= theta <= pi/p): lamination pieces with their holes (cavity curves as
    biarcs), the upper magnets, the sleeve and the named stress sections."""
    from .geometry.locations import stress_sections

    hp = rp.half_pitch
    lam = [f.solver_loops() for f in rp.half_lamination()]
    slv = Loop(np.array([_pol(rp.r_lam, 0.0), _pol(rp.r_sleeve_out, 0.0), _pol(rp.r_sleeve_out, hp),
                         _pol(rp.r_lam, hp)]), [None, (0.0, 0.0), None, (0.0, 0.0)])
    return MechModel(lamination=lam, magnets=[m.loop for m in rp.half_magnets()],
                     sleeve=slv, r_lam=rp.r_lam, sector_deg=math.degrees(hp), sections=stress_sections(rp),
                     magnet_clearance=max(rp.p.pm_side_clearance, rp.p.pm_end_clearance))


def mech_model_ring(r_disc: float, t_sleeve: float, sector_deg: float = 22.5) -> MechModel:
    """Solid disc + sleeve (validation case against the closed-form solution)."""
    th = math.radians(sector_deg)
    lam = Loop(np.array([[0.0, 0.0], _pol(r_disc, 0.0), _pol(r_disc, th)]), [None, (0.0, 0.0), None])
    ro = r_disc + t_sleeve
    slv = Loop(np.array([_pol(r_disc, 0.0), _pol(ro, 0.0), _pol(ro, th), _pol(r_disc, th)]),
               [None, (0.0, 0.0), None, (0.0, 0.0)])
    return MechModel(lamination=[(lam, [])], magnets=[], sleeve=slv, r_lam=r_disc,
                     sector_deg=sector_deg, origin_node=True)


# ======================================================================================
# deck writer
# ======================================================================================
def _num(x: float) -> str:
    return f"{x:.10g}"


def _loop_cmds(loop: Loop, kbase: int, lbase: int, area: int, name: str) -> list[str]:
    out = [f"! ---- {name}"]
    n = len(loop)
    if 2 * n >= 1000:
        raise ValueError(f"{name}: {n} segments exceed the keypoint numbering block of 1000")
    for i, (x, y) in enumerate(loop.pts):
        out.append(f"K,{kbase + i + 1},{_num(x)},{_num(y)}")
    centres = {}
    for i, c in enumerate(loop.arcs):
        if c is not None:
            kc = kbase + n + 1 + len(centres)
            centres[i] = kc
            out.append(f"K,{kc},{_num(c[0])},{_num(c[1])}")
    out.append(f"NUMSTR,LINE,{lbase}")
    for i in range(n):
        k1, k2 = kbase + i + 1, kbase + (i + 1) % n + 1
        if loop.arcs[i] is None:
            out.append(f"L,{k1},{k2}")
        else:
            r = float(np.hypot(*(loop.pts[i] - loop.arcs[i])))
            out.append(f"LARC,{k1},{k2},{centres[i]},{_num(r)}")
    out.append(f"LSEL,S,LINE,,{lbase},{lbase + n - 1}")
    out.append(f"NUMSTR,AREA,{area}")
    out.append("AL,ALL")
    out.append("LSEL,ALL")
    return out


def _layer_lines(n: int) -> tuple[list[int], list[int], list[int]]:
    """Line numbers of a layered sleeve: radial edges at theta = 0 and at the sector angle (one per
    layer) and the n + 1 arcs (bore, interfaces, OD)."""
    return ([30000 + i for i in range(n)], [30000 + n + i for i in range(n)],
            [30000 + 2 * n + i for i in range(n + 1)])


def _layered_sleeve_cmds(radii, sector_deg: float) -> list[str]:
    """Areas 301 .. 300 + n of a layered sleeve between the radii (bore .. OD), built on shared
    keypoints and lines so that neighbouring layers share the interface arc, and hence its nodes
    (perfect bond)."""
    n = len(radii) - 1
    th = math.radians(sector_deg)
    out = [f"! ---- sleeve: {n} bonded layers, shared interface lines"]
    for i, r in enumerate(radii):
        out.append(f"K,{30001 + i},{_num(r)},0")
        out.append(f"K,{30101 + i},{_num(r * math.cos(th))},{_num(r * math.sin(th))}")
    out.append("K,30200,0,0")
    out.append("NUMSTR,LINE,30000")
    for i in range(n):
        out.append(f"L,{30001 + i},{30002 + i}")
    for i in range(n):
        out.append(f"L,{30101 + i},{30102 + i}")
    for i, r in enumerate(radii):
        out.append(f"LARC,{30001 + i},{30101 + i},30200,{_num(r)}")
    rad0, rad1, arcs = _layer_lines(n)
    for i in range(n):
        out.append(f"NUMSTR,AREA,{301 + i}")
        out.append(f"AL,{rad0[i]},{arcs[i + 1]},{rad1[i]},{arcs[i]}")
    return out


def _sleeve_material(a, mat: int, c: dict, label: str) -> None:
    """MP commands of a cylindrically orthotropic sleeve material from the thirteen constants ``c``
    (local x = radial, y = hoop, z = axial; PRxy is MAPDL's minor convention nu_xy = -eps_y/eps_x
    for stress along x)."""
    a(f"! {label}: cylindrically orthotropic, x = radial, y = hoop, z = axial")
    a(f"MP,EX,{mat},{_num(c['sleeve_E_radial'] * 1e3)}")
    a(f"MP,EY,{mat},{_num(c['sleeve_E_hoop'] * 1e3)}")
    a(f"MP,EZ,{mat},{_num(c['sleeve_E_axial'] * 1e3)}")
    a(f"MP,PRXY,{mat},{_num(c['sleeve_nu_hoop_radial'] * c['sleeve_E_radial'] / c['sleeve_E_hoop'])}")
    a(f"MP,PRYZ,{mat},{_num(c['sleeve_nu_hoop_axial'])}")
    a(f"MP,PRXZ,{mat},{_num(c['sleeve_nu_radial_axial'])}")
    a(f"MP,GXY,{mat},{_num(c['sleeve_G_hoop_radial'] * 1e3)}")
    a(f"MP,GYZ,{mat},{_num(c['sleeve_G_hoop_axial'] * 1e3)}")
    a(f"MP,GXZ,{mat},{_num(c['sleeve_G_radial_axial'] * 1e3)}")
    a(f"MP,DENS,{mat},{_num(c['sleeve_density'] * 1e-12)}")


def _sleeve_area_size(model: MechModel) -> float:
    """A third of the radial thickness of the single-layer sleeve of ``model`` (from its outline)."""
    rr = np.hypot(model.sleeve.pts[:, 0], model.sleeve.pts[:, 1])
    return float(rr.max() - rr.min()) / 3


def bore_arc_sized(p: Params, model: MechModel) -> bool:
    """True when the deck sizes the sleeve bore arc (the contact target) explicitly to ``mech_mesh``
    for a single-layer sleeve: the sleeve is thinner than 3 x ``mech_mesh``, so its area size
    (thickness / 3) is finer than the lamination's. Layered sleeves size that arc always."""
    if model.sleeve is None or p.is_layered:
        return False
    return _sleeve_area_size(model) < p.mech_mesh


def write_deck(p: Params, model: MechModel, path: str | Path, results_file: str = RESULTS_FILE,
               contact_elements: bool = False) -> Path:
    """Write the APDL input deck for one mechanical point and return its path.

    The deck solves the single load case ``p.mech_rpm`` at ``p.mech_temperature``
    (``p.mech_cases`` is handled by :func:`solve_mech`, one deck per case, see
    :meth:`Params.for_case`). ``contact_elements`` adds the export of the per-element pressure
    and status of the sleeve contact to :data:`CONTACT_FILE` (used for the lift-off speed);
    without it the deck is the same as before that option existed. With ``p.sleeve_layers`` set
    the sleeve is written layer by layer (areas 301.., materials ``MAT_LAYER0 + i``) and with a
    non-zero ``p.sleeve_layer_prestress`` the INISTATE block is added; with neither, the deck is
    the single-layer deck of before."""
    if p.magnet_contact == "frictional" and p.is_bridgeless and model.magnets:
        raise ValueError(UNSUPPORTED_FRICTIONAL)
    layered = model.sleeve is not None and p.is_layered
    if layered:
        layer_consts = p.sleeve_layer_constants()
        nl = len(layer_consts)
        layer_mats = [MAT_LAYER0 + 1 + i for i in range(nl)]
        slv_sel = f"{layer_mats[0]},{layer_mats[-1]}"   # ESEL range of all sleeve elements
        rr_s = np.hypot(model.sleeve.pts[:, 0], model.sleeve.pts[:, 1])
        layer_r = layer_radii(float(rr_s.min()), float(rr_s.max()), p.sleeve_layer_fractions())
    else:
        layer_mats, slv_sel = [MAT_SLV], str(MAT_SLV)
    prestrain = p.sleeve_layer_prestrain() if model.sleeve is not None else ()
    path = Path(path)
    omega = p.mech_rpm * 2 * math.pi / 60.0
    tol = 1e-4
    sector = model.sector_deg
    kop3 = 0 if p.plane == "stress" else 2
    h = p.mech_mesh
    hf = min(p.mech_mesh_fine, h)
    L = []
    a = L.append
    a("/BATCH")
    a("! sleeved-ipm-sweep mechanical deck (generated; units mm, t, s, N, MPa)")
    a("/NOPR")
    a("/NERR,500,10000")
    a("/TITLE,sipm mechanical point")
    a("/PREP7")
    a(f"*SET,R_LAM,{_num(model.r_lam)}")
    a(f"*SET,OMG,{_num(omega)}")
    a(f"*SET,DELTA,{_num(p.interference)}")
    a(f"ET,1,PLANE183,,,{kop3}")
    # materials ---------------------------------------------------------------------------
    a(f"MP,EX,{MAT_LAM},{_num(p.lam_E * 1e3)}")
    a(f"MP,PRXY,{MAT_LAM},{_num(p.lam_nu)}")
    a(f"MP,DENS,{MAT_LAM},{_num(p.lam_density * 1e-12)}")
    a(f"MP,EX,{MAT_PM},{_num(p.pm_E * 1e3)}")
    a(f"MP,PRXY,{MAT_PM},{_num(p.pm_nu)}")
    a(f"MP,DENS,{MAT_PM},{_num(p.pm_density * 1e-12)}")
    if layered:
        angles, fracs = p.sleeve_layer_angles(), p.sleeve_layer_fractions()
        for i, (mat, c) in enumerate(zip(layer_mats, layer_consts)):
            _sleeve_material(a, mat, c, f"sleeve layer {i + 1} of {nl} (+/-{angles[i]:g} deg, "
                                        f"{100 * fracs[i]:.4g} % of the thickness)")
    elif model.sleeve is not None:
        # local x = radial, y = hoop, z = axial. PRxy is the MAPDL "minor" convention:
        # nu_xy = -eps_y/eps_x for stress along x.
        a("! sleeve: cylindrically orthotropic, x = radial, y = hoop, z = axial")
        a(f"MP,EX,{MAT_SLV},{_num(p.sleeve_E_radial * 1e3)}")
        a(f"MP,EY,{MAT_SLV},{_num(p.sleeve_E_hoop * 1e3)}")
        a(f"MP,EZ,{MAT_SLV},{_num(p.sleeve_E_axial * 1e3)}")
        a(f"MP,PRXY,{MAT_SLV},{_num(p.sleeve_nu_radial_hoop)}")
        a(f"MP,PRYZ,{MAT_SLV},{_num(p.sleeve_nu_hoop_axial)}")
        a(f"MP,PRXZ,{MAT_SLV},{_num(p.sleeve_nu_radial_axial)}")
        a(f"MP,GXY,{MAT_SLV},{_num(p.sleeve_G_hoop_radial * 1e3)}")
        a(f"MP,GYZ,{MAT_SLV},{_num(p.sleeve_G_hoop_axial * 1e3)}")
        a(f"MP,GXZ,{MAT_SLV},{_num(p.sleeve_G_radial_axial * 1e3)}")
        a(f"MP,DENS,{MAT_SLV},{_num(p.sleeve_density * 1e-12)}")
    thermal = p.mech_delta_t != 0.0
    if thermal:
        # secant CTEs about TREF = the temperature of the stress-free fitted dimensions (where the
        # interference is defined); the sleeve's in its element system (x radial, y hoop, z axial)
        a(f"! thermal load: uniform {_num(p.mech_temperature)} degC, reference {_num(p.mech_ref_temperature)} degC")
        a(f"TREF,{_num(p.mech_ref_temperature)}")
        for mat, ax, ay, az in ((MAT_LAM, p.lam_cte, p.lam_cte, p.lam_cte),
                                (MAT_PM, p.pm_cte, p.pm_cte, p.pm_cte)):
            a(f"MP,ALPX,{mat},{_num(ax)}")
            a(f"MP,ALPY,{mat},{_num(ay)}")
            a(f"MP,ALPZ,{mat},{_num(az)}")
        if layered:
            for mat, c in zip(layer_mats, layer_consts):
                a(f"MP,ALPX,{mat},{_num(c['sleeve_cte_radial'])}")
                a(f"MP,ALPY,{mat},{_num(c['sleeve_cte_hoop'])}")
                a(f"MP,ALPZ,{mat},{_num(c['sleeve_cte_axial'])}")
        elif model.sleeve is not None:
            a(f"MP,ALPX,{MAT_SLV},{_num(p.sleeve_cte_radial)}")
            a(f"MP,ALPY,{MAT_SLV},{_num(p.sleeve_cte_hoop)}")
            a(f"MP,ALPZ,{MAT_SLV},{_num(p.sleeve_cte_axial)}")
    else:
        a("! isothermal (mech_temperature = mech_ref_temperature): no thermal expansion")
    a("LOCAL,11,1,0,0,0")
    a("CSYS,0")
    # geometry ----------------------------------------------------------------------------
    # lamination pieces: outer area 1.., holes 100.., each piece minus its own holes
    nh = 0
    for j, (outer, holes) in enumerate(model.lamination):
        L += _loop_cmds(outer, 100000 * (j + 1), 100000 * (j + 1), 1 + j, f"lamination piece {j + 1} outline")
        first = nh
        for hole in holes:
            L += _loop_cmds(hole, 1000 * (nh + 1), 1000 * (nh + 1), 100 + nh, f"lamination piece {j + 1} hole")
            nh += 1
        if nh > first:
            a(f"ASEL,S,AREA,,{100 + first},{100 + nh - 1}")
            a(f"CM,CUT{j + 1},AREA")
            a("ASEL,ALL")
            a(f"ASBA,{1 + j},CUT{j + 1}")
    a("ALLSEL")
    a(f"AATT,{MAT_LAM},1,1,0")
    for i, m in enumerate(model.magnets):
        L += _loop_cmds(m, 20000 + 1000 * i, 20000 + 1000 * i, 201 + i, f"magnet {i + 1}")
    if model.magnets:
        a(f"ASEL,S,AREA,,201,{200 + len(model.magnets)}")
        a(f"AATT,{MAT_PM},1,1,0")
    if layered:
        L += _layered_sleeve_cmds(layer_r, sector)
        for i, mat in enumerate(layer_mats):
            a(f"ASEL,S,AREA,,{301 + i}")
            a(f"AATT,{mat},1,1,11")
    elif model.sleeve is not None:
        L += _loop_cmds(model.sleeve, 30000, 30000, 301, "sleeve")
        a("ASEL,S,AREA,,301")
        a(f"AATT,{MAT_SLV},1,1,11")
    a("ALLSEL")
    # mesh --------------------------------------------------------------------------------
    a("MSHAPE,1,2D")
    a("MSHKEY,0")
    a(f"ESIZE,{_num(h)}")
    a(f"LSEL,S,LENGTH,,0,{_num(model.fine_lines)}")
    a(f"LESIZE,ALL,{_num(hf)},,,,,,,1")
    a("LSEL,ALL")
    if layered:
        # at least two quadratic elements through every layer: area size <= half the layer, and the
        # radial edges divided explicitly (they are short lines, sized hf above)
        rad0, rad1, arcs = _layer_lines(nl)
        sizes = []
        for i in range(nl):
            t_i = layer_r[i + 1] - layer_r[i]
            hs = min(h, t_i / 2)
            sizes.append((t_i, hs))
            a(f"ASEL,S,AREA,,{301 + i}")
            a(f"AESIZE,ALL,{_num(hs)}")
        a("ALLSEL")
        for i, (t_i, hs) in enumerate(sizes):
            ndiv = max(2, math.ceil(t_i / min(hs, hf) - 1e-9))
            a(f"LESIZE,{rad0[i]},,,{ndiv},,1")
            a(f"LESIZE,{rad1[i]},,,{ndiv},,1")
        # the bore arc (contact target) sized like the lamination OD, so that the two sides of the
        # sleeve contact are meshed alike (a finer target than contact side makes the pressure ripple)
        a(f"LESIZE,{arcs[0]},{_num(h)},,,,1")
    elif model.sleeve is not None:
        hs = min(h, _sleeve_area_size(model))  # at least three quadratic elements through the sleeve
        a("ASEL,S,AREA,,301")
        a(f"AESIZE,ALL,{_num(hs)}")
        a("ALLSEL")
        if bore_arc_sized(p, model):
            # thin sleeve (area size finer than mech_mesh), any rotor: the bore arc (line 30003, the
            # contact target) sized like the lamination OD, as for a layered sleeve, so that both sides of
            # the sleeve contact are meshed alike (the finer target of the area size alone made the contact
            # pressure ripple). Written only when the area size is finer than mech_mesh, so every other
            # deck, the default bridge-retained one among them (1.5 mm sleeve: equal sizes), is unchanged
            # (docs/methods.md, 3).
            a(f"LESIZE,30003,{_num(h)},,,,1")
    a("AMESH,ALL")
    # contact -----------------------------------------------------------------------------
    if model.sleeve is not None:
        a("! sleeve / lamination OD contact; interference = contact offset CNOF (ramped)")
        a("ET,10,CONTA172")
        a("ET,11,TARGE169")
        a("KEYOPT,10,2,0")     # augmented Lagrange
        a("KEYOPT,10,9,4")     # offset only, ramped
        a("KEYOPT,10,10,2")    # update contact stiffness each iteration
        # ICONT = 0.005 mm (absolute): coincident sleeve/lamination surfaces start closed, which
        # holds loose lamination pieces (bridgeless rotor) at zero interference
        a("R,10,,,1.0,0.1,-0.005")
        a("RMODIF,10,10,DELTA")
        a(f"MP,MU,10,{_num(p.mu_sleeve)}")
        a("CSYS,11")
        a(f"ESEL,S,MAT,,{slv_sel}")
        a("NSLE,S")
        a(f"NSEL,R,LOC,X,R_LAM-{tol},R_LAM+{tol}")
        a("TYPE,11")
        a("REAL,10")
        a("MAT,10")
        a("ESURF")
        a(f"ESEL,S,MAT,,{MAT_LAM}")
        a("NSLE,S")
        a(f"NSEL,R,LOC,X,R_LAM-{tol},R_LAM+{tol}")
        a("TYPE,10")
        a("ESURF")
        a("ALLSEL")
        a("CSYS,0")
    if model.magnets:
        a("ET,12,CONTA172")
        a("ET,13,TARGE169")
        if p.magnet_contact == "bonded":
            # bonded (always), MPC as the Workbench default for bonded contact. The faces with the
            # designed clearance are tied too: KEYOPT(9) = 1 ignores the initial gap without moving
            # the geometry, and the absolute pinball (PINB < 0) reaches across it.
            pinb = max(3 * model.magnet_clearance, 0.05)
            a("! magnet / pocket contact: bonded (MPC), designed clearance bridged by the bond")
            a("KEYOPT,12,2,2")     # MPC
            a("KEYOPT,12,4,2")     # nodal detection normal to the target (recommended for MPC)
            a("KEYOPT,12,9,1")     # ignore the initial gap (the clearance) and keep the geometry
            a("KEYOPT,12,12,5")    # bonded always
            a(f"R,20,,,1.0,0.1,0,-{_num(pinb)}")
        else:
            a("! magnet / pocket contact: frictional, designed clearance open at the start")
            a("KEYOPT,12,2,0")
            if model.magnet_clearance > 0:
                a("KEYOPT,12,9,0")  # include the initial geometric gap: the magnets have a designed clearance
            else:
                a("KEYOPT,12,9,1")  # ignore initial geometric gap/penetration
            a("KEYOPT,12,10,2")
            a("KEYOPT,12,15,2")    # contact stabilisation damping in all load steps (energy checked below)
            # ICONT (5th real constant) < 0: absolute band in which contact is closed at the start,
            # so the coincident magnet/pocket faces start closed instead of "just open"
            a("R,20,,,1.0,0.1,-0.005")
        a(f"MP,MU,12,{_num(p.mu_magnet)}")
        a(f"ESEL,S,MAT,,{MAT_LAM}")
        a("NSLE,S")
        a("NSEL,R,EXT")
        a("TYPE,13")
        a("REAL,20")
        a("MAT,12")
        a("ESURF")
        a(f"ESEL,S,MAT,,{MAT_PM}")
        a("NSLE,S")
        a("NSEL,R,EXT")
        a("TYPE,12")
        a("ESURF")
        a("ALLSEL")
    # boundary conditions -----------------------------------------------------------------
    a("! symmetry on both radial cut lines (nodal cylindrical axes)")
    a("CSYS,11")
    a(f"NSEL,S,LOC,Y,-{tol},{tol}")
    a(f"NSEL,R,LOC,X,{tol},1e6")
    a("NROTAT,ALL")
    a("D,ALL,UY,0")
    a(f"NSEL,S,LOC,Y,{_num(sector)}-{tol},{_num(sector)}+{tol}")
    a(f"NSEL,R,LOC,X,{tol},1e6")
    a("NROTAT,ALL")
    a("D,ALL,UY,0")
    if model.origin_node:
        a("CSYS,0")
        a(f"NSEL,S,LOC,X,-{tol},{tol}")
        a(f"NSEL,R,LOC,Y,-{tol},{tol}")
        a("D,ALL,UX,0")
        a("D,ALL,UY,0")
    a("ALLSEL")
    a("CSYS,0")
    a("*GET,NNODE,NODE,0,COUNT")
    a("*GET,NELEM,ELEM,0,COUNT")
    a("FINISH")
    # solution ----------------------------------------------------------------------------
    a("/SOLU")
    a("ANTYPE,STATIC")
    a("NLGEOM,OFF")
    a("OMEGA,0,0,OMG")
    if thermal:
        a(f"BFUNIF,TEMP,{_num(p.mech_temperature)}")
    if any(prestrain):
        # winding prestress: the stress-free hoop strain of a layer is -eps0 = -sigma0/E_hoop, given as
        # an initial elastic hoop strain +eps0 (initial stress D eps0) in the element axes (x radial,
        # y hoop), all layers at once; equilibrium is found in the load step with the other loads
        a("! winding prestress: hoop eigenstrain per sleeve layer as initial elastic strain (element axes)")
        a("INISTATE,SET,CSYS,-2")
        a("INISTATE,SET,DTYP,EPEL")
        for mat, e in zip(layer_mats, prestrain):
            if e:
                a(f"ESEL,S,MAT,,{mat}")
                a(f"INISTATE,DEFINE,ALL,ALL,,,0,{_num(e)},0,0,0,0")
        a("ALLSEL")
    a("AUTOTS,ON")
    a("NSUBST,10,200,5")
    a("OUTRES,ALL,LAST")
    a("SOLVE")
    a("*GET,CNVG,ACTIVE,0,SOLU,CNVG")
    a("FINISH")
    # post --------------------------------------------------------------------------------
    a("/POST1")
    a("SET,LAST")
    a("RSYS,0")
    keys = []

    def nsort_get(mat, comp, how, key, rsys=0):
        a(f"RSYS,{rsys}")
        a(f"ESEL,S,MAT,,{mat}")
        a("NSLE,S")
        a(f"NSORT,S,{comp}")
        a(f"*GET,{key},SORT,,{how}")
        a(f"*GET,{key}_N,SORT,,I{how}")
        a(f"*GET,{key}_X,NODE,{key}_N,LOC,X")
        a(f"*GET,{key}_Y,NODE,{key}_N,LOC,Y")
        keys.extend([key, f"{key}_X", f"{key}_Y"])

    nsort_get(MAT_LAM, "EQV", "MAX", "SVM_LAM")
    if model.sections:
        # lamination nodes (x, y, averaged von Mises) for the named stress sections
        a(f"ESEL,S,MAT,,{MAT_LAM}")
        a("NSLE,S")
        a("RSYS,0")
        a("*GET,NDMX,NODE,0,NUM,MAXD")
        for arr in ("NMSK", "NXA", "NYA", "NSA"):
            a(f"*DIM,{arr},ARRAY,NDMX")
        a("*VGET,NMSK(1),NODE,1,NSEL")
        a("*VGET,NXA(1),NODE,1,LOC,X")
        a("*VGET,NYA(1),NODE,1,LOC,Y")
        a("*VGET,NSA(1),NODE,1,S,EQV")
        a(f"*CFOPEN,{Path(NODES_FILE).stem},{Path(NODES_FILE).suffix.lstrip('.')}")
        a("*VMASK,NMSK(1)")
        a("*VWRITE,NXA(1),NYA(1),NSA(1)")
        a("(3E21.12)")
        a("*CFCLOS")
    if model.magnets:
        nsort_get(MAT_PM, "EQV", "MAX", "SVM_PM")
    if layered:
        # per layer: only that layer's elements selected, so the nodal stresses are averaged within the
        # layer (in batch /POST1 a node's result is the average over the selected elements attached to
        # it, whatever their material); the envelopes over the layers are formed in Python
        for i, mat in enumerate(layer_mats):
            tag = f"_L{i + 1}"
            nsort_get(mat, "Y", "MAX", f"SH_SLV{tag}", rsys=11)
            nsort_get(mat, "X", "MIN", f"SR_SLV{tag}", rsys=11)
            nsort_get(mat, "EQV", "MAX", f"SVM_SLV{tag}")
            nsort_get(mat, "X", "MAX", f"SRT_SLV{tag}", rsys=11)
            nsort_get(mat, "XY", "MAX", f"TAU_SLV_HI{tag}", rsys=11)
            nsort_get(mat, "XY", "MIN", f"TAU_SLV_LO{tag}", rsys=11)
        # sleeve corner nodes per layer: r, theta (deg), stresses in system 11 and the layer number
        a("RSYS,11")
        a("CSYS,11")
        a("*GET,SNMX,NODE,0,NUM,MAXD")
        for arr in ("SMSK", "SRA", "STA", "SSX", "SSY", "SSZ", "SSXY", "SLA"):
            a(f"*DIM,{arr},ARRAY,SNMX")
        for i, mat in enumerate(layer_mats):
            a(f"ESEL,S,MAT,,{mat}")
            a("NSLE,S,CORNER")
            a("*VGET,SMSK(1),NODE,1,NSEL")
            a("*VGET,SRA(1),NODE,1,LOC,X")
            a("*VGET,STA(1),NODE,1,LOC,Y")
            for arr, comp in (("SSX", "X"), ("SSY", "Y"), ("SSZ", "Z"), ("SSXY", "XY")):
                a(f"*VGET,{arr}(1),NODE,1,S,{comp}")
            a(f"*VFILL,SLA(1),RAMP,{i + 1},0")
            a(f"*CFOPEN,{Path(SLEEVE_FILE).stem},{Path(SLEEVE_FILE).suffix.lstrip('.')}"
              + (",,APPEND" if i else ""))
            a("*VMASK,SMSK(1)")
            a("*VWRITE,SRA(1),STA(1),SSX(1),SSY(1),SSZ(1),SSXY(1),SLA(1)")
            a("(7E21.12)")
            a("*CFCLOS")
        a("CSYS,0")
    if model.sleeve is not None and not layered:
        nsort_get(MAT_SLV, "Y", "MAX", "SH_SLV", rsys=11)
        nsort_get(MAT_SLV, "X", "MIN", "SR_SLV", rsys=11)
        nsort_get(MAT_SLV, "EQV", "MAX", "SVM_SLV")
        nsort_get(MAT_SLV, "X", "MAX", "SRT_SLV", rsys=11)      # max radial tension (delamination)
        nsort_get(MAT_SLV, "XY", "MAX", "TAU_SLV_HI", rsys=11)  # r-theta shear, both signs
        nsort_get(MAT_SLV, "XY", "MIN", "TAU_SLV_LO", rsys=11)
        # sleeve corner nodes: r, theta (deg) and stresses in the cylindrical system 11
        a(f"ESEL,S,MAT,,{MAT_SLV}")
        a("NSLE,S,CORNER")
        a("RSYS,11")
        a("CSYS,11")
        a("*GET,SNMX,NODE,0,NUM,MAXD")
        for arr in ("SMSK", "SRA", "STA", "SSX", "SSY", "SSZ", "SSXY"):
            a(f"*DIM,{arr},ARRAY,SNMX")
        a("*VGET,SMSK(1),NODE,1,NSEL")
        a("*VGET,SRA(1),NODE,1,LOC,X")
        a("*VGET,STA(1),NODE,1,LOC,Y")
        for arr, comp in (("SSX", "X"), ("SSY", "Y"), ("SSZ", "Z"), ("SSXY", "XY")):
            a(f"*VGET,{arr}(1),NODE,1,S,{comp}")
        a(f"*CFOPEN,{Path(SLEEVE_FILE).stem},{Path(SLEEVE_FILE).suffix.lstrip('.')}")
        a("*VMASK,SMSK(1)")
        a("*VWRITE,SRA(1),STA(1),SSX(1),SSY(1),SSZ(1),SSXY(1)")
        a("(6E21.12)")
        a("*CFCLOS")
        a("CSYS,0")
    if model.sleeve is not None:
        nsort_get(MAT_LAM, "X", "MIN", "SR_LAM", rsys=11)
        a("RSYS,0")
        a("ESEL,S,TYPE,,10")
        a("ETABLE,CPRS,CONT,PRES")
        a("ETABLE,CSTA,CONT,STAT")
        a("ESORT,ETAB,CPRS")
        a("*GET,CP_MAX,SORT,,MAX")
        a("*GET,CP_MIN,SORT,,MIN")
        if contact_elements:
            # per-element contact pressure and status of the sleeve contact (lift-off extrapolation)
            a("*GET,CEMX,ELEM,0,NUM,MAXD")
            for arr in ("CMSK", "CEID", "CEX", "CEY", "CEP", "CES"):
                a(f"*DIM,{arr},ARRAY,CEMX")
            a("*VGET,CMSK(1),ELEM,1,ESEL")
            a("*VFILL,CEID(1),RAMP,1,1")
            a("*VGET,CEX(1),ELEM,1,CENT,X")
            a("*VGET,CEY(1),ELEM,1,CENT,Y")
            a("*VGET,CEP(1),ELEM,1,ETAB,CPRS")
            a("*VGET,CES(1),ELEM,1,ETAB,CSTA")
            a(f"*CFOPEN,{Path(CONTACT_FILE).stem},{Path(CONTACT_FILE).suffix.lstrip('.')}")
            a("*VMASK,CMSK(1)")
            a("*VWRITE,CEID(1),CEX(1),CEY(1),CEP(1),CES(1)")
            a("(5E21.12)")
            a("*CFCLOS")
        a("*GET,NC_ALL,ELEM,0,COUNT")
        a("ESEL,R,ETAB,CSTA,-1,1.5")      # status 0/1 = open
        a("*GET,NC_OPEN,ELEM,0,COUNT")
        a("*SET,OPEN_FR,NC_OPEN/NC_ALL")
        keys += ["CP_MAX", "CP_MIN", "OPEN_FR"]
    a("ALLSEL")
    a("RSYS,0")
    a("ETABLE,SENE_E,SENE")
    a("SSUM")
    a("*GET,SENE_TOT,SSUM,,ITEM,SENE_E")
    if model.magnets or model.sleeve is not None:
        # AENE of contact elements = contact stabilisation energy; must stay << strain energy
        a("ESEL,S,ENAME,,172")
        a("ETABLE,STAB_E,AENE")
        a("SSUM")
        a("*GET,STAB_TOT,SSUM,,ITEM,STAB_E")
        a("ALLSEL")
        keys += ["STAB_TOT"]
    keys += ["SENE_TOT", "NNODE", "NELEM", "CNVG"]
    a(f"*CFOPEN,{Path(results_file).stem},{Path(results_file).suffix.lstrip('.') or 'txt'}")
    for k in keys:
        a(f"*VWRITE,{k}")
        a(f"('{k.lower()} = ',E20.12)")
    a("*VWRITE")
    a("('done = 1')")
    a("*CFCLOS")
    a("FINISH")
    a("/EXIT,NOSAVE")
    path.write_text("\n".join(L) + "\n", encoding="ascii")
    return path


# ======================================================================================
# runner and parser
# ======================================================================================
def parse_results(path: str | Path) -> dict:
    """Read the ``key = value`` file written by the deck."""
    out = {}
    for line in Path(path).read_text(encoding="ascii", errors="replace").splitlines():
        m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\S+)\s*$", line)
        if m:
            try:
                out[m.group(1)] = float(m.group(2))
            except ValueError:
                out[m.group(1)] = m.group(2)
    return out


def parse_sleeve_nodes(path: str | Path) -> dict[str, np.ndarray]:
    """Read the sleeve node file written by the deck: one 1-D array per :data:`SLEEVE_COLUMNS`
    entry (r in mm, theta in deg, stresses in MPa in the cylindrical system). The file of a layered
    deck has a 7th column, returned as ``"layer"`` (int, 1 = bore layer): each layer's nodes with
    the stresses averaged over that layer's elements only, so a node on a layer interface appears
    once for each of the two layers."""
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] not in (len(SLEEVE_COLUMNS), len(SLEEVE_COLUMNS) + 1):
        raise ValueError(f"{path}: expected {len(SLEEVE_COLUMNS)} (or {len(SLEEVE_COLUMNS) + 1} with the layer) "
                         f"columns, got {data.shape[1]}")
    out = {name: data[:, i].copy() for i, name in enumerate(SLEEVE_COLUMNS)}
    if data.shape[1] > len(SLEEVE_COLUMNS):
        out[SLEEVE_LAYER_COLUMN] = np.rint(data[:, -1]).astype(int)
    return out


def parse_contact_elements(path: str | Path) -> dict[str, np.ndarray]:
    """Read the sleeve contact element file (:data:`CONTACT_FILE`): one 1-D array per
    :data:`CONTACT_COLUMNS` entry, sorted by element number (``elem`` as int)."""
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] != len(CONTACT_COLUMNS):
        raise ValueError(f"{path}: expected {len(CONTACT_COLUMNS)} columns, got {data.shape[1]}")
    data = data[np.argsort(data[:, 0], kind="stable")]
    out = {name: data[:, i].copy() for i, name in enumerate(CONTACT_COLUMNS)}
    out["elem"] = np.rint(out["elem"]).astype(int)
    return out


def sleeve_failure_outputs(p: Params, nodes: dict[str, np.ndarray]) -> dict:
    """Sleeve failure columns from the sleeve node arrays (:func:`parse_sleeve_nodes`).

    The nodal stresses are taken as the homogenised sleeve stress (hoop, axial, radial, 0, r-theta
    shear, 0) and evaluated by :func:`sipm.composite.sleeve_failure` with the ply strengths of ``p``:
    with a composite spec in its +/- angle plies (thermal mismatch for ``mech_delta_t``
    included), with directly entered constants as one hoop-wound ply. Returns the maximum over the
    nodes of each failure index (``sleeve_fi_<criterion>``), the minimum Tsai-Wu strength ratio
    ``sleeve_sf`` and the location of that node (``sleeve_sf_r`` mm, ``sleeve_sf_theta`` deg).
    ``sleeve_fi_tsai_wu`` is 1 / ``sleeve_sf`` (the reciprocal Tsai-Wu strength ratio, maximum over
    the nodes), not the raw quadratic of :func:`sipm.composite.failure_indices`.
    """
    from .composite import CRITERIA, sleeve_failure

    n = len(nodes["r"])
    sig = np.zeros((n, 6))
    sig[:, 0], sig[:, 1], sig[:, 2], sig[:, 4] = nodes["s_hoop"], nodes["s_axial"], nodes["s_radial"], nodes["s_shear"]
    if p.is_layered:
        return _layered_failure_outputs(p, nodes, sig)
    fi = sleeve_failure(sig, p.sleeve_strength(), p.sleeve_laminate(), p.mech_delta_t)
    out = {f"sleeve_fi_{c}": float(np.max(fi[c])) for c in CRITERIA}
    i = int(np.argmin(fi["tsai_wu_sr"]))
    sf = float(fi["tsai_wu_sr"][i])
    # Tsai-Wu column: 1 / strength ratio (linear in load, >= 0, 1 on the envelope), not the raw
    # quadratic, whose linear terms make it negative under compression-dominated states
    out["sleeve_fi_tsai_wu"] = 1.0 / sf if sf > 0 else math.inf
    out.update(sleeve_sf=sf, sleeve_sf_r=float(nodes["r"][i]), sleeve_sf_theta=float(nodes["theta"][i]))
    return out


def _layered_failure_outputs(p: Params, nodes: dict[str, np.ndarray], sig: np.ndarray) -> dict:
    """:func:`sleeve_failure_outputs` of a layered sleeve: every node row is evaluated with the
    laminate of the layer it was exported for (a node on an interface once per layer, with that
    layer's own averaged stresses). The envelope columns are those of a single layer; in addition
    ``sleeve_sf_layer`` (layer of the governing node, 1 = bore) and per layer ``sleeve_sf_L<i>``."""
    from .composite import CRITERIA, sleeve_failure

    if SLEEVE_LAYER_COLUMN not in nodes:
        raise ValueError("a layered sleeve needs the per-layer node export (the 'layer' column)")
    layer = np.asarray(nodes[SLEEVE_LAYER_COLUMN], int)
    lams = p.sleeve_layer_laminates()
    if layer.min() < 1 or layer.max() > len(lams):
        raise ValueError(f"sleeve node export has layers {layer.min()}..{layer.max()}, the sleeve has {len(lams)}")
    fi = {k: np.full(len(layer), np.nan) for c in CRITERIA for k in (c, f"{c}_sr")}
    for k, lam in enumerate(lams, start=1):
        m = layer == k
        if m.any():
            part = sleeve_failure(sig[m], p.sleeve_strength(), lam, p.mech_delta_t)
            for key in fi:
                fi[key][m] = part[key]
    out = {f"sleeve_fi_{c}": float(np.nanmax(fi[c])) for c in CRITERIA}
    srs = fi["tsai_wu_sr"]
    i = int(np.nanargmin(srs))
    sf = float(srs[i])
    out["sleeve_fi_tsai_wu"] = 1.0 / sf if sf > 0 else math.inf
    out.update(sleeve_sf=sf, sleeve_sf_r=float(nodes["r"][i]), sleeve_sf_theta=float(nodes["theta"][i]),
               sleeve_sf_layer=int(layer[i]))
    for k in range(1, len(lams) + 1):
        m = layer == k
        out[f"sleeve_sf_L{k}"] = float(np.min(srs[m])) if m.any() else math.nan
    return out


def sleeve_interface_outputs(nodes: dict[str, np.ndarray]) -> dict:
    """Interlaminar r-theta shear of a layered sleeve from the per-layer node export.

    At each layer interface the nodes on the interface radius (within :data:`INTERFACE_TOL`) are
    exported twice, once with each layer's averaged stresses; the shear traction is continuous, so
    the two are averaged node by node (paired by theta) as the estimate of the interface value.
    Returns ``tau_sleeve_interface_max``, the largest |r-theta shear| over all interfaces, or {} for
    fewer than two layers."""
    if SLEEVE_LAYER_COLUMN not in nodes:
        return {}
    layer = np.asarray(nodes[SLEEVE_LAYER_COLUMN], int)
    nl = int(layer.max()) if len(layer) else 0
    if nl < 2:
        return {}
    r, th, tau = nodes["r"], nodes["theta"], nodes["s_shear"]
    best = 0.0
    for k in range(1, nl):
        lo, hi = layer == k, layer == k + 1
        if not (lo.any() and hi.any()):
            continue
        r_if = 0.5 * (r[lo].max() + r[hi].min())
        a = np.flatnonzero(lo & (np.abs(r - r_if) < INTERFACE_TOL))
        b = np.flatnonzero(hi & (np.abs(r - r_if) < INTERFACE_TOL))
        a, b = a[np.argsort(th[a], kind="stable")], b[np.argsort(th[b], kind="stable")]
        if len(a) == len(b) and np.allclose(th[a], th[b], atol=1e-6):
            vals = 0.5 * (tau[a] + tau[b])
        else:   # not paired (should not happen with shared nodes): the larger of the two sides
            vals = np.concatenate([tau[a], tau[b]])
        if len(vals):
            best = max(best, float(np.max(np.abs(vals))))
    return {"tau_sleeve_interface_max": best}


def errors_from_output(path: str | Path, limit: int = 5) -> list[str]:
    """Text of the first *** ERROR *** blocks of a MAPDL output file."""
    p = Path(path)
    if not p.exists():
        return []
    lines = p.read_text(encoding="latin-1", errors="replace").splitlines()
    errs = []
    for i, ln in enumerate(lines):
        if "*** ERROR ***" in ln:
            block = [ln.strip()]
            for nxt in lines[i + 1:i + 6]:
                if not nxt.strip():
                    break
                block.append(nxt.strip())
            errs.append(" ".join(block))
            if len(errs) >= limit:
                break
    return errs


class MapdlError(RuntimeError):
    pass


def run_deck(deck: Path, workdir: Path, cfg: Config, job: str = "sipm", on_start=None) -> Path:
    """Run MAPDL in batch on ``deck`` inside ``workdir``; return the results-file path."""
    exe = cfg.resolved_mapdl()
    workdir.mkdir(parents=True, exist_ok=True)
    for f in (workdir / RESULTS_FILE, workdir / NODES_FILE, workdir / SLEEVE_FILE, workdir / CONTACT_FILE,
              workdir / f"{job}.lock"):
        if f.exists():
            f.unlink()
    out = workdir / f"{job}.out"
    cmd = [str(exe), "-b", "-smp", "-np", str(cfg.mapdl_cores), "-dir", str(workdir), "-j", job,
           "-i", str(deck), "-o", str(out)]
    proc = subprocess.Popen(cmd, cwd=workdir, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL, env=child_environment())
    if on_start is not None:
        on_start(proc.pid)  # lets the sweep's memory watchdog follow this process tree
    try:
        rc = proc.wait(timeout=cfg.mapdl_timeout)
    except BaseException:
        _kill_tree(proc.pid)
        raise
    res = workdir / RESULTS_FILE
    if not res.exists():
        errs = errors_from_output(out)
        raise MapdlError(f"MAPDL exit code {rc}; no results file. "
                         + (" | ".join(errs) if errs else f"see {out}"))
    return res


def _kill_tree(pid: int) -> None:
    try:
        import psutil

        parent = psutil.Process(pid)
        for c in parent.children(recursive=True):
            c.kill()
        parent.kill()
    except Exception:
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(pid)], capture_output=True)


def _cleanup(workdir: Path, job: str) -> None:
    for ext in (".rst", ".db", ".esav", ".full", ".emat", ".mntr", ".ldhi", ".r001", ".stat", ".page",
                ".lock", ".BCS", ".DSP", ".err", ".log", "0.err", "0.log", ".osav", ".rdb", ".cnd"):
        f = workdir / f"{job}{ext}"
        if f.exists():
            try:
                f.unlink()
            except OSError:
                pass
    for f in workdir.glob(f"{job}*.r0*"):
        try:
            f.unlink()
        except OSError:
            pass


def solve_mech(p: Params, cfg: Config | None = None, workdir: str | Path | None = None,
               model: MechModel | None = None, on_start=None) -> dict:
    """Build, run and post-process one mechanical point. Returns a flat result dict.

    Stresses in MPa (nodal von Mises averaged within each body):

    * ``svm_lam_max`` - the headline lamination stress: peak over the named load-carrying
      sections of :mod:`sipm.geometry.locations` (each section's own peak in ``svm_<name>``,
      e.g. ``svm_bridge1``; the section that holds the peak in ``svm_lam_where``). A model
      without sections (validation disc) reports the global lamination peak here.
    * ``svm_lam_max_all`` - global lamination peak over all nodes. It sits at sharp corners of
      the outline or at magnet contact edges, which are singular: mesh-dependent.
    * ``svm_pm_max`` - magnet peak (at contact edges, mesh-dependent).
    * ``svm_rotor_max`` - the design quantity used by the plots, = ``svm_lam_max``.
    * ``sh_sleeve_max`` - max hoop stress in the sleeve; ``sr_sleeve_min`` its most compressive
      radial stress.
    * ``sr_sleeve_max`` - max radial stress in the sleeve (radial tension, the delamination
      driver; <= 0 when the sleeve is compressed through its thickness everywhere).
    * ``tau_sleeve_max`` - max absolute radial-hoop shear stress in the sleeve.
    * ``sleeve_fi_max_stress``, ``sleeve_fi_tsai_wu``, ``sleeve_fi_hashin_fibre``,
      ``sleeve_fi_hashin_matrix`` - max ply failure index over the sleeve nodes (>= 1 = failure;
      the Tsai-Wu one is 1 / ``sleeve_sf``, linear in load);
      ``sleeve_sf`` - min Tsai-Wu strength ratio (safety factor) at ``sleeve_sf_r``,
      ``sleeve_sf_theta`` (see :func:`sleeve_failure_outputs`).

    * layered sleeve only: ``sh_sleeve_max_L<i>`` and ``sleeve_sf_L<i>`` per layer (1 = bore),
      ``sleeve_sf_layer`` (layer of ``sleeve_sf``) and ``tau_sleeve_interface_max`` (largest
      interlaminar r-theta shear, :func:`sleeve_interface_outputs`); the other sleeve columns are
      extremes over the layers.

    Sleeve stresses are nodal averages in the cylindrical system (x radial, y hoop), within each
    layer for a layered sleeve. With a
    sleeve, ``r["_sleeve_nodes"]`` holds the per-node arrays of :func:`parse_sleeve_nodes`
    (corner nodes; the file is deleted with the other solver files unless
    ``cfg.keep_mapdl_files``). Like ``_raw`` it is not written to the sweep cache.

    With ``p.mech_cases`` set the point is a load-case set: see :func:`solve_mech_cases`.
    """
    cfg = cfg or Config()
    t0 = time.time()
    if model is None:
        rp = build_rotor(p)
        model = mech_model_from_rotor(rp)
    wd = Path(workdir) if workdir else cfg.workdir / "mech_single"
    if p.mech_cases:
        return solve_mech_cases(p, cfg, wd, model, on_start)
    return _solve_single(p, cfg, wd, model, on_start, t0=t0)


def solve_mech_cases(p: Params, cfg: Config | None = None, workdir: str | Path | None = None,
                     model: MechModel | None = None, on_start=None, liftoff_search: bool | None = None) -> dict:
    """Solve the load cases of ``p.mech_cases`` and combine them (:func:`sipm.loadcases.combine`).

    Each case is its own MAPDL run of the single-case deck (:meth:`Params.for_case`: speed and
    temperature of the case), in ``<workdir>/<case>``, solved in the canonical order of
    :data:`sipm.params.MECH_CASES`; cases with the same speed and temperature are solved once.
    Every case is loaded from the fitted state in one load step, as a single point is, so a
    case's result does not depend on which other cases are requested or on their order (with or
    without friction), and a case that fails to converge cannot spoil the start of another. The
    decks export the per-element sleeve contact pressures for the lift-off speed. A failing case
    fails the whole point (MapdlError naming the case).

    Lift-off: the two-point extrapolation of :mod:`sipm.loadcases`, except for a rotor whose
    contact is not linear in speed (:func:`sipm.loadcases.search_applies`, i.e. sleeved_bridgeless
    with ``liftoff_search_tol`` > 0), where :func:`sipm.loadcases.liftoff_search` solves extra
    single-case runs at trial speeds (``<workdir>/liftoff_<label>_<n>``, the deck of the cases at
    another speed). ``liftoff_search`` True / False forces the search on / off (validation and
    tests; None = by that rule); forced on, a ``liftoff_search_tol`` of 0 is an error. The trial
    runs are listed in ``r["_liftoff_trials"]`` (label -> list of dicts; not cached)."""
    from . import loadcases

    cfg = cfg or Config()
    t0 = time.time()
    if model is None:
        model = mech_model_from_rotor(build_rotor(p))
    wd = Path(workdir) if workdir else cfg.workdir / "mech_single"
    solved: dict[tuple[float, float], dict] = {}
    results = {}
    contact = model.sleeve is not None
    for case in loadcases.canonical(p.mech_case_list):
        q = p.for_case(case)
        key = (q.mech_rpm, q.mech_temperature)
        if key not in solved:
            try:
                solved[key] = _solve_single(q, cfg, wd / case, model, on_start, contact_elements=contact)
            except MapdlError as e:
                raise MapdlError(f"load case {case} ({q.mech_rpm:g} rpm, {q.mech_temperature:g} degC): {e}") from e
        results[case] = solved[key]
    search = loadcases.search_applies(p) if liftoff_search is None else liftoff_search
    if not (search and contact):
        return loadcases.combine(p, results, wall_s=time.time() - t0)
    searched, trials, by_temp = {}, {}, {}
    for label, temp in (("cold", p.mech_temperature_cold), ("hot", p.mech_temperature)):
        if temp not in by_temp:   # equal temperatures: one search for both labels
            n = [0]

            def solve_at(rpm: float, temp=temp, label=label, n=n) -> dict:
                n[0] += 1
                q = p.replace(mech_cases="", mech_rpm=rpm, mech_temperature=temp)
                try:
                    return _solve_single(q, cfg, wd / f"liftoff_{label}_{n[0]}", model, on_start,
                                         contact_elements=True)["_contact_elems"]
                except MapdlError as e:
                    raise MapdlError(f"lift-off search run at {rpm:.1f} rpm, {temp:g} degC: {e}") from e

            entries = loadcases.liftoff_entries(p, results, temp)
            rpm, basis, tr = loadcases.liftoff_search(entries, solve_at, p.liftoff_search_tol)
            by_temp[temp] = ((rpm, basis), tr)
        searched[label], trials[label] = by_temp[temp]
    out = loadcases.combine(p, results, wall_s=time.time() - t0, searched=searched)
    out["_liftoff_trials"] = trials
    return out


def _solve_single(p: Params, cfg: Config, wd: Path, model: MechModel, on_start=None,
                  contact_elements: bool = False, t0: float | None = None) -> dict:
    """One MAPDL run of the single-case deck of ``p`` (``p.mech_cases`` is not looked at)."""
    t0 = time.time() if t0 is None else t0
    wd.mkdir(parents=True, exist_ok=True)
    deck = write_deck(p, model, wd / "sipm.inp", contact_elements=contact_elements)
    res_path = run_deck(deck, wd, cfg, on_start=on_start)
    raw = parse_results(res_path)
    if raw.get("done") != 1:
        raise MapdlError(f"incomplete results file {res_path}")
    if raw.get("cnvg", 0) != 1:
        errs = errors_from_output(wd / "sipm.out")
        raise MapdlError("MAPDL solution did not converge" + (": " + " | ".join(errs) if errs else ""))
    r = {}
    if model.sections:
        from .geometry.locations import section_peaks

        nodes = np.loadtxt(wd / NODES_FILE, ndmin=2)
        peaks = section_peaks(model.sections, nodes[:, :2], nodes[:, 2])
        for name, pk in peaks.items():
            r[f"svm_{name}"] = pk["value"]
        where = max(peaks, key=lambda n: peaks[n]["value"])
        r["svm_lam_max"] = peaks[where]["value"]
        r["svm_lam_where"] = where
        r["svm_lam_r"] = math.hypot(peaks[where]["x"], peaks[where]["y"])
        r["svm_lam_theta"] = math.degrees(math.atan2(peaks[where]["y"], peaks[where]["x"]))
    else:
        r["svm_lam_max"] = raw["svm_lam"]
        r["svm_lam_r"] = math.hypot(raw["svm_lam_x"], raw["svm_lam_y"])
        r["svm_lam_theta"] = math.degrees(math.atan2(raw["svm_lam_y"], raw["svm_lam_x"]))
    r["svm_lam_max_all"] = raw["svm_lam"]
    r["svm_lam_all_r"] = math.hypot(raw["svm_lam_x"], raw["svm_lam_y"])
    r["svm_lam_all_theta"] = math.degrees(math.atan2(raw["svm_lam_y"], raw["svm_lam_x"]))
    if "svm_pm" in raw:
        r["svm_pm_max"] = raw["svm_pm"]
    r["svm_rotor_max"] = r["svm_lam_max"]
    nl = 0
    while f"sh_slv_l{nl + 1}" in raw:
        nl += 1
    if nl:
        # layered sleeve: envelopes over the per-layer values (each averaged within its layer)
        ks = range(1, nl + 1)
        r["sh_sleeve_max"] = max(raw[f"sh_slv_l{k}"] for k in ks)
        r["sr_sleeve_min"] = min(raw[f"sr_slv_l{k}"] for k in ks)
        r["sr_sleeve_max"] = max(raw[f"srt_slv_l{k}"] for k in ks)
        r["tau_sleeve_max"] = max(max(abs(raw[f"tau_slv_hi_l{k}"]), abs(raw[f"tau_slv_lo_l{k}"])) for k in ks)
        r["svm_sleeve_max"] = max(raw[f"svm_slv_l{k}"] for k in ks)
    if "sh_slv" in raw:
        r["sh_sleeve_max"] = raw["sh_slv"]
        r["sr_sleeve_min"] = raw["sr_slv"]
        r["sr_sleeve_max"] = raw["srt_slv"]
        r["tau_sleeve_max"] = max(abs(raw["tau_slv_hi"]), abs(raw["tau_slv_lo"]))
        r["svm_sleeve_max"] = raw["svm_slv"]
    if "sr_lam" in raw:
        r["sr_lam_min"] = raw["sr_lam"]
        r["contact_p_max"] = raw["cp_max"]
        r["contact_p_min"] = raw["cp_min"]
        # CONTA172 status: 0/1 open, 2 sliding, 3 sticking
        r["sleeve_open_fraction"] = raw["open_fr"]
    for k in range(1, nl + 1):
        r[f"sh_sleeve_max_L{k}"] = raw[f"sh_slv_l{k}"]
    if "stab_tot" in raw:
        r["stab_energy_ratio"] = raw["stab_tot"] / raw["sene_tot"] if raw["sene_tot"] else float("nan")
    r["n_nodes"] = int(raw["nnode"])
    r["n_elements"] = int(raw["nelem"])
    r["mech_wall_s"] = time.time() - t0
    r["_raw"] = raw
    if model.sleeve is not None:
        r["_sleeve_nodes"] = parse_sleeve_nodes(wd / SLEEVE_FILE)
        r.update(sleeve_failure_outputs(p, r["_sleeve_nodes"]))
        r.update(sleeve_interface_outputs(r["_sleeve_nodes"]))
        if contact_elements:
            r["_contact_elems"] = parse_contact_elements(wd / CONTACT_FILE)
    if not cfg.keep_mapdl_files:
        _cleanup(wd, "sipm")
        (wd / NODES_FILE).unlink(missing_ok=True)
        (wd / SLEEVE_FILE).unlink(missing_ok=True)
        (wd / CONTACT_FILE).unlink(missing_ok=True)
    return r


def copy_deck_only(p: Params, path: str | Path) -> Path:
    """Write the deck for a design without running it (inspection / manual runs). With
    ``mech_cases`` set this is the deck of ``mech_rpm`` at ``mech_temperature``; the deck of one
    case is ``copy_deck_only(p.for_case(case), path)``."""
    rp = build_rotor(p)
    model = mech_model_from_rotor(rp)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return write_deck(p, model, path)


__all__ = ["MechModel", "mech_model_from_rotor", "mech_model_ring", "write_deck", "parse_results",
           "parse_sleeve_nodes", "parse_contact_elements", "sleeve_failure_outputs", "sleeve_interface_outputs",
           "run_deck", "solve_mech",
           "solve_mech_cases", "errors_from_output", "MapdlError", "copy_deck_only",
           "SLEEVE_FILE", "SLEEVE_COLUMNS", "SLEEVE_LAYER_COLUMN", "CONTACT_FILE", "CONTACT_COLUMNS", "MAT_LAYER0"]
