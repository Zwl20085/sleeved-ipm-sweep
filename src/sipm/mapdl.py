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
  (element axes x = radial, y = hoop through ESYS of a local cylindrical system);
* sleeve/lamination: surface-to-surface contact, interference applied as a contact
  offset (CNOF) ramped over the load step, friction ``mu_sleeve``;
* magnet/pocket (``magnet_contact``): 'bonded' (default, as the original Workbench model):
  bonded MPC contact that also ties the faces across the designed clearance (the gap is
  kept in the geometry and bridged by the bond); 'frictional': surface-to-surface contact
  with friction ``mu_magnet``, faces that touch start closed, faces with the designed
  clearance start open, contact stabilisation damping keeps the loose magnets from drifting
  (not supported for the bridgeless rotor);
* lamination bore traction-free (no shaft modelled); units mm, t, s, N, MPa;
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
from .params import Params

RESULTS_FILE = "sipm_results.txt"
NODES_FILE = "sipm_lam_nodes.txt"   # lamination nodes: x, y, von Mises (named sections)
UNSUPPORTED_FRICTIONAL = (
    "magnet_contact='frictional' is not supported for sleeved_bridgeless: the loose magnets in "
    "open-clearance contact do not converge reliably in that rotor (docs/validation.md); use "
    "magnet_contact='bonded'")

MAT_LAM, MAT_PM, MAT_SLV = 1, 2, 3


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


def write_deck(p: Params, model: MechModel, path: str | Path, results_file: str = RESULTS_FILE) -> Path:
    """Write the APDL input deck for one mechanical point and return its path."""
    if p.magnet_contact == "frictional" and p.is_bridgeless and model.magnets:
        raise ValueError(UNSUPPORTED_FRICTIONAL)
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
    if model.sleeve is not None:
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
    if model.sleeve is not None:
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
    if model.sleeve is not None:
        rr = np.hypot(model.sleeve.pts[:, 0], model.sleeve.pts[:, 1])
        hs = min(h, (rr.max() - rr.min()) / 3)  # at least three quadratic elements through the sleeve
        a("ASEL,S,AREA,,301")
        a(f"AESIZE,ALL,{_num(hs)}")
        a("ALLSEL")
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
        a(f"ESEL,S,MAT,,{MAT_SLV}")
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
    if model.sleeve is not None:
        nsort_get(MAT_SLV, "Y", "MAX", "SH_SLV", rsys=11)
        nsort_get(MAT_SLV, "X", "MIN", "SR_SLV", rsys=11)
        nsort_get(MAT_SLV, "EQV", "MAX", "SVM_SLV")
        nsort_get(MAT_LAM, "X", "MIN", "SR_LAM", rsys=11)
        a("RSYS,0")
        a("ESEL,S,TYPE,,10")
        a("ETABLE,CPRS,CONT,PRES")
        a("ETABLE,CSTA,CONT,STAT")
        a("ESORT,ETAB,CPRS")
        a("*GET,CP_MAX,SORT,,MAX")
        a("*GET,CP_MIN,SORT,,MIN")
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
    for f in (workdir / RESULTS_FILE, workdir / NODES_FILE, workdir / f"{job}.lock"):
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
    * ``sh_sleeve_max`` - max hoop stress in the sleeve.
    """
    cfg = cfg or Config()
    t0 = time.time()
    if model is None:
        rp = build_rotor(p)
        model = mech_model_from_rotor(rp)
    wd = Path(workdir) if workdir else cfg.workdir / "mech_single"
    wd.mkdir(parents=True, exist_ok=True)
    deck = write_deck(p, model, wd / "sipm.inp")
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
    if "sh_slv" in raw:
        r["sh_sleeve_max"] = raw["sh_slv"]
        r["sr_sleeve_min"] = raw["sr_slv"]
        r["svm_sleeve_max"] = raw["svm_slv"]
        r["sr_lam_min"] = raw["sr_lam"]
        r["contact_p_max"] = raw["cp_max"]
        r["contact_p_min"] = raw["cp_min"]
        # CONTA172 status: 0/1 open, 2 sliding, 3 sticking
        r["sleeve_open_fraction"] = raw["open_fr"]
    if "stab_tot" in raw:
        r["stab_energy_ratio"] = raw["stab_tot"] / raw["sene_tot"] if raw["sene_tot"] else float("nan")
    r["n_nodes"] = int(raw["nnode"])
    r["n_elements"] = int(raw["nelem"])
    r["mech_wall_s"] = time.time() - t0
    r["_raw"] = raw
    if not cfg.keep_mapdl_files:
        _cleanup(wd, "sipm")
        (wd / NODES_FILE).unlink(missing_ok=True)
    return r


def copy_deck_only(p: Params, path: str | Path) -> Path:
    """Write the deck for a design without running it (inspection / manual runs)."""
    rp = build_rotor(p)
    model = mech_model_from_rotor(rp)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return write_deck(p, model, path)


__all__ = ["MechModel", "mech_model_from_rotor", "mech_model_ring", "write_deck", "parse_results",
           "run_deck", "solve_mech", "errors_from_output", "MapdlError", "copy_deck_only"]
