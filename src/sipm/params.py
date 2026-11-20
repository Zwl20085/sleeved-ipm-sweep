"""Machine, material and analysis parameters.

All lengths are in mm, angles in degrees, moduli in GPa, densities in kg/m^3,
speeds in rpm and currents in A unless the field name says otherwise.

Every field carries a ``source`` tag:

* ``maxwell``      - the value of the author's original Maxwell 2D models at their saved
                     state (sleeved bridged ``01_ICEM_01_simplerotor``, bridgeless
                     ``02_ICEM_01_teslamotor1``); the
                     geometry is constructed the way those models construct it
                     (docs/methods.md, section 1). ``maxwell_name`` gives the model variable;
* ``paper``        - value published in Table I / Table II / Section III of the
                     companion paper (Zhang, Wu, Hua, ICEMS 2026);
* ``illustrative`` - a sensible default chosen for this toolchain. It is NOT
                     the machine of the paper; the paper does not publish it.

Rotor geometry follows the Maxwell construction: each magnet layer is a straight slot fixed
by its d-axis crossing, angle and width; the magnet *length* is derived (see
:mod:`sipm.geometry.construction`), so it is not a parameter.

``docs/parameters.md`` is generated from this module (``sipm params --markdown``).
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass, field, fields
from typing import Any

TOPOLOGIES = ("sleeved_bridged", "sleeved_bridgeless")
# Topologies of earlier versions that the package no longer builds -> reason given to the user.
REMOVED_TOPOLOGIES = {
    "conventional": "the conventional (unsleeved) rotor was removed from the package; its geometry is not "
                    "distributed and cannot be rebuilt. Use 'sleeved_bridged' or 'sleeved_bridgeless'",
}
SOURCES = ("maxwell", "paper", "illustrative")


def _f(default, unit, source, doc, maxwell_name="", **kw):
    return field(default=default, metadata=dict(unit=unit, source=source, doc=doc, maxwell_name=maxwell_name, **kw))


@dataclass
class Params:
    """Flat parameter set for one design point.

    Use :func:`default_params` to obtain topology-specific defaults, and
    :meth:`replace` to derive variants (the result is validated).
    """

    topology: str = _f("sleeved_bridged", "-", "paper",
                       "Rotor topology: sleeved_bridged or sleeved_bridgeless")

    # ---- stator (Maxwell model; paper Table I where marked) ----------------------------------
    stator_od: float = _f(210.0, "mm", "maxwell",
                          "Stator outer diameter (model Rso = 105 mm; the paper's Table I lists 220 mm)",
                          "2*Rso")
    stator_id: float = _f(147.84, "mm", "maxwell", "Stator inner (bore) diameter (paper Table I agrees)", "2*Rsi")
    slots: int = _f(48, "-", "paper", "Number of stator slots")
    poles: int = _f(8, "-", "paper", "Number of rotor poles")
    slot_opening: float = _f(1.9, "mm", "maxwell", "Stator slot opening width (paper Table I agrees)",
                             "Wslot_open")
    slot_tip_depth: float = _f(1.22, "mm", "maxwell",
                               "Radial length of the parallel-sided slot opening, measured from the point where "
                               "the opening side meets the bore", "Wtip")
    slot_tip_angle: float = _f(70.0, "deg", "maxwell",
                               "Angle between the slot centre line and the tapered tooth-tip side", "tip_angle")
    slot_width: float = _f(5.27, "mm", "maxwell", "Width of the parallel-sided slot body", "Wslot")
    slot_depth: float = _f(15.75, "mm", "maxwell",
                           "Slot bottom radius minus bore radius (flat bottom at Rsi + slot_depth)", "Hslot")
    slot_fillet: float = _f(0.3, "mm", "maxwell", "Fillet radius at every slot corner (incl. the bore corners)")
    airgap: float = _f(1.32, "mm", "maxwell",
                       "Physical single-side air gap (sleeve OD to stator bore); held fixed (paper agrees)",
                       "g - Wcf")
    cond_width: float = _f(4.3, "mm", "maxwell", "Conductor width (tangential)", "Wcoil")
    cond_height: float = _f(1.57, "mm", "maxwell", "Conductor height (radial)", "Lcoil")
    cond_start_radius: float = _f(76.8, "mm", "maxwell", "Radius of the slot-opening side of the innermost conductor",
                                  "Rcoil_start")
    cond_count: int = _f(8, "-", "maxwell",
                         "Conductor rectangles drawn per slot (touching, merged into one coil region per slot)",
                         "DuplicateAlongLine NumClones")
    hairpin_layers: int = _f(7, "-", "maxwell", "Conductors (turns) per slot assigned to the coil (paper Table I "
                             "agrees)", "Nac")
    stack_length: float = _f(130.0, "mm", "maxwell", "Active stack length (2D model depth)", "lstk")
    parallel_paths: int = _f(2, "-", "maxwell", "Parallel paths per phase")

    # ---- rotor magnet layers (Maxwell RD_ construction) --------------------------------------
    pm1_d_crossing: float = _f(56.7305951647269, "mm", "maxwell",
                               "Layer 1 (outer): distance from the shaft axis to the point where the magnet "
                               "centre line crosses the d-axis; fixes the radial position of the V",
                               "RD_xpm_lay1")
    pm2_d_crossing: float = _f(45.0864009239181, "mm", "maxwell", "As pm1_d_crossing for layer 2 (inner)",
                               "RD_xpm_lay2")
    pm1_angle: float = _f(97.007958474619, "deg", "maxwell",
                          "Layer-1 included V angle = 2 x angle between magnet axis and d-axis (paper: 97)",
                          "2*RD_thetapm_lay1")
    pm2_angle: float = _f(95.0289047910392, "deg", "maxwell", "Layer-2 included V angle (paper: 95.02)",
                          "2*RD_thetapm_lay2")
    pm1_width: float = _f(3.64900294777709, "mm", "maxwell",
                          "Layer-1 slot width = nominal magnet thickness along magnetisation (paper: 3.65)",
                          "RD_Wpm_lay1")
    pm2_width: float = _f(4.30590729759734, "mm", "maxwell", "Layer-2 slot width (paper: 4.31)", "RD_Wpm_lay2")
    pm1_end_gap: float = _f(1.0, "mm", "maxwell",
                            "Layer 1: distance along the pole-side slot wall from the slot's inner corner (on the "
                            "rib line) to the magnet's inner end", "RD_dpm_lay1")
    pm2_end_gap: float = _f(2.0, "mm", "maxwell", "As pm1_end_gap for layer 2", "RD_dpm_lay2")
    pm1_length_ratio: float = _f(0.899999999211849, "-", "maxwell",
                                 "Layer 1: magnet length / free slot length (pole-side wall from the magnet's inner "
                                 "end to the bridge circle); the magnet length is derived from it", "RD_K_Lpm_lay1")
    pm2_length_ratio: float = _f(0.9, "-", "maxwell", "As pm1_length_ratio for layer 2", "RD_K_Lpm_lay2")
    pm1_rib_half_width: float = _f(1.0, "mm", "maxwell",
                                   "Layer 1: HALF the width of the centre rib on the d-axis; the slot's inner end "
                                   "is the straight line y = this value (+ pm1_inner_offset)", "RD_FB_lay1")
    pm2_web_width: float = _f(2.31010869033253, "mm", "maxwell",
                              "Layer 2: thickness of each of the two steel webs between the tilted inner end of "
                              "the slot and the central cavity, normal to the tilted end; the slot's inner "
                              "pole-side corner sits this value + pm2_inner_offset above the d-axis",
                              "RD_FB_lay2")
    pm1_inner_offset: float = _f(0.0, "mm", "maxwell",
                                 "Layer 1: extra height of the slot's inner end above the rib", "RD_INNO_lay1")
    pm2_inner_offset: float = _f(1.05729199744521, "mm", "maxwell",
                                 "Layer 2: extra height of the slot's inner end above the web", "RD_INNO_lay2")
    pm2_inner_tilt: float = _f(9.14233419891509, "deg", "maxwell",
                               "Layer 2: tilt of the slot's inner end against the d-axis", "RD_thetaINNO_lay2")
    pm1_bridge_thickness: float = _f(1.08008400158923, "mm", "maxwell",
                                     "Layer 1: radial thickness of the outer bridge, measured inwards from the "
                                     "lamination OD; the slot is closed by the circle r_lam - this value (in the "
                                     "bridgeless rotor it only sets the free slot length, i.e. the magnet length)",
                                     "RD_WO_lay1")
    pm2_bridge_thickness: float = _f(1.0, "mm", "maxwell", "As pm1_bridge_thickness for layer 2", "RD_WO_lay2")
    pm1_slot_start: float = _f(1.0, "mm", "maxwell",
                               "Sleeved rotors: distance along the magnet axis from the d-axis crossing to the "
                               "square inner end of the layer-1 slot (strip)", "Rectangle1 XStart")
    pm2_slot_start: float = _f(5.0, "mm", "maxwell",
                               "As pm1_slot_start for layer 2 (model: strip from 0 mm, inner edge moved 5 mm)",
                               "Rectangle11 + MoveEdges")

    # ---- magnet fit (sleeved rotors) ------------------------------------------------------------
    pm_side_clearance: float = _f(0.1, "mm", "maxwell",
                                  "Magnet thinner than the slot: its q-side face is moved this far from the "
                                  "q-side slot wall (the pole-side face stays on the wall)", "MoveEdges -0.1mm")
    pm_end_clearance: float = _f(0.1, "mm", "maxwell",
                                 "Magnet's inner end face moved this far towards the outer end (bridgeless "
                                 "model: 0)", "MoveEdges -0.1mm")
    pm_corner_radius: float = _f(1.0, "mm", "maxwell", "Radius of the four magnet corners (0 = sharp)", "Fillet 1mm")

    # ---- sleeve, holes, shaft --------------------------------------------------------------
    sleeve_thickness: float = _f(1.5, "mm", "paper",
                                 "Carbon-fibre sleeve radial thickness (> 0); the lamination OD is "
                                 "stator_id/2 - airgap - this value. Paper compares 1.5 mm (bridged) and 2.0 mm "
                                 "(bridgeless)", "Wcf")
    shaft_diameter: float = _f(60.0, "mm", "maxwell", "Lamination bore diameter")
    q_hole_inner_radius: float = _f(42.0, "mm", "maxwell",
                                    "Lightening hole on each q-axis: radius of its inner arc (centred on the "
                                    "shaft axis)")
    q_hole_apex_radius: float = _f(56.0, "mm", "maxwell", "Lightening hole: radius of its outer apex on the q-axis")
    q_hole_half_angle: float = _f(20.0, "deg", "maxwell",
                                  "Lightening hole: half the angular span of its inner arc, about the q-axis")
    q_hole_fillet: float = _f(3.0, "mm", "maxwell", "Lightening hole: radius of its three corner fillets")
    d_hole_diameter: float = _f(4.0, "mm", "maxwell", "Round hole on the d-axis near the OD")
    d_hole_depth: float = _f(4.1, "mm", "maxwell",
                             "Lamination OD radius minus the d-axis hole centre radius (model: centre at "
                             "68.5 mm - Wcf)")

    # ---- bridgeless rotor ----------------------------------------------------------------
    lip_depth: float = _f(1.0, "mm", "maxwell",
                          "Bridgeless: radial depth of the pole-tip lips that cover the slot mouths under the OD",
                          "w_tip")
    pm1_opening_angle: float = _f(1.0, "deg", "maxwell",
                                  "Bridgeless: angular width of the OD opening of the layer-1 slot", "rib2_open")
    pm2_opening_angle: float = _f(1.0, "deg", "maxwell",
                                  "Bridgeless: angular width of the OD opening of the layer-2 slot", "rib1_open")
    pm1_opening_shift: float = _f(0.0, "deg", "illustrative",
                                  "Bridgeless: the q-side edge of the layer-1 opening lies this angle towards the "
                                  "d-axis from the opening reference angle (where the magnet centre line meets the "
                                  "OD). Model: 1.5 deg, compensating its layer-1 reference typo; 0 with the true "
                                  "reference (docs/methods.md, 1.4)", "1.5deg in Angle_rib3")
    pm2_opening_shift: float = _f(0.0, "deg", "maxwell", "As pm1_opening_shift for layer 2")

    # ---- feasibility thresholds (illustrative) ----------------------------------------------
    min_web: float = _f(0.5, "mm", "illustrative",
                        "Minimum steel between pockets/holes, and between pockets and the bore (feasibility check)")
    min_q_web: float = _f(1.0, "mm", "illustrative",
                          "Minimum steel width on the q-axis between adjacent poles (feasibility check)")

    # ---- electromagnetic materials -------------------------------------------------------
    lam_bh: str = _f("generic", "-", "illustrative",
                     "Lamination BH curve: 'generic' (labelled generic NO-steel curve), a path to a "
                     "two-column CSV (H [A/m], B [T]) or 'aedt:<library material name>'. The paper's "
                     "27SW1400 data are not shipped")
    pm_br: float = _f(1.40, "T", "illustrative",
                      "PM remanence (linear magnet). Paper: G52UH, 'higher than 1.40 T'")
    pm_mur: float = _f(1.05, "-", "illustrative", "PM recoil permeability")
    cu_conductivity: float = _f(5.8e7, "S/m", "illustrative", "Copper conductivity (AC-loss runs)")

    # ---- electromagnetic operating point / solver -----------------------------------------
    current_rms: float = _f(420.0, "A", "paper", "Phase current (rms); paper's full-load value")
    gamma: float = _f(130.0, "deg(e)", "illustrative",
                      "Current advance angle measured from the +d axis (90 = pure q-axis current)")
    em_rpm: float = _f(3000.0, "rpm", "illustrative", "Speed for torque runs (no eddy regions)")
    em_span: float = _f(60.0, "deg(e)", "illustrative", "Simulated electrical angle per torque run")
    em_step: float = _f(2.0, "deg(e)", "illustrative", "Time step expressed in electrical degrees")
    mesh_rotor: float = _f(1.2, "mm", "illustrative", "Max element length in rotor iron and magnets")
    mesh_stator: float = _f(1.5, "mm", "illustrative", "Max element length in stator iron and conductors")
    mesh_gap: float = _f(0.3, "mm", "illustrative", "Max element length in the air-gap layers")

    # ---- mechanical materials (paper Table II unless marked) --------------------------------
    lam_density: float = _f(7850.0, "kg/m3", "paper", "Lamination density")
    lam_E: float = _f(200.0, "GPa", "paper", "Lamination Young's modulus")
    lam_nu: float = _f(0.30, "-", "paper", "Lamination Poisson ratio")
    pm_density: float = _f(7400.0, "kg/m3", "paper", "PM density")
    pm_E: float = _f(160.0, "GPa", "paper", "PM Young's modulus")
    pm_nu: float = _f(0.24, "-", "paper", "PM Poisson ratio")
    sleeve_density: float = _f(1800.0, "kg/m3", "paper", "Sleeve density")
    sleeve_E_hoop: float = _f(125.0, "GPa", "paper",
                              "Sleeve modulus along the fibre (hoop). Table II lists 125 GPa as 'radial'; "
                              "see docs/methods.md")
    sleeve_E_radial: float = _f(8.8, "GPa", "paper", "Sleeve transverse (radial) modulus")
    sleeve_E_axial: float = _f(8.8, "GPa", "illustrative", "Sleeve axial modulus (transversely isotropic)")
    sleeve_nu_hoop_radial: float = _f(0.28, "-", "paper",
                                      "Major Poisson ratio -eps_r/eps_hoop under hoop stress. The minor "
                                      "ratio follows from symmetry (nu_hr*E_r/E_h)")
    sleeve_nu_hoop_axial: float = _f(0.28, "-", "illustrative", "Major Poisson ratio hoop-axial")
    sleeve_nu_radial_axial: float = _f(0.40, "-", "illustrative", "Transverse Poisson ratio radial-axial")
    sleeve_G_hoop_radial: float = _f(4.5, "GPa", "illustrative", "In-plane shear modulus")
    sleeve_G_hoop_axial: float = _f(4.5, "GPa", "illustrative", "Shear modulus hoop-axial")
    sleeve_G_radial_axial: float = _f(3.1, "GPa", "illustrative", "Transverse shear modulus")

    # ---- mechanical load case / solver -----------------------------------------------------
    mech_rpm: float = _f(20000.0, "rpm", "paper", "Speed for the stress analysis")
    interference: float = _f(0.10, "mm", "illustrative",
                             "Radial interference between lamination OD and sleeve bore")
    mu_sleeve: float = _f(0.0, "-", "illustrative", "Friction coefficient sleeve/lamination")
    mu_magnet: float = _f(0.0, "-", "illustrative", "Friction coefficient magnet/pocket")
    plane: str = _f("stress", "-", "illustrative", "2D assumption: 'stress' or 'strain'")
    mech_mesh: float = _f(0.5, "mm", "illustrative", "Global element size (mechanical)")
    mech_mesh_fine: float = _f(0.15, "mm", "illustrative",
                               "Element size on short lines (fillets, bridges, lugs)")
    magnet_contact: str = _f("bonded", "-", "illustrative",
                             "Magnet/pocket contact in the stress model: 'bonded' (as the original Workbench "
                             "model; the designed clearance is bridged by the bond, see docs/methods.md) or "
                             "'frictional' (open clearance, mu_magnet; not supported for sleeved_bridgeless)")
    stress_corner_margin: float = _f(0.5, "mm", "illustrative",
                                     "Named stress sections (svm_lam_max) leave out this length at both ends of "
                                     "each section, where it meets a sharp corner of the outline")

    # ------------------------------------------------------------------------------------
    def __post_init__(self):
        self.validate()

    def replace(self, **changes) -> Params:
        unknown = set(changes) - set(field_names())
        if unknown:
            raise KeyError(f"Unknown parameter(s): {sorted(unknown)}")
        return dataclasses.replace(self, **changes)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)

    # derived quantities -------------------------------------------------------------------
    @property
    def pole_pairs(self) -> int:
        return self.poles // 2

    @property
    def r_stator_in(self) -> float:
        return self.stator_id / 2

    @property
    def r_sleeve_out(self) -> float:
        return self.r_stator_in - self.airgap

    @property
    def r_lam(self) -> float:
        """Lamination outer radius: the air gap is fixed, so the sleeve eats into the rotor."""
        return self.r_sleeve_out - self.sleeve_thickness

    @property
    def sleeve_nu_radial_hoop(self) -> float:
        """Minor Poisson ratio, from compliance symmetry."""
        return self.sleeve_nu_hoop_radial * self.sleeve_E_radial / self.sleeve_E_hoop

    @property
    def is_bridgeless(self) -> bool:
        return self.topology == "sleeved_bridgeless"

    # validation ---------------------------------------------------------------------------
    def validate(self) -> None:
        err = []
        if self.topology in REMOVED_TOPOLOGIES:
            raise ValueError(f"topology {self.topology!r}: {REMOVED_TOPOLOGIES[self.topology]}")
        if self.topology not in TOPOLOGIES:
            err.append(f"topology must be one of {TOPOLOGIES}, got {self.topology!r}")
        for f in fields(self):
            v = getattr(self, f.name)
            if f.type in ("float", float) or isinstance(v, float):
                if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                    err.append(f"{f.name} must be a finite number, got {v!r}")
        nonneg = ["sleeve_thickness", "interference", "mu_sleeve", "mu_magnet", "current_rms",
                  "pm1_inner_offset", "pm2_inner_offset", "pm1_bridge_thickness", "pm2_bridge_thickness",
                  "pm1_end_gap", "pm2_end_gap", "pm_side_clearance", "pm_end_clearance", "pm_corner_radius",
                  "slot_fillet", "q_hole_fillet", "pm1_rib_half_width", "pm2_web_width", "min_web", "min_q_web",
                  "pm1_slot_start", "pm2_slot_start", "lip_depth"]
        positive = ["stator_od", "stator_id", "slot_opening", "airgap", "slot_tip_depth", "slot_tip_angle",
                    "slot_width", "slot_depth", "cond_width", "cond_height", "cond_start_radius",
                    "stack_length", "pm1_width", "pm2_width", "pm1_d_crossing", "pm2_d_crossing",
                    "pm1_length_ratio", "pm2_length_ratio",
                    "shaft_diameter", "q_hole_inner_radius", "q_hole_apex_radius", "q_hole_half_angle",
                    "d_hole_diameter", "d_hole_depth", "pm1_opening_angle", "pm2_opening_angle",
                    "pm_br", "pm_mur", "em_rpm",
                    "em_span", "em_step", "mesh_rotor", "mesh_stator", "mesh_gap", "lam_density",
                    "lam_E", "pm_density", "pm_E", "sleeve_density", "sleeve_E_hoop", "sleeve_E_radial",
                    "sleeve_E_axial", "sleeve_G_hoop_radial", "sleeve_G_hoop_axial",
                    "sleeve_G_radial_axial", "mech_mesh", "mech_mesh_fine", "cu_conductivity", "stress_corner_margin"]
        for n in nonneg:
            v = getattr(self, n)
            if isinstance(v, (int, float)) and v < 0:
                err.append(f"{n} must be >= 0, got {v}")
        for n in positive:
            v = getattr(self, n)
            if isinstance(v, (int, float)) and v <= 0:
                err.append(f"{n} must be > 0, got {v}")
        for n in ("slots", "poles", "hairpin_layers", "parallel_paths", "cond_count"):
            v = getattr(self, n)
            if not isinstance(v, int) or isinstance(v, bool) or v < 1:
                err.append(f"{n} must be a positive integer, got {v!r}")
        if not err:
            if self.poles % 2:
                err.append("poles must be even")
            elif self.slots % (3 * self.poles):
                err.append("only integer-slot windings (slots divisible by 3*poles) are supported")
            if self.stator_od <= self.stator_id:
                err.append("stator_od must exceed stator_id")
            for k in (1, 2):
                a = getattr(self, f"pm{k}_angle")
                if not 20.0 < a < 180.0:
                    err.append(f"pm{k}_angle must be in (20, 180) deg, got {a}")
            if not -45.0 < self.pm2_inner_tilt < 45.0:
                err.append("pm2_inner_tilt must be in (-45, 45) deg")
            if not 0 < self.slot_tip_angle < 90:
                err.append("slot_tip_angle must be in (0, 90) deg")
            if self.slot_width <= self.slot_opening:
                err.append("slot_width must exceed slot_opening")
            if self.sleeve_thickness <= 0:
                err.append(f"{self.topology} needs sleeve_thickness > 0")
            if self.is_bridgeless and self.lip_depth <= 0:
                err.append("sleeved_bridgeless needs lip_depth > 0")
            if self.pm_side_clearance >= min(self.pm1_width, self.pm2_width):
                err.append("pm_side_clearance must be smaller than the slot widths")
            if self.magnet_contact not in ("bonded", "frictional"):
                err.append(f"magnet_contact must be 'bonded' or 'frictional', got {self.magnet_contact!r}")
            if self.plane not in ("stress", "strain"):
                err.append("plane must be 'stress' or 'strain'")
            for n in ("lam_nu", "pm_nu", "sleeve_nu_hoop_radial", "sleeve_nu_hoop_axial",
                      "sleeve_nu_radial_axial"):
                if not 0 <= getattr(self, n) < 0.5:
                    err.append(f"{n} must be in [0, 0.5)")
            if self.sleeve_nu_radial_hoop * self.sleeve_nu_hoop_radial >= 1:
                err.append("sleeve Poisson ratios give a non-positive-definite compliance")
            if self.r_lam <= self.shaft_diameter / 2:
                err.append("sleeve too thick: lamination OD is inside the bore")
            if self.em_span < self.em_step:
                err.append("em_span must be >= em_step")
        if err:
            raise ValueError("Invalid parameters:\n  - " + "\n  - ".join(err))


# topology-specific defaults (values differing from the dataclass defaults) -----------------
# Each set is the saved state of the corresponding Maxwell design.
TOPOLOGY_DEFAULTS: dict[str, dict[str, Any]] = {
    "sleeved_bridged": dict(
        # model saved at Wcf = 2 mm; the paper's bridged design point is 1.5 mm
        sleeve_thickness=1.5,
    ),
    "sleeved_bridgeless": dict(
        sleeve_thickness=0.5,
        pm1_rib_half_width=0.1, pm2_web_width=0.1, pm2_inner_offset=0.0,
        pm1_bridge_thickness=0.5, pm2_bridge_thickness=0.5,
        pm_end_clearance=0.0,
    ),
}

# Fields whose value or source differs per topology.
TOPOLOGY_SOURCE_NOTES = {
    "sleeved_bridged": "saved state of 01_ICEM_01_simplerotor except the sleeve: 1.5 mm as in the paper's "
                       "Section III (model saved at 2.0 mm)",
    "sleeved_bridgeless": "saved state of 02_ICEM_01_teslamotor1 (sleeve 0.5 mm; the paper's bridgeless "
                          "comparison uses 2.0 mm)",
}


def default_params(topology: str = "sleeved_bridged", **overrides) -> Params:
    """Defaults for a topology, optionally with overrides (validated)."""
    if topology in REMOVED_TOPOLOGIES:
        raise ValueError(f"topology {topology!r}: {REMOVED_TOPOLOGIES[topology]}")
    if topology not in TOPOLOGIES:
        raise ValueError(f"topology must be one of {TOPOLOGIES}, got {topology!r}")
    kw = dict(TOPOLOGY_DEFAULTS[topology])
    kw.update(overrides)
    unknown = set(kw) - set(field_names())
    if unknown:
        raise KeyError(f"Unknown parameter(s): {sorted(unknown)}")
    return Params(topology=topology, **kw)


def field_names() -> list[str]:
    return [f.name for f in fields(Params)]


def coerce(name: str, value: Any) -> Any:
    """Convert a value read from a CSV/TOML table to the field's type."""
    f = {f.name: f for f in fields(Params)}[name]
    default = f.default
    if isinstance(default, bool):
        return bool(value)
    if isinstance(default, int):
        fv = float(value)
        if fv != int(fv):
            raise ValueError(f"{name} must be an integer, got {value!r}")
        return int(fv)
    if isinstance(default, float):
        return float(value)
    return str(value)


def parameter_table() -> list[dict[str, Any]]:
    rows = []
    for f in fields(Params):
        rows.append(dict(name=f.name, default=f.default, unit=f.metadata.get("unit", ""),
                         source=f.metadata.get("source", ""), doc=f.metadata.get("doc", ""),
                         maxwell_name=f.metadata.get("maxwell_name", "")))
    return rows


def parameter_markdown() -> str:
    lines = ["| Parameter | Default | Unit | Source | Maxwell variable | Meaning |", "|---|---|---|---|---|---|"]
    for r in parameter_table():
        mx = f"`{r['maxwell_name']}`" if r["maxwell_name"] else ""
        lines.append(f"| `{r['name']}` | {r['default']} | {r['unit']} | {r['source']} | {mx} | {r['doc']} |")
    lines.append("")
    lines.append("Topology-specific overrides of the defaults above:")
    lines.append("")
    for t, kw in TOPOLOGY_DEFAULTS.items():
        kv = ", ".join(f"`{k}`={v}" for k, v in kw.items())
        lines.append(f"* **{t}**: {kv} ({TOPOLOGY_SOURCE_NOTES[t]})")
    return "\n".join(lines) + "\n"
