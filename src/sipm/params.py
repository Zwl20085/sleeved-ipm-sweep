"""Machine, material and analysis parameters.

All lengths are in mm, angles in degrees, moduli in GPa, densities in kg/m^3,
speeds in rpm, currents in A, temperatures in degC and thermal expansion coefficients in
1/K unless the field name says otherwise.

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
CTE_FIELDS = ("lam_cte", "pm_cte", "sleeve_cte_hoop", "sleeve_cte_radial", "sleeve_cte_axial")
CTE_RANGE = (-2.0e-5, 1.0e-4)             # 1/K; engineering solids lie well inside
TEMPERATURE_RANGE = (-273.15, 500.0)      # degC, (exclusive, inclusive]
# Mechanical load cases of ``mech_cases``, in their canonical order (the order in which they are
# solved and reported, and which breaks ties of the envelope): speed x temperature.
MECH_CASES = ("standstill_cold", "standstill_hot", "rated_cold", "rated_hot", "overspeed_cold", "overspeed_hot")
CASE_SEPARATOR = "__"                     # per-case result column: <output>__<case>
# layered sleeve (sleeve_layers): limits that keep the deck and its mesh sane
MAX_SLEEVE_LAYERS = 10
MIN_LAYER_THICKNESS = 0.05                # mm, thinnest layer accepted
SLEEVE_SIGMA_MAX = 1.0e8                  # S/m; sleeve conductivities above copper are rejected as unit errors
# string fields that a TOML table may give as an array (joined with ',')
LIST_STRING_FIELDS = ("mech_cases", "sleeve_layers", "sleeve_layer_prestress")


def parse_mech_cases(spec: str) -> tuple[str, ...]:
    """Case names of a ``mech_cases`` string, in the order given ('' -> ()).

    Comma-separated, surrounding blanks ignored; every name must be one of :data:`MECH_CASES`
    and appear once."""
    if not isinstance(spec, str):
        raise ValueError(f"mech_cases must be a comma-separated string, got {spec!r}")
    if not spec.strip():
        return ()
    names = tuple(s.strip() for s in spec.split(","))
    bad = [n for n in names if n not in MECH_CASES]
    if bad:
        raise ValueError(f"mech_cases: unknown case(s) {bad}; allowed: {', '.join(MECH_CASES)}")
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise ValueError(f"mech_cases: case(s) listed more than once: {dup}")
    return names


def _number(text: str, what: str) -> float:
    try:
        v = float(text)
    except ValueError:
        raise ValueError(f"{what}: {text.strip()!r} is not a number") from None
    if not math.isfinite(v):
        raise ValueError(f"{what}: {text.strip()!r} is not finite")
    return v


def parse_sleeve_layers(spec: str) -> tuple[tuple[float, float], ...]:
    """(angle in deg, weight) of each layer of a ``sleeve_layers`` string, from the bore outward
    ('' -> ()).

    Comma-separated ``angle:weight`` items, blanks ignored; angle = balanced +/- winding angle
    from the hoop direction in [0, 90], weight > 0 (normalised later to fractions of the sleeve
    thickness); at most :data:`MAX_SLEEVE_LAYERS` items."""
    if not isinstance(spec, str):
        raise ValueError(f"sleeve_layers must be a string of 'angle:weight' items, got {spec!r}")
    if not spec.strip():
        return ()
    out = []
    for item in spec.split(","):
        parts = item.split(":")
        if len(parts) != 2:
            raise ValueError(f"sleeve_layers: item {item.strip()!r} is not of the form 'angle:weight'")
        ang, wt = _number(parts[0], "sleeve_layers angle"), _number(parts[1], "sleeve_layers weight")
        if not 0.0 <= ang <= 90.0:
            raise ValueError(f"sleeve_layers: angle {ang:g} deg is outside [0, 90]")
        if not wt > 0.0:
            raise ValueError(f"sleeve_layers: weight {wt:g} must be > 0")
        out.append((ang, wt))
    if len(out) > MAX_SLEEVE_LAYERS:
        raise ValueError(f"sleeve_layers: {len(out)} layers, at most {MAX_SLEEVE_LAYERS} are supported")
    return tuple(out)


def parse_layer_prestress(spec: str) -> tuple[float, ...]:
    """Values (MPa) of a ``sleeve_layer_prestress`` string ('' -> ()); comma-separated, each finite
    and >= 0."""
    if not isinstance(spec, str):
        raise ValueError(f"sleeve_layer_prestress must be a string of comma-separated values, got {spec!r}")
    if not spec.strip():
        return ()
    vals = tuple(_number(s, "sleeve_layer_prestress") for s in spec.split(","))
    neg = [v for v in vals if v < 0]
    if neg:
        raise ValueError(f"sleeve_layer_prestress: values must be >= 0 (winding tension), got {neg}")
    return vals


def layer_radii(r_in: float, r_out: float, fractions) -> tuple[float, ...]:
    """Radii of the layer boundaries, bore to OD (len(fractions) + 1 values; the last is ``r_out``
    exactly)."""
    t = r_out - r_in
    out, acc = [float(r_in)], 0.0
    for f in fractions[:-1]:
        acc += f
        out.append(r_in + acc * t)
    out.append(float(r_out))
    return tuple(out)


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
    sleeve_sigma_fibre: float = _f(0.0, "S/m", "illustrative",
                                   "Sleeve ply electrical conductivity ALONG the fibre (sigma_L). With "
                                   "sleeve_sigma_transverse it gives the axial conductivity of the sleeve "
                                   "(sleeve_conductivity_axial = sigma_L sin^2 t + sigma_T cos^2 t per layer at "
                                   "winding angle t from hoop, layers thickness-weighted), the only one a 2D model "
                                   "with axial current uses. 0 and 0 (default) = non-conducting sleeve, model and "
                                   "results unchanged; > 0 axial conductivity adds the 2D sleeve eddy-current loss "
                                   "to EM runs (docs/methods.md, 2.1). Illustrative class statement only, not a "
                                   "material property: carbon/epoxy is often quoted at 1e3-1e5 S/m along the fibre "
                                   "and orders of magnitude lower across it; measure your laminate")
    sleeve_sigma_transverse: float = _f(0.0, "S/m", "illustrative",
                                        "Sleeve ply electrical conductivity ACROSS the fibre (sigma_T), set by "
                                        "fibre-to-fibre contacts through the matrix; for a hoop-wound sleeve "
                                        "(winding angle 0) it is the axial conductivity. Illustrative class "
                                        "statement only: typically several orders of magnitude below sigma_L and "
                                        "strongly process-dependent; measure your laminate")

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

    # ---- composite sleeve specification (optional; sipm.composite) ---------------------------
    sleeve_fibre: str = _f("", "-", "illustrative",
                           "Sleeve fibre class from sipm.composite.FIBRES (carbon_sm, carbon_im, carbon_hm, "
                           "e_glass), or '' (default) to use the sleeve constants as entered. When set, the "
                           "thirteen sleeve constants (sleeve_E_*, sleeve_nu_*, sleeve_G_*, sleeve_density, "
                           "sleeve_cte_*) are DERIVED from fibre, matrix, vf, micromechanics model and winding "
                           "angle whenever the parameters are created or replaced, overwriting any values "
                           "passed for them")
    sleeve_matrix: str = _f("epoxy", "-", "illustrative",
                            "Sleeve matrix class from sipm.composite.MATRICES (epoxy, peek); used only with "
                            "sleeve_fibre")
    sleeve_vf: float = _f(0.55, "-", "illustrative", "Sleeve fibre volume fraction, in (0, 1); used only with "
                          "sleeve_fibre")
    sleeve_micromechanics: str = _f("chamis", "-", "illustrative",
                                    "Micromechanics model for the transverse and shear ply constants: chamis, "
                                    "halpin_tsai or rom; used only with sleeve_fibre")
    sleeve_wind_angle: float = _f(0.0, "deg", "illustrative",
                                  "Winding angle of a balanced +/- angle-ply sleeve, measured from the hoop "
                                  "direction (0 = hoop winding), in [0, 90]. Non-zero needs sleeve_fibre (the "
                                  "ply cannot be recovered from homogenised constants)")
    sleeve_layers: str = _f("", "-", "illustrative",
                            "Layered sleeve: comma-separated 'angle:weight' items from the bore outward, e.g. "
                            "'0:2,20:1,0:1' (angle = balanced +/- winding angle from hoop in [0, 90] deg, "
                            "weights normalised to fractions of sleeve_thickness); '' (default) = one homogeneous "
                            "layer at sleeve_wind_angle. Needs sleeve_fibre (one fibre, matrix and vf for all "
                            "layers) and sleeve_wind_angle = 0; at most 10 layers, each at least 0.05 mm thick. "
                            "In the deck each layer is its own area and material, perfectly bonded (shared "
                            "nodes). The thirteen sleeve constants then hold the membrane homogenisation of the "
                            "stack (layers share the hoop and axial strain and the radial stress; "
                            "sipm.composite.layer_stack), i.e. the stiffness and free thermal strain of an "
                            "equivalent thin homogeneous sleeve: for inspection and the cache key only, the deck "
                            "and the failure evaluation use each layer's own constants (docs/methods.md, 4.6)")
    # unidirectional ply strengths for the sleeve failure indices (sipm.composite.failure_indices)
    sleeve_Xt: float = _f(2000.0, "MPa", "illustrative",
                          "Sleeve ply strength, fibre-direction tension (default: STRENGTHS['carbon_sm_epoxy'], "
                          "a class-typical round number)")
    sleeve_Xc: float = _f(1200.0, "MPa", "illustrative", "Sleeve ply strength, fibre-direction compression")
    sleeve_Yt: float = _f(50.0, "MPa", "illustrative", "Sleeve ply strength, transverse tension")
    sleeve_Yc: float = _f(200.0, "MPa", "illustrative", "Sleeve ply strength, transverse compression")
    sleeve_S: float = _f(70.0, "MPa", "illustrative",
                         "Sleeve ply shear strength (also used for the transverse shear strength in Hashin's "
                         "matrix-compression mode)")

    # ---- thermal expansion (secant CTE about mech_ref_temperature) -----------------------------
    lam_cte: float = _f(12.0e-6, "1/K", "illustrative",
                        "Lamination thermal expansion coefficient (electrical steel, typically 11-13e-6 /K)")
    pm_cte: float = _f(2.0e-6, "1/K", "illustrative",
                       "PM thermal expansion coefficient, one isotropic value. Sintered NdFeB is strongly "
                       "anisotropic (about +5e-6 /K along the magnetisation, about -1e-6 /K across it); the "
                       "default is the mean of the two in-plane directions")
    sleeve_cte_hoop: float = _f(0.0, "1/K", "illustrative",
                                "Sleeve thermal expansion along the fibre (hoop); carbon/epoxy is about zero "
                                "(typically -1e-6 .. +0.5e-6 /K)")
    sleeve_cte_radial: float = _f(30.0e-6, "1/K", "illustrative",
                                  "Sleeve transverse (radial) thermal expansion, matrix-dominated (typically "
                                  "25-35e-6 /K for carbon/epoxy)")
    sleeve_cte_axial: float = _f(30.0e-6, "1/K", "illustrative",
                                 "Sleeve axial thermal expansion (transversely isotropic: = radial); enters plane "
                                 "strain only")

    # ---- mechanical load case / solver -----------------------------------------------------
    mech_rpm: float = _f(20000.0, "rpm", "paper", "Speed for the stress analysis")
    interference: float = _f(0.10, "mm", "illustrative",
                             "Radial interference between lamination OD and sleeve bore")
    sleeve_layer_prestress: str = _f("", "MPa", "illustrative",
                                     "Winding prestress of the sleeve layers: one value per layer (sleeve_layers "
                                     "order, from the bore) or one value for all, each >= 0 ('' = none). Value "
                                     "sigma0 is a hoop eigenstrain: the layer's stress-free hoop strain is "
                                     "-sigma0/E_hoop of that layer (sigma0 = the hoop stress the layer carries "
                                     "when held at its fitted radius with no radial stress), applied to all "
                                     "layers at once before elastic equilibrium (MAPDL INISTATE initial elastic "
                                     "strain); the outputs are the relaxed state. Acts together with "
                                     "interference (which may be 0), temperature and speed. A simultaneous-"
                                     "eigenstrain model, not the layer-by-layer build-up of real winding; "
                                     "`sipm winding` converts a winding tension into equivalent values "
                                     "(docs/methods.md, 4.6)")
    mech_temperature: float = _f(20.0, "degC", "illustrative",
                                 "Uniform rotor temperature (lamination, magnets and sleeve) of the stress "
                                 "analysis")
    mech_ref_temperature: float = _f(20.0, "degC", "illustrative",
                                     "Reference temperature: the temperature at which the parts have their "
                                     "stress-free fitted dimensions, i.e. at which interference is defined. "
                                     "Equal to mech_temperature (default) = isothermal analysis")
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

    # ---- mechanical load-case set (sipm.mapdl.solve_mech_cases) ----------------------------
    mech_cases: str = _f("", "-", "illustrative",
                         "Comma-separated mechanical load cases solved for each point, from standstill_cold, "
                         "standstill_hot, rated_cold, rated_hot, overspeed_cold, overspeed_hot ('' = the single "
                         "case mech_rpm at mech_temperature, as before). Speeds: standstill 0, rated mech_rpm, "
                         "overspeed overspeed_factor x mech_rpm; temperatures: hot mech_temperature, cold "
                         "mech_temperature_cold. Each output gets one column per case (e.g. "
                         "sh_sleeve_max__overspeed_hot), plus worst-case envelope columns and the lift-off speeds "
                         "(docs/methods.md, 3.1)")
    overspeed_factor: float = _f(1.2, "-", "illustrative",
                                 "Overspeed load cases: speed = this factor x mech_rpm (>= 1)")
    mech_temperature_cold: float = _f(20.0, "degC", "illustrative",
                                      "Uniform rotor temperature of the '_cold' load cases (the '_hot' cases use "
                                      "mech_temperature); interference is defined at mech_ref_temperature")
    liftoff_search_tol: float = _f(100.0, "rpm", "illustrative",
                                   "Lift-off speed of rotors whose sleeve contact is not linear in speed "
                                   "(sleeved_bridgeless): found with extra MAPDL runs at trial speeds until "
                                   "the bracket (contact pressure positive everywhere / zero somewhere) is this "
                                   "narrow; the value is its midpoint. 0 = no search, lift-off NaN "
                                   "(docs/methods.md, 3.1)")

    # ------------------------------------------------------------------------------------
    def __post_init__(self):
        self._derive_sleeve()
        self.validate()

    def _composite_spec_ok(self) -> bool:
        from . import composite

        try:
            parse_sleeve_layers(self.sleeve_layers)
        except ValueError:
            return False
        return (self.sleeve_fibre in composite.FIBRES and self.sleeve_matrix in composite.MATRICES
                and self.sleeve_micromechanics in composite.MODELS
                and isinstance(self.sleeve_vf, (int, float)) and 0.0 < self.sleeve_vf < 1.0
                and isinstance(self.sleeve_wind_angle, (int, float)) and 0.0 <= self.sleeve_wind_angle <= 90.0)

    def _derive_sleeve(self) -> None:
        """With a composite spec, overwrite the thirteen sleeve constants by the derived ones (an
        invalid spec is left for :meth:`validate` to report). A layered sleeve gives the membrane
        homogenisation of its layers (:func:`sipm.composite.layer_stack`)."""
        if self.sleeve_fibre and self._composite_spec_ok():
            from .composite import sleeve_params

            for k, v in sleeve_params(self.sleeve_laminate()).items():
                setattr(self, k, v)

    def _ply(self):
        from .composite import FIBRES, MATRICES, lamina

        return lamina(FIBRES[self.sleeve_fibre], MATRICES[self.sleeve_matrix], self.sleeve_vf,
                      self.sleeve_micromechanics)

    def sleeve_laminate(self):
        """The sleeve as a :class:`sipm.composite.Laminate` (balanced +/- ``sleeve_wind_angle``
        winding) from the composite spec, or None when the constants are entered directly
        (``sleeve_fibre`` = ''). With ``sleeve_layers`` set it is the
        :class:`sipm.composite.LayerStack` of the layers (the source of the thirteen constants; it
        has no plies, the layers are :meth:`sleeve_layer_laminates`)."""
        if not self.sleeve_fibre:
            return None
        from .composite import angle_ply, layer_stack

        if self.is_layered:
            return layer_stack(self.sleeve_layer_laminates(), self.sleeve_layer_fractions())
        return angle_ply(self._ply(), self.sleeve_wind_angle)

    # layers (a single-layer sleeve is one layer) ----------------------------------------------
    @property
    def is_layered(self) -> bool:
        """``sleeve_layers`` is set (each layer its own area and material in the deck)."""
        return bool(self.sleeve_layers.strip())

    @property
    def sleeve_layer_count(self) -> int:
        return len(parse_sleeve_layers(self.sleeve_layers)) or 1

    def sleeve_layer_angles(self) -> tuple[float, ...]:
        """Winding angle of each layer from the bore outward (one layer: ``sleeve_wind_angle``)."""
        layers = parse_sleeve_layers(self.sleeve_layers)
        return tuple(a for a, _ in layers) if layers else (float(self.sleeve_wind_angle),)

    def sleeve_layer_fractions(self) -> tuple[float, ...]:
        """Thickness fraction of each layer from the bore outward (sum 1)."""
        layers = parse_sleeve_layers(self.sleeve_layers)
        if not layers:
            return (1.0,)
        tot = sum(w for _, w in layers)
        return tuple(w / tot for _, w in layers)

    def sleeve_layer_radii(self) -> tuple[float, ...]:
        """Layer boundary radii of the rotor sleeve, mm (``r_lam`` .. ``r_sleeve_out``)."""
        return layer_radii(self.r_lam, self.r_sleeve_out, self.sleeve_layer_fractions())

    def sleeve_layer_laminates(self) -> list | None:
        """One :class:`sipm.composite.Laminate` per layer from the bore outward (a single-layer sleeve:
        ``[sleeve_laminate()]``), or None when the constants are entered directly."""
        if not self.sleeve_fibre:
            return None
        from .composite import angle_ply

        ply = self._ply()
        return [angle_ply(ply, a) for a in self.sleeve_layer_angles()]

    def sleeve_layer_constants(self) -> list[dict[str, float]]:
        """The thirteen sleeve constants (:data:`sipm.composite.SLEEVE_KEYS`) of each layer; with
        directly entered constants, the one layer's values from the fields."""
        from .composite import SLEEVE_KEYS, sleeve_params

        lams = self.sleeve_layer_laminates()
        if lams is None:
            return [{k: getattr(self, k) for k in SLEEVE_KEYS}]
        return [sleeve_params(lam) for lam in lams]

    def sleeve_layer_prestress_mpa(self) -> tuple[float, ...]:
        """Winding prestress of each layer, MPa (zeros when ``sleeve_layer_prestress`` is '')."""
        vals, n = parse_layer_prestress(self.sleeve_layer_prestress), self.sleeve_layer_count
        if not vals:
            return (0.0,) * n
        return vals * n if len(vals) == 1 else vals

    def sleeve_layer_prestrain(self) -> tuple[float, ...]:
        """Hoop eigenstrain magnitude of each layer, eps0 = sigma0 / E_hoop (the layer's stress-free
        hoop strain is -eps0)."""
        return tuple(s / (c["sleeve_E_hoop"] * 1e3)
                     for s, c in zip(self.sleeve_layer_prestress_mpa(), self.sleeve_layer_constants()))

    def sleeve_strength(self):
        """The sleeve ply strengths as a :class:`sipm.composite.Strength` (St = S)."""
        from .composite import Strength

        return Strength(Xt=self.sleeve_Xt, Xc=self.sleeve_Xc, Yt=self.sleeve_Yt, Yc=self.sleeve_Yc, S=self.sleeve_S)

    def sleeve_compliance(self):
        """6x6 sleeve compliance [1/GPa] from the sleeve constants, cylinder Voigt order (hoop,
        axial, radial, axial-radial, hoop-radial, hoop-axial), engineering shear strains,
        S[a, b] = -nu_ab / E_a with nu_ab = -eps_b/eps_a under stress along a."""
        import numpy as np

        Eh, Ea, Er = self.sleeve_E_hoop, self.sleeve_E_axial, self.sleeve_E_radial
        S = np.zeros((6, 6))
        S[0, 0], S[1, 1], S[2, 2] = 1 / Eh, 1 / Ea, 1 / Er
        S[0, 1] = S[1, 0] = -self.sleeve_nu_hoop_axial / Eh
        S[0, 2] = S[2, 0] = -self.sleeve_nu_hoop_radial / Eh
        S[2, 1] = S[1, 2] = -self.sleeve_nu_radial_axial / Er
        S[3, 3] = 1 / self.sleeve_G_radial_axial
        S[4, 4] = 1 / self.sleeve_G_hoop_radial
        S[5, 5] = 1 / self.sleeve_G_hoop_axial
        return S

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
    def sleeve_conductivity_axial(self) -> float:
        """Axial electrical conductivity of the sleeve, S/m: each layer's ply value along z,
        sigma_fibre sin^2 t + sigma_transverse cos^2 t at its winding angle t from hoop, thickness-
        weighted over the layers (:func:`sipm.composite.axial_conductivity`). The EM model sees the
        sleeve only through this number."""
        from .composite import axial_conductivity

        return axial_conductivity(self.sleeve_sigma_fibre, self.sleeve_sigma_transverse,
                                  self.sleeve_layer_angles(), self.sleeve_layer_fractions())

    @property
    def sleeve_nu_radial_hoop(self) -> float:
        """Minor Poisson ratio, from compliance symmetry."""
        return self.sleeve_nu_hoop_radial * self.sleeve_E_radial / self.sleeve_E_hoop

    @property
    def mech_delta_t(self) -> float:
        """Uniform temperature change of the stress analysis from the fitted (reference) state, K."""
        return self.mech_temperature - self.mech_ref_temperature

    @property
    def is_bridgeless(self) -> bool:
        return self.topology == "sleeved_bridgeless"

    @property
    def mech_case_list(self) -> tuple[str, ...]:
        """The load cases of ``mech_cases`` in the order given (() = the single default case)."""
        return parse_mech_cases(self.mech_cases)

    def case_speed_temperature(self, case: str) -> tuple[float, float]:
        """(speed in rpm, temperature in degC) of a load case of :data:`MECH_CASES`."""
        if case not in MECH_CASES:
            raise ValueError(f"unknown load case {case!r}; allowed: {', '.join(MECH_CASES)}")
        speed, temp = case.split("_")
        rpm = {"standstill": 0.0, "rated": self.mech_rpm, "overspeed": self.overspeed_factor * self.mech_rpm}[speed]
        return rpm, (self.mech_temperature if temp == "hot" else self.mech_temperature_cold)

    def for_case(self, case: str) -> Params:
        """Single-case parameters of one load case: ``mech_rpm`` and ``mech_temperature`` set to the
        case's speed and temperature, ``mech_cases`` cleared."""
        rpm, temp = self.case_speed_temperature(case)
        return self.replace(mech_cases="", mech_rpm=rpm, mech_temperature=temp)

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
                  "pm1_slot_start", "pm2_slot_start", "lip_depth", "sleeve_sigma_fibre", "sleeve_sigma_transverse"]
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
                    "sleeve_G_radial_axial", "mech_mesh", "mech_mesh_fine", "cu_conductivity", "stress_corner_margin",
                    "sleeve_Xt", "sleeve_Xc", "sleeve_Yt", "sleeve_Yc", "sleeve_S"]
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
            for n in ("lam_nu", "pm_nu"):
                if not 0 <= getattr(self, n) < 0.5:
                    err.append(f"{n} must be in [0, 0.5)")
            # sleeve: the full orthotropic compliance must be symmetric positive definite (the
            # condition MAPDL checks too); a +/- angle-ply winding legitimately has Poisson ratios
            # outside [0, 0.5), so the ratios are not bounded individually
            import numpy as np

            if np.linalg.eigvalsh(self.sleeve_compliance()).min() <= 0:
                err.append("sleeve elastic constants give a compliance that is not positive definite "
                           "(check the Poisson ratios against the moduli)")
            if self.sleeve_fibre:
                from . import composite

                if self.sleeve_fibre not in composite.FIBRES:
                    err.append(f"sleeve_fibre must be '' or one of {tuple(composite.FIBRES)}, "
                               f"got {self.sleeve_fibre!r}")
                if self.sleeve_matrix not in composite.MATRICES:
                    err.append(f"sleeve_matrix must be one of {tuple(composite.MATRICES)}, got {self.sleeve_matrix!r}")
                if self.sleeve_micromechanics not in composite.MODELS:
                    err.append(f"sleeve_micromechanics must be one of {composite.MODELS}, "
                               f"got {self.sleeve_micromechanics!r}")
            elif self.sleeve_wind_angle != 0:
                err.append("sleeve_wind_angle != 0 needs a composite spec (sleeve_fibre): the ply cannot be "
                           "recovered from homogenised constants")
            try:
                layers = parse_sleeve_layers(self.sleeve_layers)
            except ValueError as e:
                err.append(str(e))
                layers = None
            if layers:
                if not self.sleeve_fibre:
                    err.append("sleeve_layers needs a composite spec (sleeve_fibre): the layer constants are "
                               "derived from fibre, matrix and vf")
                if self.sleeve_wind_angle != 0:
                    err.append("with sleeve_layers set, sleeve_wind_angle must be 0 (each layer's angle is given "
                               f"in sleeve_layers), got {self.sleeve_wind_angle}")
                tot = sum(w for _, w in layers)
                thin = min(w / tot for _, w in layers) * self.sleeve_thickness
                if thin < MIN_LAYER_THICKNESS:
                    err.append(f"sleeve_layers: the thinnest layer is {thin:.4g} mm, below the {MIN_LAYER_THICKNESS} "
                               "mm minimum (fewer layers, other weights or a thicker sleeve)")
            try:
                pre = parse_layer_prestress(self.sleeve_layer_prestress)
            except ValueError as e:
                err.append(str(e))
                pre = ()
            if pre and layers is not None:
                n = len(layers) or 1
                if len(pre) not in (1, n):
                    err.append(f"sleeve_layer_prestress: {len(pre)} values for {n} sleeve layer(s); give one value "
                               "per layer or one for all")
            if not 0.0 < self.sleeve_vf < 1.0:
                err.append(f"sleeve_vf must be in (0, 1), got {self.sleeve_vf}")
            if not 0.0 <= self.sleeve_wind_angle <= 90.0:
                err.append(f"sleeve_wind_angle must be in [0, 90] deg, got {self.sleeve_wind_angle}")
            for n in CTE_FIELDS:
                if not CTE_RANGE[0] <= getattr(self, n) <= CTE_RANGE[1]:
                    err.append(f"{n} must be in [{CTE_RANGE[0]:g}, {CTE_RANGE[1]:g}] 1/K (a value in "
                               f"1e-6/K entered without the factor?), got {getattr(self, n)}")
            try:
                parse_mech_cases(self.mech_cases)
            except ValueError as e:
                err.append(str(e))
            if self.overspeed_factor < 1.0:
                err.append(f"overspeed_factor must be >= 1, got {self.overspeed_factor}")
            if self.liftoff_search_tol < 0:
                err.append(f"liftoff_search_tol must be >= 0 rpm (0 = no search), got {self.liftoff_search_tol}")
            for n in ("mech_temperature", "mech_ref_temperature", "mech_temperature_cold"):
                if not TEMPERATURE_RANGE[0] < getattr(self, n) <= TEMPERATURE_RANGE[1]:
                    err.append(f"{n} must be in ({TEMPERATURE_RANGE[0]:g}, {TEMPERATURE_RANGE[1]:g}] degC, "
                               f"got {getattr(self, n)}")
            if self.r_lam <= self.shaft_diameter / 2:
                err.append("sleeve too thick: lamination OD is inside the bore")
            if self.em_span < self.em_step:
                err.append("em_span must be >= em_step")
            for n in ("sleeve_sigma_fibre", "sleeve_sigma_transverse"):
                if getattr(self, n) > SLEEVE_SIGMA_MAX:
                    err.append(f"{n} = {getattr(self, n):g} S/m exceeds {SLEEVE_SIGMA_MAX:g} S/m (above copper; "
                               "wrong unit?)")
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
    if name in LIST_STRING_FIELDS and isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)   # TOML may give the cases / layers / prestresses as an array
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
