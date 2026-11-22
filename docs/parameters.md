# Parameters

Generated with `sipm params --markdown` (constants section below written by hand from
`sipm.geometry.rotor`). Source tags: **maxwell** = value of the author's original Maxwell 2D
models at their saved state (the "Maxwell variable" column names the model variable or
operation); **paper** = value published in the companion paper (Tables I and II, Section III);
**illustrative** = default chosen for this toolchain, not a dimension of the paper's machine.

The rotor geometry follows the Maxwell construction (docs/methods.md, section 1): each magnet
layer is fixed by its d-axis crossing, angle and width, and the **magnet length is derived**
(`rp.metrics["pm1_length"]`, `["pm2_length"]`), so it is not a parameter. Two names do not
measure what similar names measured before the rebase:

* `pm1_rib_half_width` (`RD_FB_lay1`) is **half** the layer-1 centre rib: the slot's inner end
  is the line y = this value, so the rib is twice as wide.
* `pm2_web_width` (`RD_FB_lay2`) is the thickness of each of the two webs between the tilted
  inner end of the layer-2 slot and the central cavity; there is no single layer-2 centre rib.

| Parameter | Default | Unit | Source | Maxwell variable | Meaning |
|---|---|---|---|---|---|
| `topology` | sleeved_bridged | - | paper |  | Rotor topology: sleeved_bridged or sleeved_bridgeless |
| `stator_od` | 210.0 | mm | maxwell | `2*Rso` | Stator outer diameter (model Rso = 105 mm; the paper's Table I lists 220 mm) |
| `stator_id` | 147.84 | mm | maxwell | `2*Rsi` | Stator inner (bore) diameter (paper Table I agrees) |
| `slots` | 48 | - | paper |  | Number of stator slots |
| `poles` | 8 | - | paper |  | Number of rotor poles |
| `slot_opening` | 1.9 | mm | maxwell | `Wslot_open` | Stator slot opening width (paper Table I agrees) |
| `slot_tip_depth` | 1.22 | mm | maxwell | `Wtip` | Radial length of the parallel-sided slot opening, measured from the point where the opening side meets the bore |
| `slot_tip_angle` | 70.0 | deg | maxwell | `tip_angle` | Angle between the slot centre line and the tapered tooth-tip side |
| `slot_width` | 5.27 | mm | maxwell | `Wslot` | Width of the parallel-sided slot body |
| `slot_depth` | 15.75 | mm | maxwell | `Hslot` | Slot bottom radius minus bore radius (flat bottom at Rsi + slot_depth) |
| `slot_fillet` | 0.3 | mm | maxwell |  | Fillet radius at every slot corner (incl. the bore corners) |
| `airgap` | 1.32 | mm | maxwell | `g - Wcf` | Physical single-side air gap (sleeve OD to stator bore); held fixed (paper agrees) |
| `cond_width` | 4.3 | mm | maxwell | `Wcoil` | Conductor width (tangential) |
| `cond_height` | 1.57 | mm | maxwell | `Lcoil` | Conductor height (radial) |
| `cond_start_radius` | 76.8 | mm | maxwell | `Rcoil_start` | Radius of the slot-opening side of the innermost conductor |
| `cond_count` | 8 | - | maxwell | `DuplicateAlongLine NumClones` | Conductor rectangles drawn per slot (touching, merged into one coil region per slot) |
| `hairpin_layers` | 7 | - | maxwell | `Nac` | Conductors (turns) per slot assigned to the coil (paper Table I agrees) |
| `stack_length` | 130.0 | mm | maxwell | `lstk` | Active stack length (2D model depth) |
| `parallel_paths` | 2 | - | maxwell |  | Parallel paths per phase |
| `pm1_d_crossing` | 56.7305951647269 | mm | maxwell | `RD_xpm_lay1` | Layer 1 (outer): distance from the shaft axis to the point where the magnet centre line crosses the d-axis; fixes the radial position of the V |
| `pm2_d_crossing` | 45.0864009239181 | mm | maxwell | `RD_xpm_lay2` | As pm1_d_crossing for layer 2 (inner) |
| `pm1_angle` | 97.007958474619 | deg | maxwell | `2*RD_thetapm_lay1` | Layer-1 included V angle = 2 x angle between magnet axis and d-axis (paper: 97) |
| `pm2_angle` | 95.0289047910392 | deg | maxwell | `2*RD_thetapm_lay2` | Layer-2 included V angle (paper: 95.02) |
| `pm1_width` | 3.64900294777709 | mm | maxwell | `RD_Wpm_lay1` | Layer-1 slot width = nominal magnet thickness along magnetisation (paper: 3.65) |
| `pm2_width` | 4.30590729759734 | mm | maxwell | `RD_Wpm_lay2` | Layer-2 slot width (paper: 4.31) |
| `pm1_end_gap` | 1.0 | mm | maxwell | `RD_dpm_lay1` | Layer 1: distance along the pole-side slot wall from the slot's inner corner (on the rib line) to the magnet's inner end |
| `pm2_end_gap` | 2.0 | mm | maxwell | `RD_dpm_lay2` | As pm1_end_gap for layer 2 |
| `pm1_length_ratio` | 0.899999999211849 | - | maxwell | `RD_K_Lpm_lay1` | Layer 1: magnet length / free slot length (pole-side wall from the magnet's inner end to the bridge circle); the magnet length is derived from it |
| `pm2_length_ratio` | 0.9 | - | maxwell | `RD_K_Lpm_lay2` | As pm1_length_ratio for layer 2 |
| `pm1_rib_half_width` | 1.0 | mm | maxwell | `RD_FB_lay1` | Layer 1: HALF the width of the centre rib on the d-axis; the slot's inner end is the straight line y = this value (+ pm1_inner_offset) |
| `pm2_web_width` | 2.31010869033253 | mm | maxwell | `RD_FB_lay2` | Layer 2: thickness of each of the two steel webs between the tilted inner end of the slot and the central cavity, normal to the tilted end; the slot's inner pole-side corner sits this value + pm2_inner_offset above the d-axis |
| `pm1_inner_offset` | 0.0 | mm | maxwell | `RD_INNO_lay1` | Layer 1: extra height of the slot's inner end above the rib |
| `pm2_inner_offset` | 1.05729199744521 | mm | maxwell | `RD_INNO_lay2` | Layer 2: extra height of the slot's inner end above the web |
| `pm2_inner_tilt` | 9.14233419891509 | deg | maxwell | `RD_thetaINNO_lay2` | Layer 2: tilt of the slot's inner end against the d-axis |
| `pm1_bridge_thickness` | 1.08008400158923 | mm | maxwell | `RD_WO_lay1` | Layer 1: radial thickness of the outer bridge, measured inwards from the lamination OD; the slot is closed by the circle r_lam - this value (in the bridgeless rotor it only sets the free slot length, i.e. the magnet length) |
| `pm2_bridge_thickness` | 1.0 | mm | maxwell | `RD_WO_lay2` | As pm1_bridge_thickness for layer 2 |
| `pm1_slot_start` | 1.0 | mm | maxwell | `Rectangle1 XStart` | Sleeved rotors: distance along the magnet axis from the d-axis crossing to the square inner end of the layer-1 slot (strip) |
| `pm2_slot_start` | 5.0 | mm | maxwell | `Rectangle11 + MoveEdges` | As pm1_slot_start for layer 2 (model: strip from 0 mm, inner edge moved 5 mm) |
| `pm_side_clearance` | 0.1 | mm | maxwell | `MoveEdges -0.1mm` | Magnet thinner than the slot: its q-side face is moved this far from the q-side slot wall (the pole-side face stays on the wall) |
| `pm_end_clearance` | 0.1 | mm | maxwell | `MoveEdges -0.1mm` | Magnet's inner end face moved this far towards the outer end (bridgeless model: 0) |
| `pm_corner_radius` | 1.0 | mm | maxwell | `Fillet 1mm` | Radius of the four magnet corners (0 = sharp) |
| `sleeve_thickness` | 1.5 | mm | paper | `Wcf` | Carbon-fibre sleeve radial thickness (> 0); the lamination OD is stator_id/2 - airgap - this value. Paper compares 1.5 mm (bridged) and 2.0 mm (bridgeless) |
| `shaft_diameter` | 60.0 | mm | maxwell |  | Lamination bore diameter |
| `q_hole_inner_radius` | 42.0 | mm | maxwell |  | Lightening hole on each q-axis: radius of its inner arc (centred on the shaft axis) |
| `q_hole_apex_radius` | 56.0 | mm | maxwell |  | Lightening hole: radius of its outer apex on the q-axis |
| `q_hole_half_angle` | 20.0 | deg | maxwell |  | Lightening hole: half the angular span of its inner arc, about the q-axis |
| `q_hole_fillet` | 3.0 | mm | maxwell |  | Lightening hole: radius of its three corner fillets |
| `d_hole_diameter` | 4.0 | mm | maxwell |  | Round hole on the d-axis near the OD |
| `d_hole_depth` | 4.1 | mm | maxwell |  | Lamination OD radius minus the d-axis hole centre radius (model: centre at 68.5 mm - Wcf) |
| `lip_depth` | 1.0 | mm | maxwell | `w_tip` | Bridgeless: radial depth of the pole-tip lips that cover the slot mouths under the OD |
| `pm1_opening_angle` | 1.0 | deg | maxwell | `rib2_open` | Bridgeless: angular width of the OD opening of the layer-1 slot |
| `pm2_opening_angle` | 1.0 | deg | maxwell | `rib1_open` | Bridgeless: angular width of the OD opening of the layer-2 slot |
| `pm1_opening_shift` | 0.0 | deg | illustrative | `1.5deg in Angle_rib3` | Bridgeless: the q-side edge of the layer-1 opening lies this angle towards the d-axis from the opening reference angle (where the magnet centre line meets the OD). Model: 1.5 deg, compensating its layer-1 reference typo; 0 with the true reference (docs/methods.md, 1.4) |
| `pm2_opening_shift` | 0.0 | deg | maxwell |  | As pm1_opening_shift for layer 2 |
| `min_web` | 0.5 | mm | illustrative |  | Minimum steel between pockets/holes, and between pockets and the bore (feasibility check) |
| `min_q_web` | 1.0 | mm | illustrative |  | Minimum steel width on the q-axis between adjacent poles (feasibility check) |
| `lam_bh` | generic | - | illustrative |  | Lamination BH curve: 'generic' (labelled generic NO-steel curve), a path to a two-column CSV (H [A/m], B [T]) or 'aedt:<library material name>'. The paper's 27SW1400 data are not shipped |
| `pm_br` | 1.4 | T | illustrative |  | PM remanence (linear magnet). Paper: G52UH, 'higher than 1.40 T' |
| `pm_mur` | 1.05 | - | illustrative |  | PM recoil permeability |
| `cu_conductivity` | 58000000.0 | S/m | illustrative |  | Copper conductivity (AC-loss runs) |
| `sleeve_sigma_fibre` | 0.0 | S/m | illustrative |  | Sleeve ply electrical conductivity ALONG the fibre (sigma_L). With sleeve_sigma_transverse it gives the axial conductivity of the sleeve (sleeve_conductivity_axial = sigma_L sin^2 t + sigma_T cos^2 t per layer at winding angle t from hoop, layers thickness-weighted), the only one a 2D model with axial current uses. 0 and 0 (default) = non-conducting sleeve, model and results unchanged; > 0 axial conductivity adds the 2D sleeve eddy-current loss to EM runs (docs/methods.md, 2.1). Illustrative class statement only, not a material property: carbon/epoxy is often quoted at 1e3-1e5 S/m along the fibre and orders of magnitude lower across it; measure your laminate |
| `sleeve_sigma_transverse` | 0.0 | S/m | illustrative |  | Sleeve ply electrical conductivity ACROSS the fibre (sigma_T), set by fibre-to-fibre contacts through the matrix; for a hoop-wound sleeve (winding angle 0) it is the axial conductivity. Illustrative class statement only: typically several orders of magnitude below sigma_L and strongly process-dependent; measure your laminate |
| `current_rms` | 420.0 | A | paper |  | Phase current (rms); paper's full-load value |
| `gamma` | 130.0 | deg(e) | illustrative |  | Current advance angle measured from the +d axis (90 = pure q-axis current) |
| `em_rpm` | 3000.0 | rpm | illustrative |  | Speed for torque runs (no eddy regions) |
| `em_span` | 60.0 | deg(e) | illustrative |  | Simulated electrical angle per torque run |
| `em_step` | 2.0 | deg(e) | illustrative |  | Time step expressed in electrical degrees |
| `mesh_rotor` | 1.2 | mm | illustrative |  | Max element length in rotor iron and magnets |
| `mesh_stator` | 1.5 | mm | illustrative |  | Max element length in stator iron and conductors |
| `mesh_gap` | 0.3 | mm | illustrative |  | Max element length in the air-gap layers |
| `lam_density` | 7850.0 | kg/m3 | paper |  | Lamination density |
| `lam_E` | 200.0 | GPa | paper |  | Lamination Young's modulus |
| `lam_nu` | 0.3 | - | paper |  | Lamination Poisson ratio |
| `pm_density` | 7400.0 | kg/m3 | paper |  | PM density |
| `pm_E` | 160.0 | GPa | paper |  | PM Young's modulus |
| `pm_nu` | 0.24 | - | paper |  | PM Poisson ratio |
| `sleeve_density` | 1800.0 | kg/m3 | paper |  | Sleeve density |
| `sleeve_E_hoop` | 125.0 | GPa | paper |  | Sleeve modulus along the fibre (hoop). Table II lists 125 GPa as 'radial'; see docs/methods.md |
| `sleeve_E_radial` | 8.8 | GPa | paper |  | Sleeve transverse (radial) modulus |
| `sleeve_E_axial` | 8.8 | GPa | illustrative |  | Sleeve axial modulus (transversely isotropic) |
| `sleeve_nu_hoop_radial` | 0.28 | - | paper |  | Major Poisson ratio -eps_r/eps_hoop under hoop stress. The minor ratio follows from symmetry (nu_hr*E_r/E_h) |
| `sleeve_nu_hoop_axial` | 0.28 | - | illustrative |  | Major Poisson ratio hoop-axial |
| `sleeve_nu_radial_axial` | 0.4 | - | illustrative |  | Transverse Poisson ratio radial-axial |
| `sleeve_G_hoop_radial` | 4.5 | GPa | illustrative |  | In-plane shear modulus |
| `sleeve_G_hoop_axial` | 4.5 | GPa | illustrative |  | Shear modulus hoop-axial |
| `sleeve_G_radial_axial` | 3.1 | GPa | illustrative |  | Transverse shear modulus |
| `sleeve_fibre` |  | - | illustrative |  | Sleeve fibre class from sipm.composite.FIBRES (carbon_sm, carbon_im, carbon_hm, e_glass), or '' (default) to use the sleeve constants as entered. When set, the thirteen sleeve constants (sleeve_E_*, sleeve_nu_*, sleeve_G_*, sleeve_density, sleeve_cte_*) are DERIVED from fibre, matrix, vf, micromechanics model and winding angle whenever the parameters are created or replaced, overwriting any values passed for them |
| `sleeve_matrix` | epoxy | - | illustrative |  | Sleeve matrix class from sipm.composite.MATRICES (epoxy, peek); used only with sleeve_fibre |
| `sleeve_vf` | 0.55 | - | illustrative |  | Sleeve fibre volume fraction, in (0, 1); used only with sleeve_fibre |
| `sleeve_micromechanics` | chamis | - | illustrative |  | Micromechanics model for the transverse and shear ply constants: chamis, halpin_tsai or rom; used only with sleeve_fibre |
| `sleeve_wind_angle` | 0.0 | deg | illustrative |  | Winding angle of a balanced +/- angle-ply sleeve, measured from the hoop direction (0 = hoop winding), in [0, 90]. Non-zero needs sleeve_fibre (the ply cannot be recovered from homogenised constants) |
| `sleeve_layers` |  | - | illustrative |  | Layered sleeve: comma-separated 'angle:weight' items from the bore outward, e.g. '0:2,20:1,0:1' (angle = balanced +/- winding angle from hoop in [0, 90] deg, weights normalised to fractions of sleeve_thickness); '' (default) = one homogeneous layer at sleeve_wind_angle. Needs sleeve_fibre (one fibre, matrix and vf for all layers) and sleeve_wind_angle = 0; at most 10 layers, each at least 0.05 mm thick. In the deck each layer is its own area and material, perfectly bonded (shared nodes). The thirteen sleeve constants then hold the membrane homogenisation of the stack (layers share the hoop and axial strain and the radial stress; sipm.composite.layer_stack), i.e. the stiffness and free thermal strain of an equivalent thin homogeneous sleeve: for inspection and the cache key only, the deck and the failure evaluation use each layer's own constants (docs/methods.md, 4.6) |
| `sleeve_Xt` | 2000.0 | MPa | illustrative |  | Sleeve ply strength, fibre-direction tension (default: STRENGTHS['carbon_sm_epoxy'], a class-typical round number) |
| `sleeve_Xc` | 1200.0 | MPa | illustrative |  | Sleeve ply strength, fibre-direction compression |
| `sleeve_Yt` | 50.0 | MPa | illustrative |  | Sleeve ply strength, transverse tension |
| `sleeve_Yc` | 200.0 | MPa | illustrative |  | Sleeve ply strength, transverse compression |
| `sleeve_S` | 70.0 | MPa | illustrative |  | Sleeve ply shear strength (also used for the transverse shear strength in Hashin's matrix-compression mode) |
| `lam_cte` | 1.2e-05 | 1/K | illustrative |  | Lamination thermal expansion coefficient (electrical steel, typically 11-13e-6 /K) |
| `pm_cte` | 2e-06 | 1/K | illustrative |  | PM thermal expansion coefficient, one isotropic value. Sintered NdFeB is strongly anisotropic (about +5e-6 /K along the magnetisation, about -1e-6 /K across it); the default is the mean of the two in-plane directions |
| `sleeve_cte_hoop` | 0.0 | 1/K | illustrative |  | Sleeve thermal expansion along the fibre (hoop); carbon/epoxy is about zero (typically -1e-6 .. +0.5e-6 /K) |
| `sleeve_cte_radial` | 3e-05 | 1/K | illustrative |  | Sleeve transverse (radial) thermal expansion, matrix-dominated (typically 25-35e-6 /K for carbon/epoxy) |
| `sleeve_cte_axial` | 3e-05 | 1/K | illustrative |  | Sleeve axial thermal expansion (transversely isotropic: = radial); enters plane strain only |
| `mech_rpm` | 20000.0 | rpm | paper |  | Speed for the stress analysis |
| `interference` | 0.1 | mm | illustrative |  | Radial interference between lamination OD and sleeve bore |
| `sleeve_layer_prestress` |  | MPa | illustrative |  | Winding prestress of the sleeve layers: one value per layer (sleeve_layers order, from the bore) or one value for all, each >= 0 ('' = none). Value sigma0 is a hoop eigenstrain: the layer's stress-free hoop strain is -sigma0/E_hoop of that layer (sigma0 = the hoop stress the layer carries when held at its fitted radius with no radial stress), applied to all layers at once before elastic equilibrium (MAPDL INISTATE initial elastic strain); the outputs are the relaxed state. Acts together with interference (which may be 0), temperature and speed. A simultaneous-eigenstrain model, not the layer-by-layer build-up of real winding; `sipm winding` converts a winding tension into equivalent values (docs/methods.md, 4.6) |
| `mech_temperature` | 20.0 | degC | illustrative |  | Uniform rotor temperature (lamination, magnets and sleeve) of the stress analysis |
| `mech_ref_temperature` | 20.0 | degC | illustrative |  | Reference temperature: the temperature at which the parts have their stress-free fitted dimensions, i.e. at which interference is defined. Equal to mech_temperature (default) = isothermal analysis |
| `mu_sleeve` | 0.0 | - | illustrative |  | Friction coefficient sleeve/lamination |
| `mu_magnet` | 0.0 | - | illustrative |  | Friction coefficient magnet/pocket |
| `plane` | stress | - | illustrative |  | 2D assumption: 'stress' or 'strain' |
| `mech_mesh` | 0.5 | mm | illustrative |  | Global element size (mechanical) |
| `mech_mesh_fine` | 0.15 | mm | illustrative |  | Element size on short lines (fillets, bridges, lugs) |
| `magnet_contact` | bonded | - | illustrative |  | Magnet/pocket contact in the stress model: 'bonded' (as the original Workbench model; the designed clearance is bridged by the bond, see docs/methods.md) or 'frictional' (open clearance, mu_magnet; not supported for sleeved_bridgeless) |
| `stress_corner_margin` | 0.5 | mm | illustrative |  | Named stress sections (svm_lam_max) leave out this length at both ends of each section, where it meets a sharp corner of the outline |
| `mech_cases` |  | - | illustrative |  | Comma-separated mechanical load cases solved for each point, from standstill_cold, standstill_hot, rated_cold, rated_hot, overspeed_cold, overspeed_hot ('' = the single case mech_rpm at mech_temperature, as before). Speeds: standstill 0, rated mech_rpm, overspeed overspeed_factor x mech_rpm; temperatures: hot mech_temperature, cold mech_temperature_cold. Each output gets one column per case (e.g. sh_sleeve_max__overspeed_hot), plus worst-case envelope columns and the lift-off speeds (docs/methods.md, 3.1) |
| `overspeed_factor` | 1.2 | - | illustrative |  | Overspeed load cases: speed = this factor x mech_rpm (>= 1) |
| `mech_temperature_cold` | 20.0 | degC | illustrative |  | Uniform rotor temperature of the '_cold' load cases (the '_hot' cases use mech_temperature); interference is defined at mech_ref_temperature |

Topology-specific overrides of the defaults above:

* **sleeved_bridged**: `sleeve_thickness`=1.5 (saved state of 01_ICEM_01_simplerotor except the sleeve: 1.5 mm as in the paper's Section III (model saved at 2.0 mm))
* **sleeved_bridgeless**: `sleeve_thickness`=0.5, `pm1_rib_half_width`=0.1, `pm2_web_width`=0.1, `pm2_inner_offset`=0.0, `pm1_bridge_thickness`=0.5, `pm2_bridge_thickness`=0.5, `pm_end_clearance`=0.0 (saved state of 02_ICEM_01_teslamotor1 (sleeve 0.5 mm; the paper's bridgeless comparison uses 2.0 mm))

## Fixed shape constants

Polygons of the Maxwell models that are design features but not parameters (mm). Each lug or
relief polygon is given in a frame whose origin is a corner of the nominal magnet and whose x
axis points to a second corner (y = x turned +90 deg); `W/2` = half the slot width. Lugs are
added to the lamination, reliefs removed from it.

| Feature | Layer | Origin -> x axis | Polygon |
|---|---|---|---|
| `layer2_inner_lug` | 2 | P9 -> P10 | (0.0, -0.6), (0.0, -W/2), (-2.0, 0.8), (-0.8, 0.4) |
| `layer2_outer_lug` | 2 | P10 -> P11 | (0.4, 0.0), (1.6, 0.0), (-0.6, 2.2), (-0.6, 0.8) |
| `layer1_inner_lug` | 1 | P9 -> P8 | (0.6, 0.0), (W/2, 0.0), (-0.4, -1.8), (-0.6, -0.8) |
| `layer1_outer_lug` | 1 | P10 -> P11 | (0.5, 0.0), (1.5, 0.0), (-1.0, 2.0), (-1.0, 1.0) |
| `layer2_inner_relief` | 2 | P8 -> P11 | (1.0, 0.0), (0.0, -0.2), (-0.4, 0.0) |
| `layer2_outer_relief` | 2 | P11 -> P10 | (0.0, 0.8), (-0.2, 0.0), (0.0, -0.2) |
| `layer1_inner_relief` | 1 | P8 -> P11 | (0.6, 0.0), (0.0, -0.1), (-0.1, 0.0) |
| `layer1_outer_relief` | 1 | P11 -> P8 | (-0.2, 0.0), (0.0, 0.2), (0.8, 0.0) |

Central cavity (`CAVITY`, layer-2 points FB and FBd, mm): `tip` = 1.0, `step` = 2.0, `tangent_tip` = 3.0, `tangent_tip_steep` = 6.0, `tangent_fb` = 2.0, `tangent_fbd` = 6.0.
