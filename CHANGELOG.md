# Changelog

## Unreleased

### Review fixes

* **Thin sleeves: contact target meshed like the lamination.** With a single-layer sleeve thinner
  than 3 x `mech_mesh` (the bridgeless rotor's default 0.5 mm sleeve, a 1.0 mm bridge-retained
  sleeve, or the 1.5 mm default at a 0.8 mm mesh) the sleeve area size (thickness / 3) made the
  bore arc, the target of the sleeve contact, finer than the lamination OD; it is now sized
  `mech_mesh`, as for a layered sleeve, on both rotors (only when the area size is finer). The
  bridgeless default deck changes by that one line (pinned hash updated); headline, global and
  magnet peaks move by < 0.01 MPa, sleeve hoop 236.6 to 236.5 MPa, lips by up to 0.4 MPa,
  `tau_sleeve_max` 6.0 to 3.8 MPa (docs/validation.md sections 2, 3, 11). The default
  bridge-retained deck (1.5 mm sleeve, 0.5 mm mesh: equal sizes) is byte-identical (pinned hash
  unchanged). Re-solved: the 1.0 mm rows of `examples/data/sample_mech_results.csv` and of
  docs/validation.md section 7 (section peaks within 0.03 MPa, sleeve hoop within 0.2 MPa, most
  compressive sleeve radial stress up to 12 %, peak contact pressure up to 3.4 %; the 1.5 and
  2.0 mm rows re-solved bit-identical), `docs/img/stress_map.png`, and the 0.8 mm bridge-retained
  row of the mesh study (section 3: `bridge2` 97.4 to 97.3, sleeve hoop 222.9 to 222.6 MPa).
* **Cached points changed by a solver rule are re-solved.** The point key cannot see a change of
  the deck rules, so cached thin-sleeve mech points (and bridgeless load-case points, now with a
  lift-off value) would have been reused with their earlier results. `sipm.sweep.point_rules`
  names the rules that apply to a point (`bore-arc-mesh`, `liftoff-search`); the cache record of
  such a point lists them, and a cached record without them is solved again under the same key
  (`results()` shows it as `pending` until then). Keys are unchanged for every point (pinned
  default keys and the keys of `examples/data/sample_mech_results.csv`), and points no rule
  applies to (every EM point, the default bridge-retained mech point, sleeves at least
  3 x `mech_mesh` thick, layered sleeves) keep their records and reuse.
* **Bridgeless lift-off speed by search.** `liftoff_rpm_cold` / `_hot` of `sleeved_bridgeless`
  were NaN (contact not linear in speed). They are now found directly with the definition of the
  line (lowest speed at which an element's contact pressure reaches zero; open at standstill = 0):
  bracket from the solved cases, extra single-case MAPDL runs at trial speeds placed from the
  local extrapolation through the closed solutions (bisection as fallback) until the bracket is at
  most `liftoff_search_tol` wide (new parameter, 100 rpm, out of the point key and the EM key at
  its default; 0 = no search, NaN as before), at most 12 extra runs per temperature; the basis
  says "found by iteration" with the bracket, the cases and the number of extra runs. Default
  bridgeless rotor: 24 414 rpm cold (2 extra runs; the old line gave 23 591 rpm, where the
  pressure was still positive), 0 hot (open at standstill). Forced on the bridge-retained rotor
  it gives 39 512 rpm against the line's 39 357 rpm (+0.39 %, within the line's 0.5 % accuracy;
  docs/validation.md section 11). The bridge-retained rotor keeps the two-point extrapolation
  (same runs, same values, same decks); `solve_mech_cases(..., liftoff_search=True/False)`
  forces the search on or off.
* **Open-circuit and no-load field runs model the sleeve non-conducting** (as `ac_loss`): they
  start from a static field, so a conducting sleeve added its eddy start-up transient. With
  1e4 S/m axial the bridge-retained open-circuit EMF is the non-conducting 99.240 / 99.239 /
  99.244 V. Offline tests with a recording session for these runs and for `ac_loss`.
* README "Proposed next": the entries on 3D sleeve eddy loss, axial effects, more detailed
  temperature, 3D failure criteria and calibration were removed (out of scope); the limits in
  section 6 stay.
* **`sipm validate sleeve-loss`** (`sipm.validation.sleeve_loss_check`): the FE sleeve loss against
  the thin-shell estimate of docs/validation.md section 13 from the package (re-run: 0.9929).
* `MaxwellSession.ac_loss` models the sleeve non-conducting, so a conducting sleeve's loss is
  not counted as copper loss in `kac`.
* `failure_indices` rejects `f12_star` outside (-1, 1) (open Tsai-Wu envelope).
* Citations checked against publisher / DOI / NTRS records: Chamis cited as NASA TM-83320
  (1983); issue numbers added; unconfirmed chapter and equation numbers removed.
* ansys-marked tests for the bridgeless rotor with a layered, prestressed sleeve and with a
  conducting sleeve, and for `sleeve_loss_check`.

### Sleeve eddy-current loss (2D, axial current)

* **Sleeve conductivity.** New parameters `sleeve_sigma_fibre` and `sleeve_sigma_transverse`
  (S/m, ply values along and across the fibre; default 0 and 0 = non-conducting sleeve, model,
  cache keys and results unchanged). Derived `Params.sleeve_conductivity_axial`: per layer
  sigma_L sin^2 t + sigma_T cos^2 t (t = winding angle from hoop; a hoop winding gives the
  across-fibre value), layers thickness-weighted (in parallel); helpers
  `sipm.composite.ply_axial_conductivity` and `axial_conductivity`. No conductivity is shipped
  as a material property.
* **Cache keys.** The two fields are left out of the point key at their default
  (`KEY_OMITTED_AT_DEFAULT`) and always out of `em_model_key`, which adds the derived axial
  conductivity when it is > 0: every non-conducting sleeve keeps its key (pinned against the
  previous code), winding angle and layers reach the EM key only through the axial conductivity.
* **EM runs with a conducting sleeve** (`physics = "em"`, axial conductivity > 0): sleeve
  conductivity and eddy effects on (set directly: PyAEDT's `eddy_effects_on` skips materials
  below 1e5 S/m), no net-current constraint (the odd-periodic one-pole model makes the ring total
  zero; Maxwell leaves a passive solid's net current free), same mesh as before; the run covers
  a lead-in plus two 60-degree loss periods (180 electrical degrees by default). New columns
  `sleeve_loss_w` (time average of `SolidLoss` over the last loss period, whole machine, stack
  length), `sleeve_loss_rel_change` (change from the period before), `sleeve_loss_frac` (over
  shaft power), `sleeve_sigma_axial`, `sleeve_skin_depth_mm` and `sleeve_skin_ratio` (at the
  slot-passing frequency), `sleeve_mesh_mm`, `sleeve_steps_per_slot`; torque is averaged over the
  last `em_span` degrees. A run with too few time steps per slot period or a mesh coarser than
  half the skin depth is refused. Non-conducting runs and their columns are unchanged.
* **`sipm.eddy`**: loss period, skin depth, `period_average`, and the thin-shell estimate of the
  resistance-limited loss from the non-conducting field (`shell_loss_from_potential`, harmonic
  form `shell_loss_from_br`, with the sleeve's J x B torque, `rotor_frame`).
  `MaxwellSession.sleeve_lines` exports field quantities on arcs through the sleeve over a run,
  `parse_line_report_csv` reads them and `sipm.validation.sleeve_shell_estimate` turns a
  non-conducting export into the estimate.
* **CLI / config**: `sipm composite --sigma-fibre S --sigma-transverse S` prints the axial
  conductivity (per layer with `--layers`); `examples/configs/em_sleeve_loss.toml`; plot labels.
* **Validation** (docs/validation.md section 13) and the limits of a 2D axial-current loss
  (docs/methods.md 2.1): for hoop windings it is the across-fibre, axial-current part of the
  loss only; for helical layers it over-counts by leaving out the end-return resistance.

### Layered sleeve, winding prestress, disc + N rings closed form

* **Layered sleeve.** New parameter `sleeve_layers` (default '' = the single-layer sleeve, deck
  and results unchanged): `angle:weight` items from the bore outward, e.g. `"0:2,20:1,0:1"`,
  weights normalised to fractions of `sleeve_thickness`; needs a composite spec
  (`sleeve_fibre`, one ply for all layers) and `sleeve_wind_angle = 0`; at most 10 layers of at
  least 0.05 mm. In the deck every layer is its own area and material on shared lines (perfect
  bond), with at least two quadratic elements through each layer and the contact target meshed
  like the lamination. The thirteen sleeve constants then hold the membrane homogenisation of
  the stack (`sipm.composite.layer_stack`, `LayerStack`), for inspection and the cache key.
* **Per-layer post-processing.** Stresses are sorted and exported layer by layer, because
  MAPDL averages nodal results across material boundaries (an interface node showed the mean of
  the two layers' hoop stresses, 311 MPa between 445 and 177 MPa); the sleeve node file gets a
  seventh column `layer` (`parse_sleeve_nodes` returns it). Every node is evaluated for failure
  with the laminate of its layer. The sleeve columns keep their meaning as extremes over the
  whole sleeve; new columns for layered sleeves only: `sh_sleeve_max_L<i>`, `sleeve_sf_L<i>`,
  `sleeve_sf_layer`, `tau_sleeve_interface_max` (interlaminar r-theta shear; also as a
  worst-case envelope over load cases). Plot labels for them.
* **Winding prestress.** New parameter `sleeve_layer_prestress` (MPa, one value per layer or one
  for all, default ''): a hoop eigenstrain eps0 = sigma0 / E_hoop of each layer applied to all
  layers at once (MAPDL `INISTATE` initial elastic strain in the element axes), together with
  interference (which may be 0), speed, temperature and every load case of `mech_cases`; the
  outputs are the relaxed state. A simultaneous-eigenstrain model, not the layer-by-layer
  build-up. `sipm winding --tension T` (`sipm.analytic.sequential_winding`,
  `sipm.validation.winding_prestress`) converts a winding stress laid layer by layer onto the
  closed-form model into equivalent per-layer values that reproduce every layer's hoop force and
  the interface pressures of the build-up.
* **Closed form for N rings.** `sipm.analytic.disc_rings` (`RingStack`), `rings_summary`,
  `rings_coupled_pressure`, `rings_liftoff_speed`, `free_strain_terms`, `hoop_modulus`: disc +
  N bonded orthotropic rings with per-ring thermal strain and hoop eigenstrain; N = 1 is
  bit-identical to `disc_ring`, whose API and numbers are unchanged.
* **Validation** (docs/validation.md section 12): 3-layer ring (hoop / +/-30 deg / hoop,
  interference, prestress, 20000 rpm, +100 K) against the closed form within 0.15 % for every
  layer, the interfaces and the contact pressure, lift-off speed within 3e-6, `sipm validate
  ring-multilayer`; single-layer prestress on the ring within 0.12 %; rotor study of a hoop,
  a hoop/helical/hoop and a prestressed hoop sleeve. `sipm composite --layers` prints every
  layer's constants; `examples/configs/mech_multilayer_sleeve.toml`. `sipm deck/params --set`
  keep string fields strings. New point keys only when a new parameter leaves its default.

### Load cases as outputs and the lift-off speed

* **Load-case set.** New parameters `mech_cases` (default '' = the single case `mech_rpm` at
  `mech_temperature`, deck and columns unchanged), `overspeed_factor` (1.2) and
  `mech_temperature_cold` (20 degC). `mech_cases` lists cases from `standstill_cold`,
  `standstill_hot`, `rated_cold`, `rated_hot`, `overspeed_cold`, `overspeed_hot`; each is its
  own MAPDL run of the single-case deck (`Params.for_case`, `sipm.mapdl.solve_mech_cases`), so
  a case's result does not depend on the other cases or their order, also with friction.
  Columns: `<output>__<case>` for every output and case; the un-suffixed columns hold the case
  named in `mech_case` (`rated_hot` if requested, else the first listed); worst-case envelope
  `<output>__worst` with the governing case in `<output>__worst_case` for `svm_lam_max`,
  `sh_sleeve_max`, `sr_sleeve_max`, the sleeve failure indices, `sleeve_sf`,
  `sleeve_open_fraction` and `contact_p_min` (`sipm.loadcases`). TOML `[base]` may give the
  cases as an array; `examples/configs/mech_load_cases.toml`; plot labels for the new columns
  and stress-line plots of the worst case in `sipm plot`.
* **Lift-off speed.** `liftoff_rpm_cold` / `liftoff_rpm_hot` (with `liftoff_basis_*`): the
  speed at which the first sleeve contact element's pressure reaches zero, extrapolated element
  by element (p = p0 - c omega^2, exact while the contact is closed and frictionless and the
  magnets bonded) from two closed cases at the same temperature; 0 if open at standstill, NaN
  without two closed cases, with friction or for `sleeved_bridgeless` (its contact is not
  linear in speed), `inf` if the pressure does not fall with speed. With three or more closed
  speeds the middle ones are checked against the line and the basis is flagged above 1 % of
  the standstill pressure.
  Decks of a load-case set export the per-element contact pressure
  (`sipm_contact_elems.txt`, `parse_contact_elements`); single-case decks are byte-identical
  to before. Closed form `sipm.analytic.liftoff_speed` (with temperature) and
  `coupled_pressure`; `sipm validate ring-liftoff` / `liftoff`. Validation (docs/validation.md
  section 11): FE against closed form on the ring 3e-6; on the bridge-retained rotor a direct run
  at the extrapolated speed opens the first element (about 0.5 %); on the bridgeless rotor the
  linearity fails (6.6 % deviation, edges open at standstill when hot).
* `write_deck` gains the keyword `contact_elements` (default False). New point keys only when
  a new parameter leaves its default, so existing caches stay valid.

### Composite sleeve: thermal load case, sleeve stress field, ply failure

* **Thermal load case.** New parameters `lam_cte`, `pm_cte`, `sleeve_cte_hoop`,
  `sleeve_cte_radial`, `sleeve_cte_axial` (1/K), `mech_temperature` and `mech_ref_temperature`
  (degC; the interference is defined at the reference temperature). The deck writes `TREF`,
  `BFUNIF,TEMP` and `MP,ALPX/ALPY/ALPZ` (the sleeve's in its element axes) only when the two
  temperatures differ; the defaults are equal, so every existing result is unchanged
  (re-run, bit-identical). `sipm.analytic` gains the thermal particular term for the
  orthotropic ring (valid through k = 1) and `ring_check("thermal")` validates the FE path
  against it (agreement 0.12 %, docs/validation.md section 1).
* **Sleeve stress field.** The deck writes the sleeve corner nodes (r, theta and the stresses
  in the cylindrical system) to `sipm_sleeve_nodes.txt`; `sipm.mapdl.parse_sleeve_nodes` reads
  it and `solve_mech` returns the arrays as `r["_sleeve_nodes"]` (not cached). New result
  columns `sr_sleeve_max` (max radial stress: radial tension, the delamination driver) and
  `tau_sleeve_max` (max |r-theta shear|).
* **Composite sleeve specification** (`sipm.composite`): micromechanics (Chamis, Halpin-Tsai,
  rule of mixtures; Schapery CTEs), balanced +/- winding homogenisation and ply failure
  criteria. New parameters `sleeve_fibre` (default '' = constants as entered, as before),
  `sleeve_matrix`, `sleeve_vf`, `sleeve_micromechanics`, `sleeve_wind_angle`: with a fibre set,
  the thirteen sleeve constants are derived from the spec and overwrite passed values. Presets
  are illustrative class-typical numbers, not commercial grades. `sipm composite --fibre ...
  [--toml]` prints the derived constants without Ansys; `examples/06_composite_sleeve.py`
  tabulates them, `examples/configs/mech_composite_sleeve.toml` sweeps winding angle and
  temperature.
* **Sleeve failure outputs.** Ply strengths `sleeve_Xt`, `sleeve_Xc`, `sleeve_Yt`, `sleeve_Yc`,
  `sleeve_S` (illustrative carbon_sm/epoxy defaults). New columns `sleeve_fi_max_stress`,
  `sleeve_fi_tsai_wu`, `sleeve_fi_hashin_fibre`, `sleeve_fi_hashin_matrix` and the Tsai-Wu
  safety factor `sleeve_sf` with its location `sleeve_sf_r`, `sleeve_sf_theta`
  (`sleeve_fi_tsai_wu` = 1 / `sleeve_sf`, linear in load and never negative), evaluated in
  the +/- plies (or one hoop ply for directly entered constants) on both transverse planes.
* **Sleeve constant validation** now requires a positive definite 6x6 sleeve compliance (what
  MAPDL requires) instead of bounding each Poisson ratio to [0, 0.5): angle-ply constants
  with nu > 1 or nu < 0 are accepted. Sets accepted before are still accepted unless their
  compliance was indefinite.
* **Sweep caches stay valid.** The new parameters are left out of the point key while they sit
  at their defaults, so existing EM and mech caches are reused. Old cached mech points lack
  the new columns (`sr_sleeve_max`, `tau_sleeve_max`, the sleeve failure columns); use a new
  sweep name to get them. None of the new parameters enters the EM model key.
* The Table II default constants are unchanged, including the published density of
  1800 kg/m3 (the composite spec gives 1522 kg/m3 for carbon_sm/epoxy at vf = 0.536, where E1
  matches Table II's 125 GPa; docs/methods.md section 4.5).

### README and small fixes

* README reorganised as a reader's guide: what it answers, a five-minute try-out without
  Ansys, running your own design, running the solvers, reading the results, assumptions and
  limits, citation.
* A shaft bore that comes within `min_web` of the lightening holes' inner arc is now rejected
  with a clear GeometryError (it used to fail inside the outline evaluation, or build a rotor
  whose holes open into the bore).
* Geometry figure: callout leaders no longer cross or stop short of their feature.

### Conventional rotor removed, bonded magnets, stress at named sections (breaking)

* **The conventional (unsleeved) rotor is removed** from the package and its geometry is no
  longer distributed: its barrier-pocket construction (the Hermite `RD_x13..x29` chains), its
  parameters and topology defaults, its reference fixture, its part of the fixture generator
  (`--highpower` input), its example config and its figures. `topology="conventional"` now
  raises a ValueError that says so. The stator reference check is kept with a stator-only
  fixture (`tests/fixtures/maxwell_ref/stator_saved.json`, the same outline). The only remnant
  is `examples/data/sample_em_conventional.csv`, reduced to operating-point and torque columns
  and labelled `conventional_baseline`: a torque curve computed earlier with that rotor, drawn
  as a reference baseline in `docs/img/torque_vs_gamma.png`. `docs/img/geometry_topologies.png`
  shows the two sleeved rotors.
* **Magnet contact option** `magnet_contact` = `"bonded"` (default, as the original Workbench
  model: bonded MPC contact that bridges the designed 0.1 mm clearance) or `"frictional"` (the
  previous open-clearance contact with stabilisation damping; refused for the bridgeless
  rotor, where it does not converge reliably). With bonded magnets the bridgeless stress model
  converges at its defaults and across its example sweep.
* **Headline stress at named sections**: `svm_lam_max` (= `svm_rotor_max`) is now the peak von
  Mises stress over named load-carrying sections (`sipm.geometry.locations`): bridge-retained
  `bridge1`, `bridge2`, `rib1`, `web2`, `qhole`; bridgeless `qhole`, `lip1_pole`, `lip1_q`,
  `lip2_pole`, `lip2_q`; each reported as `svm_<name>`, plus `svm_lam_where`. Each section
  leaves out `stress_corner_margin` (0.5 mm) at its sharp end corners. All bridge-retained
  sections and the bridgeless `qhole` converge within 2 % between the two finest meshes; the
  bridgeless lips do not (3-10 %, docs/validation.md section 3). The global lamination peak
  stays as `svm_lam_max_all` (mesh-dependent). Plots name the quantity ("Max von Mises at named
  lamination sections").
* Sample mechanical results re-solved with bonded magnets (`sample_mech_results.csv`, new
  `sample_mech_bridgeless.csv`); `stress_map.png` redrawn; validation.md mechanical sections
  rewritten.

What a user of the old names must change:

| old | new |
|---|---|
| `topology="conventional"`, `default_params("conventional")` | removed (ValueError); use a sleeved topology |
| `TOPOLOGIES` had three entries | two: `("sleeved_bridged", "sleeved_bridgeless")`; `REMOVED_TOPOLOGIES` explains the third |
| `Params.has_sleeve` | removed (every topology has a sleeve) |
| `pm1_barrier_factor`, `pm2_barrier_factor` (`if_RD_lay*`) | removed (conventional barrier curves only) |
| `pm1_inner_tilt` | removed (the sleeved layer-1 slot ends on the flat rib line; the tilt was required to be 0) |
| `contact_exclusion`, deck components `NEAR*`, results key `svm_lamx` | removed; replaced by the named sections and `stress_corner_margin` |
| `svm_lam_max` = lamination peak outside the magnet contact zones | = peak over named sections; per-section `svm_<name>`, `svm_lam_where` |
| magnets in open-clearance frictional contact (always) | `magnet_contact="bonded"` by default; `"frictional"` for the old model (bridge-retained only) |
| `MechModel.magnet_frames` | `MechModel.sections` |
| `LayerInput.if_rd`; construction points `P13`..`P29`; scalars `distance`, `distance4` | removed |
| `RotorPole.sleeve` could be `None` | always a `Loop` |
| `examples/configs/em_conventional.toml`, `docs/img/construction_conventional.png` | removed |
| `make_maxwell_fixtures.py --highpower` | removed; the stator fixture is no longer regenerated |
| fixture `conventional_saved.json` | `stator_saved.json` (stator and coils only) |

### Geometry rebased on the author's Maxwell models (breaking)

The rotor and stator are now constructed the way the author's original Maxwell 2D models
construct them (sleeved bridge-retained `01_ICEM_01_simplerotor`, sleeved bridgeless
`02_ICEM_01_teslamotor1`; a conventional rotor of this release was removed again, see above),
with the models' saved
states as defaults. Each magnet layer is a straight slot fixed by the d-axis crossing of its
magnet centre line, its angle and width; the magnet length is derived
(`length_ratio x (free slot length - end gap)`), so a thicker sleeve shortens the magnets
instead of moving them. See docs/methods.md, section 1.

* `sipm.geometry` is a package: `primitives` (exact line/arc/Hermite chains, fillets, biarcs,
  `Loop`), `csg` (exact Booleans), `construction` (RD_ points, derived magnet length), `rotor`
  (sleeved topologies), `checks` (feasibility checks, derived metrics), `stator`, `polygons`.
  `RotorPole` now carries `lamination` (pieces with holes, exact chains), `half_lamination()`
  and `layers` (construction points); `sector`/`cutters` of the old interface are replaced.
* New features: q-axis lightening holes, d-axis hole, central cavity, magnet lugs and corner
  reliefs, magnet clearance and corner radius, pole-tip lips with angular openings, three
  lamination pieces in the bridgeless rotor, Hermite cavity curves (sent to both solvers as
  biarcs within 5e-4 mm).
* New defaults: stator OD 210 mm (model; paper Table I: 220 mm), Maxwell slot (1.22 mm tip,
  70 deg taper, 5.27 mm wide, 15.75 mm deep, 0.3 mm fillets), one coil region per slot
  (8 x 1.57 x 4.3 mm, 7 turns), stack 130 mm, shaft 60 mm, bridgeless sleeve 0.5 mm.
* `maxwell.py`: several lamination pieces united into `Rotor`, one coil object per slot with
  `hairpin_layers` conductors; AC-loss (solid conductor) runs raise `NotImplementedError`
  until the strand geometry is defined. `mapdl.py`: lamination pieces each minus their own
  holes; magnet contact includes the initial gap (KEYOPT(12,9) = 0) when the magnets have a
  designed clearance.
* Reference fixtures `tests/fixtures/maxwell_ref` (2 saved states + 13 variable sets, and the
  stator) and
  their generator `tools/make_maxwell_fixtures.py` with the replay oracle `tools/maxwell_ref`;
  `tests/test_maxwell_reference.py` checks every vertex within 1e-6 mm.
* `docs/make_construction_figures.py` and `docs/img/construction_*.png`: construction diagrams.
* Bridgeless layer-1 OD openings deviate from the model on purpose: the model's opening
  reference `angle_temp_lay1` uses `RD_xpm_lay2` in its x term (a typo; 13.616 deg instead of
  the true 12.394 deg at the saved state), offset by a hand-set 1.5 deg shift. The package
  uses the true angle where the layer-1 magnet centre line meets the OD for both layers and
  `pm1_opening_shift` now defaults to 0 (was 1.5), so the layer-1 opening moves 0.28 deg away
  from the d-axis at the saved state (docs/methods.md, 1.4). The fixtures are unchanged; the
  reference test checks this deviation and the rib cut exactly.

| old | new | Maxwell variable | note |
|---|---|---|---|
| `pm1_length`, `pm2_length` | removed | `RD_Lpm_lay*` | derived: `rp.metrics["pm1_length"]` |
| `pm1_bridge_top`, `pm2_bridge_top` | `pm1_bridge_thickness`, `pm2_bridge_thickness` | `RD_WO_lay*` | same measure (radial, from the OD) |
| `pm1_bridge_bottom` | `pm1_rib_half_width` | `RD_FB_lay1` | HALF the rib: old 1.0 mm rib = new 0.5 |
| `pm2_bridge_bottom` | `pm2_web_width` | `RD_FB_lay2` | different meaning: web beside the central cavity |
| `pm1_inner_gap`, `pm2_inner_gap` | `pm1_end_gap`, `pm2_end_gap` | `RD_dpm_lay*` | different measure: along the pole-side wall from the slot corner P5 |
| `pm1_barrier`, `pm2_barrier` | `pm1_d_crossing`, `pm2_d_crossing` | `RD_xpm_lay*` | the layer position is now a parameter |
| `pm1_lug`, `pm2_lug`, `seat_radius` | `pm_side_clearance`, `pm_end_clearance`, `pm_corner_radius` | MoveEdges, Fillet | lugs/reliefs are fixed polygons (docs/parameters.md) |
| `fillet_radius` | removed | - | |
| `rotor_slot_opening` | `pm1_opening_angle`, `pm2_opening_angle`, `lip_depth` | `rib2_open`, `rib1_open`, `w_tip` | openings in deg |
| `slot_opening_depth`, `slot_wedge_depth` | `slot_tip_depth`, `slot_tip_angle` | `Wtip`, `tip_angle` | |
| `cond_gap`, `slot_liner` | `cond_start_radius`, `cond_count` | `Rcoil_start`, NumClones | |
| - | `pm1/2_length_ratio`, `pm1/2_inner_offset`, `pm2_inner_tilt`, `pm1/2_slot_start`, `pm1/2_opening_shift`, `q_hole_*`, `d_hole_*`, `slot_depth`, `slot_fillet` | `RD_K_Lpm`, `RD_INNO`, `RD_thetaINNO`, `if_RD`, ... | new |

`min_web` 1.0 -> 0.5 mm and `min_q_web` 2.0 -> 1.0 mm (illustrative thresholds; every corner of
the original sweep ranges builds).

### Results recomputed with the rebased geometry

(EM results stand; the mechanical results and limitations of this step are superseded by the
section above.)

* `examples/data/sample_em_results.csv`, `sample_em_conventional.csv` (now the full results
  table of a real run instead of values typed in from the validation notes), new
  `sample_em_bridgeless.csv` and `sample_mech_results.csv`: re-solved from `examples/configs`
  with Ansys Electronics Desktop 2024 R1 and Mechanical APDL 2024 R1. Outputs of this
  toolchain, not data of the companion paper.
* `docs/img/torque_vs_gamma.png` and `docs/img/stress_map.png` redrawn from them.
* `docs/validation.md` sections 1-8 recomputed: the Maxwell driver runs end to end on all three
  topologies, the MAPDL driver on the conventional and bridge-retained rotors, open-circuit phase/symmetry checks for every topology, an ampere-turn
  check of the 7-turn coil region, mesh studies, sweeps. Section 9 (AC loss) still describes
  the earlier geometry; AC-loss runs stay disabled.
* Known limitations found by these runs (not fixed here): the lamination stress peak of the
  sleeved rotors is not mesh-converged (sharp corners of the Maxwell geometry, magnet bearing
  on the outer lug), and the bridgeless mechanical model does not converge reliably (loose
  magnets in open-clearance contact; default 0.5 mm sleeve diverges).

### Figures and README (earlier in this release)

Self-explanatory, publication-style figures and a rewritten README.

* `sipm._style` (new): house style for every plot - serif fonts, axis titles larger than tick
  labels, series told apart by line style and open markers (readable in grayscale), legends
  and notes placed where they cover no data, figure widths for a double-column paper.
* `sipm.post`: all plots restyled; `mtpa_lines`, `stress_lines`, `stress_map` and
  `kac_lines` take `note=` (operating conditions written inside the axes), `mtpa_lines` and
  `kac_lines` take `labels=`; new `series_label`, `conditions_note` and `topology_figure`.
* `sipm._geomfig` (new): `geometry_plot` labels the features of the pole with callouts
  (`heading=`, `callouts=`, `legend=`, `note=`, drawing into a caller's axes);
  `topology_figure` draws the three topologies side by side with panel headings (a)-(c),
  optional `headings=` and `callouts=`, and one shared legend.
* `sipm.cli`: `plot_results` adds the shared operating conditions as a note and closes each
  figure after saving it.
* Robustness: plots accept empty, single-row and all-NaN tables and tables with a repeated
  index; `stress_map` raises a clear ValueError for collinear points; inputs are not mutated.
* `docs/make_figures.py` regenerates the README figures from the sample CSV files; new
  `docs/img/geometry_parametric.png` and `examples/data/sample_em_conventional.csv` (at that
  time torque values from `docs/validation.md`; replaced by a solved table, see above).
* README: "How it works" diagram, "Parametric geometry" section, captions for the example
  results, BibTeX citation.
* Tests: `tests/test_post_edge.py` (edge cases above, backend selection, figure lifetime,
  topology figure headings and callouts) and further cases in `tests/test_post_cli.py`.

Behaviour changes:

* Figures are saved at 300 dpi (was 160).
* `stress_map` draws no title unless `title=` is given (was "Feasible (hatched): ...").
* A series with zero sleeve thickness is labelled "No sleeve".
* The plotting code no longer probes matplotlib GUI backends (Qt/Tk); it switches to Agg
  unless a non-interactive backend is already selected or `SIPM_INTERACTIVE` is set.
* `labels` keys that match no plotted series raise a UserWarning (were silently ignored).

## v0.1.0

First public release.

* `sipm.params`: flat, validated parameter set; every field tagged as published in the
  companion paper or as an illustrative default.
* `sipm.geometry`: licence-free pole geometry for the conventional, sleeved bridge-retained
  and sleeved bridgeless double-V rotors (constant physical air gap, fillets, feasibility
  checks) and the hairpin stator.
* `sipm.maxwell`: Maxwell 2D transient driver through PyAEDT (one pole, odd periodicity,
  current-driven windings with a current advance angle, average torque, waveforms, no-load
  flux linkage/back-EMF, solid-conductor AC loss factor).
* `sipm.mapdl`: APDL deck writer, batch runner and results parser for the centrifugal +
  interference stress analysis (orthotropic sleeve, contact offset, magnets in contact).
* `sipm.sweep`: design tables, per-point cache, resume after interruption, failure records.
* `sipm.memory`: free-RAM gate before each point, solver RSS watchdog with optional cap,
  periodic AEDT restart.
* `sipm.post`: torque scatter, torque against current angle, stress against interference,
  interpolated stress map with strength contour and feasible region.
* `sipm.analytic` and `sipm.validation`: closed-form disc + sleeve solution and checks.
* Command line interface `sipm` (TOML sweep files) and a thin MATLAB front end (`matlab/`).
