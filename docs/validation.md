# Validation

All numbers below were measured with this code, with the geometry rebased on the author's
Maxwell models (docs/methods.md section 1), on one Windows 11 workstation with Ansys
Electronics Desktop 2024 R1 (PyAEDT 0.24.1) and Mechanical APDL 2024 R1, using the default
geometry unless stated. They describe the behaviour of the toolchain; they are not results of
the companion paper and were not tuned towards it. Re-run them with
`python examples/05_validation.py`, `pytest -m ansys -s` and the configurations in
`examples/configs`. Several solver runs shared the workstation, so wall times only indicate
the order of magnitude. Section 9 still describes the geometry before the rebase (see there).

The package builds the two sleeved rotors only; the conventional rotor of earlier versions
was removed (docs/methods.md section 1) and is not validated here.

Stress quantities (see [methods.md](methods.md) section 3): `svm_lam_max` = `svm_rotor_max`
is the peak von Mises stress over the named load-carrying sections of the lamination
(`svm_<section>` each); `svm_lam_max_all` (global lamination peak) and `svm_pm_max` (magnet
peak) sit at sharp corners and contact edges and are mesh-dependent (section 3). Magnets are
bonded to their pockets unless stated (`magnet_contact = "bonded"`, as the original Workbench
model).

## 1. Mechanical path against a closed-form solution

Case: solid steel disc (radius 71.1 mm, E = 200 GPa, nu = 0.3, 7850 kg/m^3) with a 1.5 mm
ring, 0.1 mm radial interference, 20000 rpm, plane stress, 22.5 deg sector with symmetry
conditions, frictionless contact with the interference as contact offset, element size
0.5 mm. The reference is the closed-form solution in `sipm.analytic` (rotating disc + ring,
cylindrically orthotropic ring, `u = A r^k + B r^-k + C r^3`). That solution is itself
checked against the textbook isotropic formulas (agreement 1e-13) and against a numerical
boundary-value solution of the orthotropic ring (agreement 1e-12) in
`tests/test_analytic.py`. This case does not use the rotor geometry; the rerun gave the same
numbers as before the rebase.

**Isotropic ring** (E = 70 GPa, nu = 0.33, 2700 kg/m^3):

| quantity | FE | closed form | error |
|---|---|---|---|
| disc max von Mises (MPa) | 70.8333 | 70.8323 | +0.001 % |
| ring max hoop stress (MPa) | 108.5503 | 108.5559 | -0.005 % |
| interface pressure, max over contact (MPa) | 0.9718 | 0.9719 | -0.007 % |
| interface pressure, min over contact (MPa) | 0.9717 | 0.9719 | -0.017 % |

**Orthotropic carbon-fibre ring** (default sleeve constants: E_hoop = 125 GPa,
E_radial = 8.8 GPa, nu_hoop,radial = 0.28, 1800 kg/m^3):

| quantity | FE | closed form | error |
|---|---|---|---|
| disc max von Mises (MPa) | 68.6568 | 68.6524 | +0.006 % |
| ring max hoop stress (MPa) | 192.4089 | 192.5863 | -0.092 % |
| interface pressure, max over contact (MPa) | 3.1483 | 3.1518 | -0.111 % |
| interface pressure, min over contact (MPa) | 3.1482 | 3.1518 | -0.114 % |

The orthotropic case also confirms the material orientation. The sleeve's hoop stress
(192 MPa) is that of a hoop-stiff ring. With the two moduli exchanged (radial-stiff ring) the
closed form gives 41.2 MPa and the ring lifts off the disc at speed (zero interface
pressure), so a mis-oriented element system would not pass this check.

**Thermal: orthotropic ring, whole rotor 100 K above the fitting temperature**
(`ring_check("thermal")`: the orthotropic case above with `mech_temperature` = 120 degC,
`mech_ref_temperature` = 20 degC, at which the 0.1 mm interference is defined; default CTEs:
steel 12e-6 /K, sleeve 0 /K along the fibre (hoop) and 30e-6 /K radially). The closed form
adds the thermal particular term to `sipm.analytic` (module docstring); it is checked in
`tests/test_thermal_sleeve.py` against the equilibrium ODE and the constitutive law by finite
differences, the interface and OD conditions, a numerical boundary-value solution of a free
heated orthotropic ring, an isotropic free expansion (zero stress) and the k = 1 limit.

| quantity | FE | closed form | error |
|---|---|---|---|
| disc max von Mises (MPa) | 65.5003 | 65.4924 | +0.012 % |
| ring max hoop stress (MPa) | 339.9660 | 340.3190 | -0.104 % |
| interface pressure, max over contact (MPa) | 6.3048 | 6.3118 | -0.111 % |
| interface pressure, min over contact (MPa) | 6.3046 | 6.3118 | -0.114 % |

The FE path agrees as closely as in the isothermal cases. Heating by 100 K doubles the
interface pressure (3.15 to 6.31 MPa) and raises the sleeve hoop stress by 77 % (192.6 to
340.3 MPa): the steel grows by 12e-6 x 100 x 71.1 = 0.085 mm while the fibre direction does
not, close to the 0.1 mm fitted interference. The sleeve node export (`sipm_sleeve_nodes.txt`,
`r["_sleeve_nodes"]`) reproduces the scalar results exactly in all three ring cases; its
radial tension (at most 0.003 MPa) and r-theta shear (at most 0.5 MPa, 0.15 % of the hoop
stress) are the discretisation noise of an axisymmetric problem, whose exact values are zero.

## 2. IPM rotor stress, one point per topology (20000 rpm, default mesh)

Topology defaults (saved states of the Maxwell models): sleeved bridge-retained with a 1.5 mm
sleeve, sleeved bridgeless with a 0.5 mm sleeve, 0.1 mm interference, frictionless sleeve
contact, bonded magnets (the designed 0.1 mm clearance bridged by the bond). Integration tests
(`pytest -m ansys`); both converge.

| topology | headline `svm_lam_max` (MPa) | section, location r, theta | global lamination peak (MPa), location | magnet peak (MPa) | max sleeve hoop (MPa) | open sleeve contact | nodes | wall |
|---|---|---|---|---|---|---|---|---|
| sleeved_bridged | 253.3 | `qhole`, 54.42 mm, 22.35 deg | 691.1, 61.15 mm, 2.06 deg | 199.6 | 223.9 | 0 % | 25 160 | 22.7 s |
| sleeved_bridgeless | 318.1 | `qhole`, 54.42 mm, 22.35 deg | 502.6, 60.33 mm, 1.23 deg | 220.9 | 236.5 | 0 % | 21 633 | 16.8 s |

These numbers were re-run after the thermal load case and the sleeve node export were added:
with the default (equal) temperatures the deck writes no thermal commands and every result
above is bit-identical. The bridgeless row was re-run once more after its thin sleeve's bore
arc was sized like the lamination OD (methods.md section 3, item 2): the headline, global and
magnet peaks moved by less than 0.01 MPa, the sleeve hoop stress from 236.6 to 236.5 MPa, the
contact pressure min / max from 0.850 / 4.784 to 0.832 / 4.865 MPa, 21 932 to 21 633 nodes. The
new sleeve columns at these points: max radial tension `sr_sleeve_max` 0.02 MPa
(bridge-retained) and 0.42 MPa (bridgeless; 0.23 before the bore-arc sizing), max |r-theta
shear| `tau_sleeve_max` 1.6 and 3.8 MPa (bridgeless 6.0 before). With the bridge-retained rotor 100 K above the fitting
temperature (default CTEs) the sleeve hoop stress rises from 223.9 to 356.6 MPa, the headline
`svm_lam_max` from 253.3 to 274.3 MPa and the peak contact pressure from 4.57 to 8.42 MPa
(`test_mech_point_thermal`).

Section peaks (MPa): bridge-retained `bridge1` 103.9, `bridge2` 97.9, `rib1` 13.0, `web2` 76.4,
`qhole` 253.3; bridgeless `lip1_pole` 6.7, `lip1_q` 4.0, `lip2_pole` 3.4, `lip2_q` 3.6, `qhole`
318.1. The headline peak is at the apex fillet of the q-axis lightening hole, on the q-axis,
in both rotors. The global peaks sit at sharp corners of the slot ends next to the bonded
magnets (bridge-retained: a corner of the layer-1 inner relief; bridgeless: near the inner
end of the layer-1 slot).

**Frictional magnets (bridge-retained rotor, default point).** With
`magnet_contact = "frictional"` (open 0.1 mm clearance, frictionless) the solution converges
with a stabilisation / strain energy ratio of 1.2e-6. With the stabilisation switched off
(KEYOPT(12,15) = 1) it converges too, and every reported quantity agrees with the stabilised
run within 3e-7 (relative): the stabiliser does not change the result. The loose magnets load the lamination much more than bonded ones:
`bridge1` 192.6, `bridge2` 304.0, `rib1` 146.8, `web2` 248.9, `qhole` 330.8 MPa (headline, at
the inner fillet of the lightening hole, 42.15 mm, 13.72 deg), magnet peak 719.5 MPa, sleeve
hoop stress 285.8 MPa. For the bridgeless rotor the frictional option is refused with an
error: in Phase-2 runs with this contact, two of six sweep points and the default point
diverged (a loose magnet's displacement grew without limit), also with magnet friction 0.2,
with the clearance ignored, or with a 1 mm pinball.

## 3. Mesh refinement of the named sections

Defaults, bonded magnets, 0.1 mm interference, 20000 rpm; global / fine element size.

**Sleeved bridge-retained rotor:**

| mesh (mm) | nodes | `bridge1` | `bridge2` | `rib1` | `web2` | `qhole` | global lamination peak | magnet peak | max sleeve hoop | wall |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 / 0.3 | 9 693 | 104.1 | 97.4 | 12.5 | 73.6 | 252.0 | 441.4 | 187.2 | 222.9 | 9 s |
| 0.5 / 0.15 (default) | 25 160 | 103.9 | 97.9 | 13.0 | 76.4 | 253.3 | 691.1 | 199.6 | 223.9 | 20 s |
| 0.3 / 0.08 | 73 496 | 105.0 | 99.0 | 12.9 | 79.5 | 253.7 | 669.2 | 306.0 | 224.4 | 53 s |
| 0.2 / 0.04 | 200 675 | 106.5 | 99.8 | 12.9 | 80.9 | 254.7 | 823.9 | 682.6 | 225.2 | 138 s |
| change between the two finest | | +1.4 % | +0.9 % | +0.4 % | +1.8 % | +0.4 % | +23 % | +123 % | +0.4 % | |

**Sleeved bridgeless rotor:**

| mesh (mm) | nodes | `lip1_pole` | `lip1_q` | `lip2_pole` | `lip2_q` | `qhole` | global lamination peak | magnet peak | max sleeve hoop | wall |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 / 0.3 | 9 497 | 4.7 | 3.8 | 2.9 | 3.9 | 315.9 | 359.7 | 206.9 | 235.3 | 11 s |
| 0.5 / 0.15 (default) | 21 633 | 6.7 | 4.0 | 3.4 | 3.6 | 318.1 | 502.6 | 220.9 | 236.5 | 18 s |
| 0.3 / 0.08 | 56 188 | 8.2 | 4.6 | 4.2 | 4.4 | 321.0 | 614.7 | 282.3 | 237.8 | 40 s |
| 0.2 / 0.04 | 136 117 | 8.9 | 5.1 | 4.3 | 4.2 | 323.5 | 849.2 | 655.0 | 239.0 | 64 s |
| change between the two finest | | +7.5 % | +10.4 % | +3.1 % | -4.5 % | +0.8 % | +38 % | +132 % | +0.5 % | |

(bridgeless rows re-run with the bore arc of the thin sleeve sized like the lamination OD,
methods.md section 3, item 2; the bridge-retained rows are unchanged by it)

(stresses in MPa)

* **Converged (< 3 % between the two finest meshes):** every section of the bridge-retained
  rotor (0.4-1.8 %), the bridgeless `qhole` (0.8 %), the sleeve hoop stress of both rotors,
  and with them the headline `svm_lam_max` of both rotors (`qhole` at every mesh).
* **Not converged: the four bridgeless lips** (3-10 %). Their peaks lie at the end of the
  section next to the lip root, a sharp corner, and keep rising as the mesh resolves its
  field. They are small (below 9 MPa, under 3 % of the headline value), so they do not affect
  `svm_lam_max`, but they are not reliable on their own.
* The global lamination peak and the magnet peak grow without limit (sharp corners and the
  edges of the bonded faces); they are reported but not used for design.

## 4. No-load flux density (generic BH curve)

`MaxwellSession.no_load_flux_density`: zero current, rotor d-axis on the middle tooth (t = 0),
|B| at the middle of the d-axis tooth (r = 82.7 mm) and of the yoke near the q-axis
(r = 97.3 mm), and the maximum over each lamination (a point value at saturated bridges or
corners):

| topology | mid-tooth (T) | mid-yoke (T) | max in stator (T) | max in rotor lamination (T) |
|---|---|---|---|---|
| sleeved_bridged | 1.01 | 0.68 | 1.33 | 2.49 |
| sleeved_bridgeless | 1.68 | 1.31 | 2.07 | 2.45 |

The bridgeless rotor, with no bridges or ribs and the thinnest sleeve, puts the most magnet
flux into the stator; the bridge-retained rotor's 1.5 mm sleeve widens its magnetic gap. The
stator (Maxwell model) has a 15.3 mm yoke and teeth at least 4.6 mm wide.

## 5. Electromagnetic sign conventions and periodicity

No-load runs, defaults, 3000 rpm, 360 electrical degrees in 2-degree steps:

| topology | flux-linkage phase A / B / C (deg) | axis predicted (deg) | EMF fundamental A / B / C (V peak) | half-wave asymmetry (max) | EMF THD | cogging pk-pk / mean (N m) |
|---|---|---|---|---|---|---|
| sleeved_bridged | -150.0000 / -30.0000 / 90.0000 | -150 / -30 / 90 | 99.240 / 99.239 / 99.244 | 7.5e-8 | 6.0 % | 2.87 / 0.000 |
| sleeved_bridgeless | -149.9996 / -29.9996 / 90.0004 | -150 / -30 / 90 | 179.778 / 179.783 / 179.756 | 4.5e-8 | 12.5 % | 6.63 / -0.224 |

* The flux-linkage phase of every phase equals the axis computed from the winding layout to
  within 1e-5 deg (bridge-retained) and 4e-4 deg (bridgeless), rotor d-axis at 0 at t = 0.
  This confirms the magnet polarity, the coil polarity convention and the phase axes used to
  place the current vector, in both topologies.
* The flux linkages are odd half-wave symmetric to 1e-7, as required by the master/slave
  boundaries with B_slave = -B_master over one pole.
* EMF amplitudes agree within 1.5e-4 (relative); EMF phase shifts B-A = 120.00-120.01 deg and
  C-A = 240.00 deg (positive sequence for the positive rotation). Each EMF phase
  differs by 1.0 deg from the phase of d(psi)/dt, half a time step of the backward difference.
* The EMF amplitude equals omega x psi_1 (bridge-retained: 2 pi 200 Hz x 0.07898 Wb =
  99.25 V).

## 6. Torque against current advance angle (MTPA), one point per topology, ampere-turns

Sleeved bridge-retained rotor (defaults, 1.5 mm sleeve), 420 A rms, 3000 rpm, average over 60
electrical degrees:

| gamma (deg) | 90 | 100 | 110 | 120 | 125 | 130 | 135 | 140 | 150 | 160 | 170 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| torque (N m) | 250.54 | 286.62 | 315.28 | 332.74 | 336.01 | 334.37 | 326.55 | 312.04 | 262.91 | 191.06 | 101.18 |

Torque is positive (motoring, in the direction of rotation) over the whole range. It peaks
near 125 deg, i.e. with negative d-axis current, as expected from the reluctance torque of an
IPM rotor.

One point per topology at gamma = 130 deg, 420 A rms (integration tests):

| topology | average torque (N m) | torque ripple pk-pk / average | wall time incl. model build |
|---|---|---|---|
| sleeved_bridged | 334.37 | 0.122 | 71.0 s |
| sleeved_bridgeless | 497.66 | 0.158 | 91.6 s |

Every torque waveform of the sweeps in section 7 returns to its start value after the 60 electrical
degrees of the run (largest difference 2.7e-6 of the average), as a balanced three-phase
machine must; all averages are positive.

Mesh of the torque models (AEDT mesh statistics): 17 846 (bridge-retained) and 18 044
(bridgeless) elements per pole, rotor element edges up to 1.2 mm (`mesh_rotor`), air-gap
layers 0.3 mm.

**Ampere-turns.** The coil region of each slot is the model's 8 touching rectangles
(4.3 x 1.57 mm, 54.0 mm^2) carrying `hairpin_layers` = 7 turns, with 2 parallel paths: at
420 A rms the slot carries 7 x 210 = 1470 A rms (2079 A peak), and a phase has
16 x 7 / (2 x 2) = 28 series turns. Check against the no-load flux linkage psi_1 (section 5):
with pure q-axis current (gamma = 90 deg) the magnet torque is 1.5 x p x psi_1 x I_peak:

| topology | psi_1 (Wb) | 20 A rms: FE / estimate | 420 A rms: FE / estimate |
|---|---|---|---|
| sleeved_bridged | 0.07898 | 13.46 / 13.40 N m = 1.004 | 250.54 / 281.47 N m = 0.890 |
| sleeved_bridgeless | 0.14307 | 24.16 / 24.28 N m = 0.995 | 399.11 / 509.89 N m = 0.783 |

At 20 A rms the FE torque matches the estimate within 0.5 %, so the flux linkage seen by the
windings and the ampere-turns that make the torque use the same 28 series turns (a 7-turn
current in an 8-turn flux linkage, or the reverse, would give 0.875 or 1.14). At 420 A rms
saturation lowers the torque by 11-22 % below the linear estimate. These are internal
consistency checks of the Maxwell model; `sipm.analytic` has no EM estimate. 7 turns over 8 drawn rectangles is therefore
self-consistent for the field solution: Maxwell spreads the 7 turns' current uniformly over
the 54.0 mm^2 region, so only the current density (27.2 A/mm^2 rms instead of 31.1 A/mm^2 in
seven 6.75 mm^2 conductors) and its distribution over the slot depth differ from a
seven-conductor layout, which changes slot leakage slightly but not the ampere-turns. The DC
loss formula (`dc_copper_loss`) uses seven conductors of 6.75 mm^2. The eighth rectangle has no
counterpart in the turn count; whether the model meant 8 conductors (one of them for a
different purpose) cannot be told from the projects.

## 7. Sweeps through the runner

The runner's interruption and resume (stop after N points, hard kill of the process tree)
were demonstrated with the geometry before the rebase; they do not depend on the geometry and
were not repeated here. All sweeps below ran without failures unless stated.

**EM, sleeved bridge-retained** (`em_sleeved_bridged.toml`, 420 A rms). Average torque (N m):

| gamma (deg) | 1.0 mm sleeve | 2.0 mm sleeve |
|---|---|---|
| 110 | 340.88 | 289.90 |
| 120 | 362.62 | 303.30 |
| 130 | 368.75 | 299.85 |
| 140 | 352.40 | 275.10 |
| 150 | 302.92 | 229.31 |

Wall time 73-138 s per point.

**EM, sleeved bridgeless** (`em_sleeved_bridgeless.toml`, 0.5 mm sleeve): 489.2, 502.1 and
494.3 N m at gamma = 115, 125 and 135 deg with a 0.5 mm lip depth; 479.1, 493.5 and 488.3 N m
with a 1.5 mm lip depth.

**Mechanical, sleeved bridge-retained** (`mech_sleeved_bridged.toml`, 20000 rpm, bonded
magnets, all 15 points converged). Section peaks (MPa); the headline `svm_lam_max` is `qhole`
at every point:

| sleeve (mm) | interference (mm) | `bridge1` | `bridge2` | `rib1` | `web2` | `qhole` = headline | max sleeve hoop | open sleeve contact |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.0 | 102.7 | 100.8 | 18.6 | 88.5 | 278.8 | 47.1 | 0.719 |
| 1.0 | 0.1 | 103.9 | 99.0 | 15.2 | 81.5 | 264.3 | 224.4 | 0 |
| 1.0 | 0.2 | 110.6 | 105.4 | 11.8 | 74.5 | 249.7 | 398.5 | 0 |
| 1.0 | 0.3 | 117.6 | 111.9 | 8.4 | 67.6 | 235.2 | 574.0 | 0 |
| 1.0 | 0.4 | 124.7 | 118.6 | 7.1 | 60.6 | 220.6 | 749.6 | 0 |
| 1.5 | 0.0 | 101.8 | 98.2 | 18.1 | 86.9 | 275.2 | 46.5 | 0.714 |
| 1.5 | 0.1 | 103.9 | 97.9 | 13.0 | 76.4 | 253.3 | 223.9 | 0 |
| 1.5 | 0.2 | 113.3 | 106.1 | 7.8 | 65.9 | 231.3 | 399.7 | 0 |
| 1.5 | 0.3 | 122.8 | 114.7 | 9.2 | 55.3 | 209.4 | 576.0 | 0 |
| 1.5 | 0.4 | 132.2 | 123.2 | 14.0 | 44.8 | 187.4 | 752.4 | 0 |
| 2.0 | 0.0 | 101.9 | 97.2 | 17.9 | 85.5 | 271.8 | 46.3 | 0.732 |
| 2.0 | 0.1 | 104.6 | 96.7 | 10.9 | 71.3 | 242.4 | 223.9 | 0 |
| 2.0 | 0.2 | 116.0 | 106.9 | 7.9 | 57.2 | 213.0 | 401.0 | 0 |
| 2.0 | 0.3 | 127.4 | 117.0 | 14.4 | 43.0 | 183.5 | 578.2 | 0 |
| 2.0 | 0.4 | 138.8 | 127.1 | 20.9 | 28.8 | 154.1 | 755.4 | 0 |

`qhole` and `web2` fall as the interference rises, more steeply with a thicker sleeve; the
outer bridges rise slowly with it. The sleeve hoop stress is set mainly by the interference. Without
interference about 72 % of the sleeve contact elements open at speed. With a 400 MPa example
yield strength every point meets the stress limit, so the feasible region is bounded by the
interference limit (interference <= 0.2 x sleeve thickness) alone. Global lamination peak
397-756 MPa and magnet peak 125-217 MPa (mesh-dependent). Wall time 14-26 s per point.

**Mechanical, sleeved bridgeless** (`mech_sleeved_bridgeless.toml`, 20000 rpm, bonded
magnets, all 6 points converged; `examples/data/sample_mech_bridgeless.csv`). Section peaks
(MPa); the headline is `qhole` at every point:

| sleeve (mm) | interference (mm) | `lip1_pole` | `lip1_q` | `lip2_pole` | `lip2_q` | `qhole` = headline | max sleeve hoop | open sleeve contact |
|---|---|---|---|---|---|---|---|---|
| 1.5 | 0.0 | 3.8 | 7.2 | 5.3 | 13.1 | 320.0 | 53.9 | 0.885 |
| 1.5 | 0.2 | 25.5 | 30.5 | 21.3 | 37.3 | 269.6 | 414.8 | 0 |
| 1.5 | 0.4 | 44.2 | 72.5 | 42.5 | 87.7 | 219.0 | 784.8 | 0 |
| 2.0 | 0.0 | 3.4 | 6.4 | 5.0 | 11.3 | 317.0 | 52.4 | 0.865 |
| 2.0 | 0.2 | 28.3 | 38.6 | 25.7 | 35.5 | 249.0 | 420.8 | 0 |
| 2.0 | 0.4 | 49.0 | 87.0 | 52.0 | 81.3 | 180.9 | 796.7 | 0 |

With bonded magnets the three lamination pieces are tied together by the magnets, and the
q-axis body carries the hoop load as in the bridge-retained rotor. The interference loads the
lips (sleeve pressure on the thin lips over the slot mouths), whose values are not
mesh-converged (section 3). Global lamination peak 337-533 MPa and magnet peak 130-224 MPa
(mesh-dependent). Wall time 13-20 s per point.

Figures are written to `<workdir>/<sweep>/figures/`. A selection is in the README
(`docs/img`).

## 8. Memory gate and solver memory (observed)

Peak RSS of the solver process tree started by the runner (`solver_rss_peak_gb`):

* EM points (one-pole Maxwell model, ansysedt + solver processes): 0.862-0.889 GB over the
  two sleeved EM sweeps (16 points).
* MAPDL points (default mesh, bonded magnets): 0.47-0.52 GB (bridge-retained rotor) and
  0.39-0.40 GB (bridgeless rotor).

The memory gate (4 GB free by default) and the AEDT restart every 10 points were exercised on
the earlier geometry; they did not fire in these runs.

## 9. AC copper loss factor (solid conductors)

> **Note:** this section was obtained with the geometry before the rebase on the author's
> Maxwell models (v0.1.0 parameterisation) and has not been recomputed: AC-loss runs are
> disabled until the strand geometry of the rebased stator is defined.

`examples/configs/acloss_example.toml`: sleeved bridge-retained rotor, 420 A rms,
gamma = 125 deg, 1.5 electrical periods with 90 steps per period, loss averaged over the last
period, conductor mesh 0.3 mm, copper 5.8e7 S/m (20 C), active length only:

| speed (rpm) | p_ac (W) | p_dc (W) | Kac | average torque (N m) | wall |
|---|---|---|---|---|---|
| 1000 | 4237.7 | 4094.2 | 1.035 | 435.13 | 527 s |
| 10000 | 17401.5 | 4094.2 | 4.250 | 433.31 | 406 s |
| 20000 | 51409.7 | 4094.2 | 12.557 | 432.73 | 491 s |

Kac tends to 1 at low speed, as it should. The average torque of these solid-conductor runs
matched the stranded torque run at the same angle of that geometry, a consistency check of the
excitation. The high values at speed came from that illustrative winding: seven unsegmented
2.6 mm tall conductors per slot, with a 1.8 mm skin depth in copper at 1333 Hz.

## 10. Composite sleeve

**MAPDL and angle-ply constants.** Parameter validation requires a symmetric positive
definite 6x6 sleeve compliance instead of bounding each Poisson ratio (methods.md section
4.2). MAPDL checks the same property: a deck with `sleeve_nu_hoop_radial` = 4 (forced past the
validation) stops with "The stress-strain matrix of material 3 is not positive definite,
which is required for real materials", while a +/-30 deg carbon_sm/epoxy winding
(`nu_hoop_axial` = 1.46, `nu_hoop_radial` = -0.13) solves, with no material warning in the
output of the disc + ring model; all points below (up to +/-45 deg) solve as well.

**Composite spec against the Table II constants** (bridge-retained rotor, defaults otherwise;
carbon_sm/epoxy, Chamis, vf = 0.536, hoop winding: E_hoop 124.9, E_radial 7.98 GPa,
nu 0.270, 1522 kg/m3 against Table II's 125, 8.8 GPa, 0.28, 1800 kg/m3; methods.md
section 4.5):

| quantity | Table II constants | composite spec | difference |
|---|---|---|---|
| headline `svm_lam_max` (`qhole`) (MPa) | 253.28 | 252.51 | -0.30 % |
| `bridge1` / `bridge2` (MPa) | 103.91 / 97.88 | 104.73 / 98.75 | +0.8 % / +0.9 % |
| `rib1` / `web2` (MPa) | 12.96 / 76.40 | 12.77 / 76.03 | -1.4 % / -0.5 % |
| magnet peak (MPa) | 199.60 | 199.05 | -0.27 % |
| max sleeve hoop (MPa) | 223.92 | 222.99 | -0.41 % |
| min sleeve radial (MPa) | -4.63 | -4.68 | +1.1 % |
| contact pressure max / min (MPa) | 4.57 / 2.58 | 4.62 / 2.75 | +1.1 % / +6.8 % |
| max-stress index / Tsai-Wu safety factor | 0.112 / 8.92 | 0.111 / 8.94 | -0.4 % / +0.2 % |
| Tsai-Wu index `sleeve_fi_tsai_wu` (= 1 / safety factor) | 0.1122 | 0.1119 | -0.2 % |

The lower density of the composite spec (centrifugal load on the sleeve -15 %) and its lower
transverse stiffness move the results by about 1 %, the minimum contact pressure by 7 %. At
the default point (Table II constants, isothermal) the sleeve failure columns are: maximum
stress index 0.112 (hoop stress against Xt = 2000 MPa), Hashin fibre 0.0125, Hashin matrix
1.7e-5, Tsai-Wu safety factor 8.92 at the sleeve bore on the d-axis (r = 71.1 mm,
theta = 0) and Tsai-Wu index 1/8.92 = 0.112. The column is the reciprocal strength ratio,
not the Tsai-Wu quadratic, which is negative at this point (-0.050) because compression of
the matrix directions dominates its linear terms.

**Sweep `examples/configs/mech_composite_sleeve.toml`** (bridge-retained rotor, 20000 rpm,
0.1 mm interference fitted at 20 degC, carbon_sm/epoxy, vf = 0.55, Chamis; all 8 points
converged, 12.5-20.5 s each; stresses in MPa):

| winding (deg) | rotor (degC) | `svm_lam_max` | max sleeve hoop | min / max sleeve radial | max \|tau_rt\| | contact p min / max | open sleeve contact | max-stress index | Hashin fibre / matrix | Tsai-Wu safety factor (r mm, theta deg); the Tsai-Wu index is its reciprocal |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 20 | 251.8 | 228.5 | -4.83 / 0.02 | 1.54 | 2.81 / 4.76 | 0 | 0.114 | 0.013 / 1.8e-5 | 8.72 (71.1, 0) |
| 0 | 120 | 272.6 | 364.1 | -8.64 / 0.04 | 2.49 | 5.26 / 8.57 | 0 | 0.182 | 0.033 / 4.3e-5 | 5.49 (72.6, 22.5) |
| 15 | 20 | 257.1 | 186.7 | -3.80 / 0.02 | 1.30 | 2.15 / 3.77 | 0 | 0.105 | 0.021 / 1.1e-5 | 7.12 (71.1, 0) |
| 15 | 120 | 278.1 | 318.9 | -7.49 / 0.03 | 2.06 | 4.60 / 7.43 | 0 | 0.170 | 0.034 / 2.5e-5 | 4.47 (72.6, 22.5) |
| 30 | 20 | 270.7 | 77.8 | -1.32 / 0.01 | 0.68 | 0.16 / 1.32 | 0 | 0.228 | 0.054 / 1.6e-6 | 5.24 (71.1, 0) |
| 30 | 120 | 299.7 | 145.5 | -3.11 / 0.02 | 0.91 | 1.46 / 3.08 | 0 | 0.275 | 0.084 / 2.5e-6 | 3.28 (72.6, 21.3) |
| 45 | 20 | 275.8 | 46.1 | -0.82 / 4.86 | 1.76 | 0 / 0.76 | 0.857 | 0.330 | 0.109 / 0.117 | 2.77 (71.1, 18.9) |
| 45 | 120 | 313.2 | 38.2 | -0.41 / 0.79 | 0.47 | 0 / 0.37 | 0.518 | 0.273 | 0.076 / 2.5e-4 | 4.17 (72.6, 12.8) |

* Winding away from the hoop direction softens the sleeve in hoop (E_hoop 128, 103, 42 and
  14 GPa at 0, 15, 30 and 45 deg), so it carries less hoop stress and presses less on the
  lamination; the lamination stress rises (251.8 to 275.8 MPa cold). At 45 deg the sleeve
  lifts off over most of the contact at speed (86 % open); the sleeve then carries radial
  tension (up to 4.9 MPa) and its Hashin matrix index rises to 0.12.
* Heating by 100 K raises the sleeve hoop stress and the contact pressure in the hoop and
  15-30 deg windings, because the steel expands more than the sleeve in hoop (sleeve hoop CTE
  0.24, -1.9 and -5.5e-6 /K against 12e-6 /K). At 45 deg the sleeve hoop CTE (3.2e-6 /K) is
  still below the steel's: heating closes part of the gap (open fraction 86 % to 52 %).
* The ply failure indices stay well below 1 at every point. The safety factor falls from 8.7
  (hoop, cold) to 2.8 (45 deg, cold) and 3.3 (30 deg, hot). In the angle plies the
  fibre-direction stress is no longer the hoop stress, so the max-stress index rises with the
  angle although the homogenised hoop stress falls.
* The +/-theta plies of the hot points carry the ply-to-ply thermal mismatch as well
  (methods.md section 4.4); cure residual stress is not included.

## 11. Load cases and lift-off speed

**Six load cases, default bridge-retained point** (`mech_cases` = all six; 20000 rpm rated,
24000 rpm overspeed (`overspeed_factor` 1.2), cold 20 degC = the fitting temperature, hot
`mech_temperature` = 120 degC; 0.1 mm interference, frictionless sleeve, bonded magnets;
stresses in MPa). Six MAPDL runs, 12.3-14.4 s each, 80 s for the point:

| case | `svm_lam_max` (section) | magnet peak | max sleeve hoop | min / max sleeve radial | contact p min / max | open sleeve contact | max-stress index | Tsai-Wu safety factor |
|---|---|---|---|---|---|---|---|---|
| `standstill_cold` | 37.3 (`bridge2`) | 19.3 | 176.4 | -4.40 / 0.014 | 2.620 / 4.30 | 0 | 0.088 | 11.25 |
| `standstill_hot` | 115.6 (`bridge2`) | 225.1 | 312.4 | -8.57 / 0.041 | 5.074 / 8.47 | 0 | 0.156 | 6.38 |
| `rated_cold` | 253.3 (`qhole`) | 199.6 | 223.9 | -4.63 / 0.020 | 2.576 / 4.57 | 0 | 0.112 | 8.92 |
| `rated_hot` | 274.3 (`qhole`) | 369.2 | 356.6 | -8.50 / 0.041 | 4.970 / 8.42 | 0 | 0.178 | 5.61 |
| `overspeed_cold` | 374.4 (`qhole`) | 294.6 | 247.6 | -4.85 / 0.026 | 2.234 / 4.82 | 0 | 0.124 | 8.08 |
| `overspeed_hot` | 395.4 (`qhole`) | 439.0 | 376.1 | -8.58 / 0.042 | 4.781 / 8.40 | 0 | 0.188 | 5.31 |

`rated_cold` reproduces the single default point of section 2 (223.920 MPa hoop, 253.3 MPa
headline) and `rated_hot` the 100 K point there. The envelope: `svm_lam_max__worst` 395.4,
`sh_sleeve_max__worst` 376.1, `sr_sleeve_max__worst` 0.042, every failure index and
`sleeve_sf__worst` 5.31, all governed by `overspeed_hot`; `contact_p_min__worst` 2.234 MPa
by `overspeed_cold`; `sleeve_open_fraction__worst` 0 (tie, `standstill_cold`). Heating raises
the interface pressure at every speed (the steel grows into the sleeve), so the hot cases
load the sleeve most and the cold ones are closest to lift-off. The same point with the case
list reversed gives identical columns; the cases are solved in the same fixed order, and the
only differences, in the failure indices at the 1e-12 relative level, are the run-to-run
round-off of the parallel solver in the node file.

**Lift-off, closed form against FE (disc + ring of section 1).** `ring_liftoff_check`:
orthotropic ring, cases standstill and 20000 rpm at 20 degC and at 120 degC, lift-off
extrapolated element by element (methods.md section 3.1) against `sipm.analytic.liftoff_speed`:

| temperature | FE (rpm) | closed form (rpm) | error |
|---|---|---|---|
| 20 degC (= fitting) | 56 193.5 | 56 193.7 | -3.4e-6 |
| 120 degC | 76 958.7 | 76 958.9 | -2.7e-6 |

The extrapolation is a ratio of pressures, so the 0.11 % discretisation error of the pressure
itself (section 1) cancels.

**Lift-off of the rotors, checked by direct runs.** `liftoff_check`: cases standstill, rated
and overspeed at one temperature, the extrapolated speed, then single-case runs at that speed
and at 0.95 x that speed (min contact pressure over the 56 / 53 sleeve contact elements;
`linearity` = largest deviation of a closed element's FE pressure from its fitted line,
relative to its standstill pressure):

| rotor, temperature | lift-off (rpm), cases used | p min at standstill | at lift-off: p min, elements open | at 0.95 x: p min, elements open | linearity |
|---|---|---|---|---|---|
| bridge-retained, 20 degC | 39 357, standstill + overspeed | 2.620 | 0.026 (1.0 %), 1 | 0.345, 0 | 2.8e-5 at 20000 rpm, 4.0e-4 at 0.95 x |
| bridge-retained, 120 degC | 51 553, standstill + overspeed | 5.074 | 0.044 (0.9 %), 1 | 0.590, 0 | |
| bridgeless, 20 degC | 23 591, standstill + rated (overspeed already open) | 1.136 | 0.398, 1 | 0.460, 1 | 6.6 % at 20000 rpm |
| bridgeless, 120 degC | 0: open at standstill | | | | |

Bridge-retained rotor: the linear model holds. At the extrapolated speed the governing element
(at the d-axis, theta = 0.2 deg; not the element with the lowest pressure at standstill) is the
first to open and the minimum pressure is 1 % of its standstill value; at 0.95 x every element
is closed. The prediction of a closed element's pressure from the standstill and overspeed runs
is good to 3e-5 at 20000 rpm and 4e-4 at 0.95 x the lift-off speed (0.5 % at the lift-off
speed itself, where the first element is partly open), i.e. the lift-off speed is accurate to
about 0.5 %.

**Bridgeless rotor: the linearity does not hold.** Its three lamination pieces are held only by
the sleeve and the bonded magnets, and the OD has openings between them (11.39-12.39 deg and
19.05-20.05 deg). The contact elements at the piece edges next to the openings (theta 11.2, 12.6
and 18.9 deg) are partly open: the pieces tilt under load, so the contact zone at each edge
changes with speed and temperature. Cold, every element is closed at standstill and 20000 rpm and one
edge element (12.6 deg) is partly open at 24000 rpm; the lines through standstill and 24000 rpm
miss the pressures of the closed elements at 20000 rpm by up to 6.6 % of their standstill
pressure, against 3e-5 in the bridge-retained rotor. The extrapolation from standstill and
20000 rpm gives 23 591 rpm, but that edge element is already open at 0.95 x that speed
(22 411 rpm): the true first opening lies below it. Hot (120 degC), three edge
elements (10.8, 11.2 and 18.9 deg) are already open at standstill, so the extrapolation gives 0 by
the definition (first point of the interface to open), although those edges close again with speed
(at 20000 rpm hot every element is closed, minimum 1.42 MPa) and the edge element at 12.6 deg
opens only between 24000 and 42 300 rpm. The pressure is not monotonic in speed here, which no
two-point line can represent. The bridgeless rows were re-run after the bore arc of its thin
sleeve was sized like the lamination OD (methods.md section 3, item 2; before: 23 937 rpm,
1.126 / 0.340 / 0.400 MPa, 8.7 %, two edge elements open hot at standstill), with the NaN rule
below bypassed to obtain them; `liftoff_rpm_cold` / `_hot` are now NaN for `sleeved_bridgeless`, with the reason in
`liftoff_basis_*`. For this rotor read `contact_p_min__<case>` and
`sleeve_open_fraction__<case>` instead.

**Three-point check.** When three or more speeds at one temperature have fully closed contact,
the middle ones are compared with the two-point lines element by element and the basis is
flagged above 1 % of the element's standstill pressure. In the six-case point above
(standstill, rated and overspeed all closed) the rated case deviates by at most 2.8e-5 (cold),
so no flag is set; the bridgeless rotor's 6.6 % would be flagged.

## 12. Layered sleeve and winding prestress

**Closed form** (`sipm.analytic.disc_rings`, methods.md section 4.6), checked offline in
`tests/test_multilayer_sleeve.py`: for N = 1 the coefficients, the pressure and the fields equal
those of `disc_ring` bit for bit (four load sets, one of them separated); three bonded rings
of one material reproduce one ring to 1e-9; in a 3-layer stack (the case below) the
equilibrium equation and the constitutive law with thermal and eigenstrain hold by finite
differences to 1e-6 of the stress scale, u and sigma_r are continuous at both interfaces
(1e-10), the interference condition and the traction-free OD hold, the pressure is the sum of
its interference, thermal, prestress and centrifugal parts (1e-10), and a free ring with a hoop
eigenstrain agrees with a numerical boundary-value solution (1e-9 mm, 1e-6 MPa) with zero net
hoop force.

**Which MAPDL mechanism for the prestress.** On the disc + 1.5 mm Table II ring (no
interference, standstill, hoop eigenstrain eps0 = 100 MPa / 125 GPa) both `INISTATE` forms were
run: an initial elastic hoop strain (`DTYP,EPEL`, eps0) and an initial stress (`DTYP,STRE`,
the plane-stress stiffness times eps0, computed in Python), both in the element axes
(`CSYS,-2`). The two runs gave identical results to every printed digit: contact pressure
2.0720 MPa against 2.0743 MPa closed form (-0.11 %), hoop stress 98.39 / 98.13 MPa at bore / OD
against 98.51 / 98.24 MPa (-0.12 %). The strain form is used (methods.md section 4.6).

**Averaging across the layer interface.** On the 3-layer ring below, a post-only run on the
results file read the hoop stress of the interface node at r = 71.85 mm, theta = 0 with
different element selections: both layers selected 311.28 MPa (also with everything selected),
layer 1 only 445.34 MPa, layer 2 only 177.22 MPa. MAPDL averages nodal results across the
material boundary, which is why the deck post-processes layer by layer.

**3-layer ring against the closed form.** `multilayer_ring_check` (`sipm validate
ring-multilayer`, `test_ring_multilayer_against_closed_form`): solid steel disc 71.1 mm with a
1.5 mm carbon_sm/epoxy sleeve (vf 0.55, Chamis), layers hoop / +/-30 deg / hoop with
thicknesses 0.75 / 0.375 / 0.375 mm, winding prestress 100 / 40 / 100 MPa, 0.1 mm interference,
20000 rpm, rotor 120 degC against a fitting temperature of 20 degC, element size 0.5 mm, 22 163
nodes, 12 s. FE layer values from that layer's node export (mean over the nodes at the radius):

| quantity | FE | closed form | error |
|---|---|---|---|
| disc max von Mises (MPa) | 64.5480 | 64.5390 | +0.014 % |
| interface pressure, max / min over contact (MPa) | 7.2572 / 7.2571 | 7.2653 | -0.110 % / -0.112 % |
| layer 1 (hoop): max hoop, at its bore, at its OD (MPa) | 445.352 / 444.512 / 445.345 | 445.812 / 445.006 / 445.812 | -0.103 / -0.111 / -0.105 % |
| layer 2 (+/-30 deg): max hoop, at its bore, at its OD (MPa) | 177.599 / 177.218 / 177.599 | 177.750 / 177.370 / 177.750 | -0.085 / -0.086 / -0.085 % |
| layer 3 (hoop): max hoop, at its bore, at its OD (MPa) | 447.735 / 447.019 / 447.733 | 448.194 / 447.483 / 448.194 | -0.103 / -0.104 / -0.103 % |
| radial stress at interface 1/2, 2/3 (MPa) | -2.8948 / -2.1397 | -2.8986 / -2.1424 | -0.131 / -0.128 % |

Lift-off speed from the cases standstill and 20000 rpm at 20 degC and at 120 degC, element by
element (methods.md section 3.1), against `rings_liftoff_speed`:

| temperature | FE (rpm) | closed form (rpm) | error |
|---|---|---|---|
| 20 degC (= fitting) | 70 739.0 | 70 739.1 | -2.5e-6 |
| 120 degC | 88 797.7 | 88 798.0 | -2.3e-6 |

The layered path agrees with the closed form as closely as the single ring does (section 1):
the hoop stress jumps from 445 to 177 MPa at the first interface and back to 447 MPa at the
second, and each layer's own values are reproduced. The radial stresses of the two sides of an
interface agree within 0.05 % and the interlaminar r-theta shear is 0.20 MPa (0.05 % of the hoop
stress), the discretisation noise of an axisymmetric problem. The helical layer governs the
Tsai-Wu safety factor (2.68, against 4.44 and 4.47 in the hoop layers), although it carries the
lowest hoop stress. Before the bore arc of the first layer was sized like the lamination OD
(methods.md section 4.6), the finer target mesh made the contact pressure ripple between 7.085
and 7.444 MPa around the right mean; matched, it is uniform to 2e-5.

**Single-layer prestress** (`test_ring_prestress_single_layer`): the orthotropic Table II
ring of section 1 with no interference, a 150 MPa winding prestress, 20000 rpm and 100 K hot,
against the N = 1 closed form: sleeve hoop 317.81 / 318.12 MPa (-0.098 %), contact pressure
5.8078 / 5.8144 MPa (-0.114 %), disc 65.997 / 65.990 MPa (+0.011 %).

**Rotor study** (`examples/configs/mech_multilayer_sleeve.toml`): default bridge-retained
rotor, 1.5 mm sleeve, frictionless sleeve contact, bonded magnets, cases standstill, 20000 rpm
and 24000 rpm at 20 degC (the fitting temperature), lift-off from them. Points: (a) today's
sleeve, one hoop layer with the Table II constants, 0.1 mm interference; (a2) one hoop layer of
carbon_sm/epoxy (vf 0.55, Chamis), the material of (b) and (c), 0.1 mm interference; (b) hoop /
+/-20 deg / hoop layers of 0.75 / 0.375 / 0.375 mm, same material and interference; (c) three
0.5 mm hoop layers, no interference, a uniform winding prestress of 169.7 MPa. The prestress of
(c) was set from one trial run at 100 MPa (standstill minimum contact pressure 1.5440 MPa),
scaled linearly to (a)'s 2.6196 MPa; the run at 169.7 MPa gave 2.6200 MPa, so the response is
linear in the prestress to 1e-4. All runs converged with every contact element closed; 25 160
(single layer), 27 719 (b) and 27 147 (c) nodes, 34-44 s per point for the three cases.
Stresses in MPa, at 20000 rpm unless stated:

| point | `svm_lam_max` (`qhole`), 20000 / 24000 rpm | sleeve hoop max per layer, bore to OD | contact p min, standstill / 20000 / 24000 rpm | Tsai-Wu safety factor (governing layer), 20000 / 24000 rpm | interlaminar shear | lift-off (rpm) |
|---|---|---|---|---|---|---|
| (a) Table II hoop, 0.1 mm | 253.3 / 374.4 | 223.9 | 2.620 / 2.576 / 2.234 | 8.92 / 8.08 | | 39 357 |
| (a2) carbon_sm hoop, 0.1 mm | 251.8 / 372.5 | 228.5 | 2.745 / 2.812 / 2.636 | 8.72 / 7.92 | | 45 711 |
| (b) hoop / +/-20 / hoop, 0.1 mm | 254.2 / 375.2 | 229.0 / 141.5 / 221.8 | 2.521 / 2.529 / 2.223 | 6.48 (2) / 5.93 (2) | 1.37 | 41 589 |
| (c) 3 x hoop, 169.7 MPa prestress, no interference | 252.9 / 373.6 | 218.4 / 209.8 / 214.4 | 2.620 / 2.685 / 2.463 | 9.13 (1) / 8.28 (1) | 1.42 | 44 629 |

* **Lamination.** The headline stress hardly depends on how the sleeve is made: 251.8 to
  254.2 MPa at 20000 rpm, 372.5 to 375.2 MPa at overspeed, always at the lightening-hole
  fillet. What retains the rotor is the sleeve's radial pressure, and all four sleeves give
  2.2 to 2.8 MPa.
* **Helical layer (b against a2, same material).** Replacing a quarter of the hoop thickness by
  a +/-20 deg layer (E_hoop 83 against 128 GPa) lowers the minimum contact pressure by 8-16 %
  and the lift-off speed by 9 % (45 711 to 41 589 rpm). The helical layer carries 141.5 MPa
  against 229 / 222 MPa in the hoop layers beside it, but it governs the safety factor
  (6.48 against 8.73 and 9.02 in the hoop layers), because its fibres are not aligned with the
  load; the sleeve's safety factor falls from 8.72 to 6.48.
* **Winding prestress instead of interference (c).** At the same standstill pressure as (a),
  the prestressed sleeve keeps more pressure at speed (2.685 against 2.576 MPa at 20000 rpm)
  and lifts off later (44 629 against 39 357 rpm). Part of that is the material (compare a2,
  45 711 rpm with 2.745 MPa at standstill); against a2 the prestressed sleeve is close to the
  interference-fitted one, with lower hoop stress (218.4 against 228.5 MPa) and a higher safety
  factor (9.13 against 8.72). A winding prestress sigma0 acts like an interference of about
  sigma0 r / E_hoop: 169.7 MPa x 71.1 mm / 128 GPa = 0.094 mm, against the 0.1 mm of a2.
  The three layers relax only a little (hoop maxima 170.3 / 164.6 / 168.5 MPa at standstill,
  from 169.7 MPa), because the rotor beneath is stiff.
* **Interlaminar shear.** 1.4 MPa in both layered sleeves at 20000 rpm (1.7-1.8 MPa at
  overspeed), from the non-uniform support of the lamination; compare the 70 MPa in-plane shear
  strength of the ply. The radial stress stays compressive in all sleeves except for traces
  (at most 0.03 MPa).
* The simultaneous prestress of (c) is not the build-up of real winding (methods.md section
  4.6): wound layer by layer onto this rotor, the inner layers would be relaxed by the outer
  ones. For a winding stress of 100 MPa on the closed-form model of this rotor
  (`sipm winding --tension 100 --set sleeve_fibre=carbon_sm sleeve_layers=0:1,0:1,0:1`) the
  equivalent prestresses are within 2 % of the winding stress, because the solid disc beneath
  is stiff; the slotted rotor gives way more, so the real relaxation is larger.

**Bridgeless rotor** (`test_mech_point_multilayer_bridgeless`): the default 0.5 mm sleeve as
hoop / +/-20 deg / hoop layers of 0.25 / 0.125 / 0.125 mm (carbon_sm/epoxy, vf 0.55), 0.1 mm
interference, 50 MPa winding prestress, 20000 rpm, converged with every contact element closed,
29 464 nodes, 26 s: `svm_lam_max` 316.1 MPa (`qhole`), sleeve hoop max per layer 286.2 / 197.7 /
291.9 MPa, contact pressure min 1.343 MPa, Tsai-Wu safety factor 4.63 in the helical layer
(6.98 / 6.85 in the hoop layers), interlaminar shear 3.55 MPa. Without the prestress: 318.5 MPa,
233.3 / 148.4 / 244.0 MPa, 0.818 MPa, 6.16 (layer 2), 3.26 MPa; one carbon_sm hoop layer for
comparison: 317.5 MPa, sleeve hoop 241.1 MPa, 1.122 MPa, 8.29. The interlaminar shear is larger
than on the bridge-retained rotor (1.4 MPa), from the edges of the three loose lamination pieces.

## 13. Sleeve eddy-current loss (2D, axial current)

Default bridge-retained rotor (1.5 mm sleeve), 420 A rms, gamma = 130 deg, 3000 rpm,
`em_step` = 2 deg unless stated. The sleeve's axial conductivity is set through
`sleeve_sigma_transverse` with a hoop winding (axial conductivity = sigma_T): 1e4 S/m at the
check point, 1e5 and 1e6 S/m for the scaling. These are test values that span the range, not
material data. What the 2D loss means and leaves out: methods.md section 2.1. Each conducting
run simulates 180 electrical degrees (lead-in plus two loss periods), 123-165 s per run at
`em_step` = 2 deg including the model build (254 s at 1 deg, 449 s at 0.5 deg), against 69.5 s
for the torque run. The field exports of (c) took 6-8 min per run (17 min with 721 points per
arc), mostly in the line reports.

**(a) Non-conducting path unchanged.** Default point after the change: 334.3668 N m, ripple
0.1221, 69.5 s, as in section 6 (334.37 N m, 0.122). The default model, its cache keys
(`tests/test_sleeve_eddy.py` pins `em_model_key` against the previous code, the point keys were
already pinned) and its outputs are those of earlier versions.

**Setting eddy effects.** The first conducting run gave `SolidLoss` = 0 and the torque of the
non-conducting model to 1e-7: PyAEDT's `eddy_effects_on` lists only materials of at least
1e5 S/m as conductors and silently left the 1e4 S/m sleeve out of the eddy-effect table (read
back from the saved project). The driver sets it directly (`SetEddyEffect`) and reads it back
(`GetEddyEffect`).

**Net current.** Three versions of the conducting sleeve at 1e4 S/m (`em_step` 2 deg,
120-degree runs of the same model):

| sleeve object | loss (W) | torque (N m) |
|---|---|---|
| eddy effects on, no winding (the driver's model) | 5.05602 | 334.35253 |
| solid voltage winding, 0 V | 5.05602 | 334.35253 |
| solid current winding, 0 A (zero net current in the modelled pole) | 4.58187 | 334.35388 |

A passive eddy object in Maxwell 2D carries a free net current, exactly as if shorted by a
resistance-free return (0 V). For the ring sleeve that is right: the odd periodicity makes the
current anti-periodic over a pole pitch, so the ring total is zero whatever one pole carries.
Forcing zero net current per pole removes 9.4 % of the loss. The thin-shell estimate (c)
reproduces both versions. On the conducting result itself, `Jz` exported on the mid-radius arc
over the run is anti-periodic across the pole (mismatch 8e-9 of the peak), and the net current
of the modelled pole is not zero: up to 1.30 A (rms 0.55 A) against a mean 1.57 A of the
integral of |Jz| over the pole cross-section, alternating from pole to pole, so the ring total
is zero. (Arcs on the sleeve bore and OD are unusable for `Jz`: they lie on the material
interface, where the element current density jumps.)

**Start-up and averaging.** The loss waveform repeats with the 60-degree period from the first
time step on; only Maxwell's t = 0 sample (a static solution without eddy current) differs. The
change of the period average from the second to the third period (`sleeve_loss_rel_change`) is
7e-6 to 1.1e-5 in every run below.

**Mesh and time step** (1e4 S/m):

| sleeve mesh (mm) | `em_step` (deg) | steps per slot period | loss (W) | torque (N m) |
|---|---|---|---|---|
| 1.2 (`mesh_rotor`, default) | 2 | 15 | 5.0561 | 334.3526 |
| 0.3 | 2 | 15 | 5.0577 | 334.1282 |
| 1.2 | 1 | 30 | 5.1064 | 334.3523 |
| 1.2 | 0.5 | 60 | 5.1199 | 334.3522 |

The sleeve mesh hardly matters for the loss (0.03 %); refining it moves the torque by 0.07 %
through the changed mesh alone, so a conducting model keeps the non-conducting mesh and differs
from it only by the conductivity. The time step matters more: the loss rises by 0.99 % from 2 to
1 deg and by 0.26 % from 1 to 0.5 deg, a ratio of 3.7 (about second order), which extrapolates
to 5.125 W. The default `em_step` of 2 deg (15 steps per slot-passing period, independent of
speed) is therefore about 1.3 % low, 1 deg about 0.4 %. Runs with fewer than 15 steps per slot
period are refused; use `em_step` = 1 for sub-percent loss. The skin depth at the slot-passing
frequency (2400 Hz) is 103 mm at 1e4 S/m and 10.3 mm at 1e6 S/m, 7 to 68 times the sleeve
thickness and far above the 1.2 mm mesh.

**(b) Scaling.** Loss against conductivity at 3000 rpm and against speed at 1e4 S/m, fixed
current and current angle:

| sigma axial (S/m) | speed (rpm) | loss (W) | ratio to previous row | exponent | loss / shaft power | torque (N m) |
|---|---|---|---|---|---|---|
| 0 | 3000 | - | | | | 334.3668 |
| 1e4 | 3000 | 5.0561 | | | 4.8e-5 | 334.3526 |
| 1e5 | 3000 | 49.891 | 9.868 | 0.994 | 4.8e-4 | 334.2273 |
| 1e6 | 3000 | 432.39 | 8.667 | 0.938 | 4.1e-3 | 333.1687 |
| 1e4 | 6000 | 20.199 | 3.995 (to 1e4, 3000 rpm) | 1.998 (speed) | 9.6e-5 | 334.3385 |

The loss is proportional to sigma within 1.3 % up to 1e5 S/m (exponent 0.994) and to the speed
squared (exponent 1.998): the resistance-limited regime. From 1e5 to 1e6 S/m the exponent falls
to 0.938: the sleeve currents begin to screen the harmonics (shell reaction), although the skin
depth is still 7 times the thickness. Above about 1e5 S/m at this speed the thin-shell estimate
(c), which assumes an unscreened field, is no longer valid, while the FE loss is.

**Torque.** The conducting sleeve lowers the torque by 0.0142, 0.1395 and 1.198 N m at 1e4,
1e5 and 1e6 S/m, linear in sigma like the loss and 0.88, 0.88 and 0.87 times the loss over the
mechanical speed (P / Omega = 0.0161, 0.1588, 1.376 N m). The drop need not be P / Omega: with
current-driven windings the electromagnetic power the stator delivers, sum <i dpsi/dt>, changes
too, and over a period sum <i dpsi/dt> = T Omega + P. At 1e5 S/m, with dpsi/dt from the
exported flux linkages (central differences) and the imposed currents (Maxwell's `Current`
report returns 0 for these windings), the stator power rises by 14.6 W and T Omega falls by
43.8 W, against a loss of 49.9 W: the balance closes to 8.5 W. That is within the accuracy of
this voltage-based power, which misses T Omega by 104 W (0.1 %) already for the non-conducting
model, so it supports the reading (part of the loss is supplied electrically) without proving
the 0.88. At no load, where the stator exchanges no power, the drop is P / Omega within 3 %
(below). In the harmonic form of (c) the sleeve's own J x B torque is 0.71 P / Omega at load and
0.92 P / Omega at no load: the eddy field also changes the force on the rotor iron and magnets,
so the sleeve force alone is not the rotor torque change. At realistic carbon conductivities the
torque change is far below the 0.07 % mesh sensitivity of the torque.

**(c) Independent estimate.** From the non-conducting model (same mesh, same time step, fields
saved at every step) the vector potential A_z (`Flux_Lines`) was exported on five arcs through
the sleeve (bore, quarter, middle, three quarters, OD), 0 to 45 deg, over the run
(`MaxwellSession.sleeve_lines`). The modelled pole was extended to the full ring by the odd
periodicity (seam mismatch 4e-8 of the peak), resampled into the rotor frame and the
resistance-limited thin-shell loss evaluated over the last loss period
(`sipm.validation.sleeve_shell_estimate` with `sipm.eddy.shell_loss_from_potential`, methods.md
section 2.1), with the time derivative taken as Maxwell takes it (backward difference) and zero
total current over the ring:

| case | FE loss (W) | estimate (W) | FE / estimate |
|---|---|---|---|
| load, default mesh, zero ring current (the driver's model) | 5.0561 | 5.0923 | 0.993 |
| load, default mesh, zero current per pole (FE: 0 A winding) | 4.5819 | 4.6116 | 0.994 |
| load, sleeve mesh 0.3 mm, zero ring current | 5.0577 | 5.0943 | 0.993 |
| no load, default mesh, zero ring current | 0.06244 | 0.06398 | 0.976 |

The first row is re-run from the package by `sipm validate sleeve-loss` (`--sigma`, default
1e4 S/m; `sipm.validation.sleeve_loss_check`, `test_sleeve_loss_check`): the conducting run, then
the non-conducting model with A_z on the five arcs and the estimate. Re-run: FE 5.0561 W (change
over the last period 8.1e-6), estimate 5.0923 W, FE / estimate 0.9929; variants 4.6115 W (zero
current per pole), 4.9115 W (central differences), 4.8085 W (harmonic form), seam mismatch 1.6e-8.

The FE loss and the estimate agree within 1 % at load and 2.4 % at no load (where the loss is
80 times smaller and made of higher harmonics), for both net-current conditions; the ratio is
the same for both, so the effect of the per-pole constraint is reproduced. The B_r recovered from
the exported A_z peaks at 1.19 T in the sleeve, as the `Mag_B` line report there (1.18 T), which
confirms that `Flux_Lines` is A_z in Wb/m.
Variants of the estimate from the same export (load, default mesh) show its own discretisation
spread: central instead of backward time differences 4.912 W, and the harmonic formula from B_r
at the mid radius (`shell_loss_from_br`, spectral derivatives, field taken as uniform through
the thickness) 4.809 W, 3.5 % and 5.6 % below the backward-difference value. The backward
difference is the one consistent with the solver's time stepping, which is why it is the
comparison above; the variants say that, at this time step, the estimate itself is good to a
few per cent, and the FE step study shows the loss converging upwards with finer steps.
The largest contributions of the harmonic form (rotor frame, load): order 44 at 2400 Hz (the
slot harmonic, 40 % of the loss), order 4 at 2400 Hz (15 %), order 20 at 1200 Hz (14 %), order 12
at 1200 Hz (10 %).

**Load against no load** (1e4 S/m, 3000 rpm): 5.0561 W at 420 A rms, 0.0624 W at zero current.
The no-load loss, from the slotting of the magnet field alone (orders 44 and 52 at 2400 Hz), is
1.2 % of the load loss: in this machine, with 1.9 mm slot openings, the sleeve loss comes from the
armature field, the winding harmonics and the slot harmonics of the stator current. At no load
the stator exchanges no power, so the torque must fall by exactly P / Omega = 1.99e-4 N m; the
FE gives a mean torque of -1.93e-4 N m (the cogging torque of the non-conducting model averages
to 0.000 N m, section 5).

**Bridgeless rotor** (`test_em_point_conducting_sleeve_bridgeless`, 0.5 mm sleeve, 1e4 S/m,
420 A rms, gamma = 130 deg, 3000 rpm): loss 3.352 W (change over the last period 3e-7), torque
497.6629 to 497.6536 N m, a drop of 0.0093 N m = 0.87 P / Omega, as on the bridge-retained
rotor; skin depth 205 times the sleeve thickness. The non-conducting torque is that of
section 6 (497.66 N m).
