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

## 2. IPM rotor stress, one point per topology (20000 rpm, default mesh)

Topology defaults (saved states of the Maxwell models): sleeved bridge-retained with a 1.5 mm
sleeve, sleeved bridgeless with a 0.5 mm sleeve, 0.1 mm interference, frictionless sleeve
contact, bonded magnets (the designed 0.1 mm clearance bridged by the bond). Integration tests
(`pytest -m ansys`); both converge.

| topology | headline `svm_lam_max` (MPa) | section, location r, theta | global lamination peak (MPa), location | magnet peak (MPa) | max sleeve hoop (MPa) | open sleeve contact | nodes | wall |
|---|---|---|---|---|---|---|---|---|
| sleeved_bridged | 253.3 | `qhole`, 54.42 mm, 22.35 deg | 691.1, 61.15 mm, 2.06 deg | 199.6 | 223.9 | 0 % | 25 160 | 22.7 s |
| sleeved_bridgeless | 318.1 | `qhole`, 54.42 mm, 22.35 deg | 502.6, 60.33 mm, 1.23 deg | 220.9 | 236.6 | 0 % | 21 932 | 20.2 s |

Section peaks (MPa): bridge-retained `bridge1` 103.9, `bridge2` 97.9, `rib1` 13.0, `web2` 76.4,
`qhole` 253.3; bridgeless `lip1_pole` 6.3, `lip1_q` 3.9, `lip2_pole` 3.3, `lip2_q` 3.5, `qhole`
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
| 0.8 / 0.3 | 9 751 | 4.4 | 3.8 | 2.5 | 3.8 | 315.9 | 359.7 | 206.8 | 234.1 | 9 s |
| 0.5 / 0.15 (default) | 21 932 | 6.3 | 3.9 | 3.3 | 3.5 | 318.1 | 502.6 | 220.9 | 236.6 | 18 s |
| 0.3 / 0.08 | 56 409 | 8.2 | 4.6 | 4.1 | 4.4 | 321.0 | 614.7 | 282.3 | 238.4 | 41 s |
| 0.2 / 0.04 | 136 309 | 8.8 | 5.1 | 4.4 | 4.2 | 323.5 | 849.3 | 655.0 | 239.1 | 96 s |
| change between the two finest | | +8.0 % | +10.8 % | +8.5 % | -3.2 % | +0.8 % | +38 % | +132 % | +0.3 % | |

(stresses in MPa)

* **Converged (< 3 % between the two finest meshes):** every section of the bridge-retained
  rotor (0.4-1.8 %), the bridgeless `qhole` (0.8 %), the sleeve hoop stress of both rotors,
  and with them the headline `svm_lam_max` of both rotors (`qhole` at every mesh).
* **Not converged: the four bridgeless lips** (3-11 %). Their peaks lie at the end of the
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
