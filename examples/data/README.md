# Sample results

Small results tables written by `sipm.sweep` (`results.csv` of the sweep), so that
`examples/04_plots_from_samples.py` and `sipm plot` can be tried without Ansys. They are
outputs of this toolchain, computed with Ansys Electronics Desktop 2024 R1 and Mechanical
APDL 2024 R1 from the configs in `examples/configs`, with the default geometry (the saved
states of the author's Maxwell models, docs/methods.md section 1). They are not results or
data of the companion paper. Wall-time columns (`em_wall_s`, `mech_wall_s`, `wall_s`) were
measured on one workstation while other solver runs shared it, so they only indicate the
order of magnitude.

* `sample_em_results.csv` - `examples/configs/em_sleeved_bridged.toml`: sleeved bridge-retained
  rotor, sleeve thickness 1.0 and 2.0 mm x current advance angle 110-150 deg, 420 A rms,
  3000 rpm, average over 60 electrical degrees.
* `sample_em_bridgeless.csv` - `examples/configs/em_sleeved_bridgeless.toml`: sleeved
  bridgeless rotor (0.5 mm sleeve), lip depth 0.5 and 1.5 mm x current advance angle 115, 125
  and 135 deg, 420 A rms.
* `sample_mech_results.csv` - `examples/configs/mech_sleeved_bridged.toml`: sleeved
  bridge-retained rotor at 20000 rpm, bonded magnets, sleeve thickness 1.0/1.5/2.0 mm x radial
  interference 0-0.4 mm. `svm_rotor_max` (= `svm_lam_max`) is the peak von Mises stress over
  the named lamination sections (`svm_bridge1`, `svm_bridge2`, `svm_rib1`, `svm_web2`,
  `svm_qhole`; `svm_lam_where` names the one that holds the peak). `svm_lam_max_all` and
  `svm_pm_max` are mesh-dependent (docs/validation.md section 3).
* `sample_mech_bridgeless.csv` - `examples/configs/mech_sleeved_bridgeless.toml`: sleeved
  bridgeless rotor at 20000 rpm, bonded magnets, sleeve 1.5/2.0 mm x interference 0/0.2/0.4 mm;
  sections `svm_qhole` and the four lips `svm_lip1_pole` ... `svm_lip2_q` (the lips are not
  mesh-converged).
* `sample_em_conventional.csv` - **reference baseline only**: average torque against current
  advance angle (110-160 deg, 420 A rms, 3000 rpm) computed earlier with this toolchain for
  the author's conventional bridge-retained rotor without a sleeve. That rotor's geometry is
  not part of the package, and the package cannot rebuild it; the table keeps only the
  operating point and torque columns (`topology` = `conventional_baseline`, a display label).
  It is drawn as the reference curve in `docs/img/torque_vs_gamma.png`.
