"""Post-processing from the committed sample results (licence free).

    python examples/04_plots_from_samples.py [output_dir]

The CSV files in examples/data were produced by the sweep runner with the default
illustrative geometry (see examples/data/README.md); they are not results of the paper.
"""

import sys
from pathlib import Path

import pandas as pd

from sipm import post
from sipm.sweep import mtpa

data = Path(__file__).parent / "data"
out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("figures")

em = pd.read_csv(data / "sample_em_results.csv")
em_note = post.conditions_note(em, "em")  # e.g. "Sleeved bridge-retained rotor, 420 A rms, 3000 rpm"
post.torque_scatter(em, ["gamma", "sleeve_thickness"], out=out / "sample_torque_scatter.png")
post.mtpa_lines(em, "sleeve_thickness", note=em_note, out=out / "sample_torque_vs_gamma.png")
print(mtpa(em)[["sleeve_thickness", "gamma", "torque_avg"]].to_string(index=False))

mech = pd.read_csv(data / "sample_mech_results.csv")
mech_note = post.conditions_note(mech, "mech")  # e.g. "Sleeved bridge-retained rotor, 20 000 rpm"
post.stress_lines(mech, "svm_rotor_max", note=mech_note, out=out / "sample_rotor_stress.png")
post.stress_lines(mech, "sh_sleeve_max", note=mech_note, out=out / "sample_sleeve_stress.png")
# 400 MPa = example lamination yield strength; replace with your material's value
post.stress_map(mech, sigma_allow=400.0, interference_fraction=0.2, note=mech_note,
                out=out / "sample_stress_map.png")
feas = post.feasible_points(mech, 400.0, 0.2)
print(feas[["sleeve_thickness", "interference", "svm_rotor_max"]].to_string(index=False))
print("figures in", out.resolve())
