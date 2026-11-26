"""Small stress sweeps at 20000 rpm (needs Ansys Mechanical APDL).

    python examples/03_mech_sweeps.py [workdir]
"""

import sys
from pathlib import Path

from sipm.cli import load_config, plot_results, sweep_from_config

here = Path(__file__).parent / "configs"
for name in ("mech_sleeved_bridged", "mech_sleeved_bridgeless"):
    cfg = load_config(here / f"{name}.toml")
    if len(sys.argv) > 1:
        cfg["sweep"]["workdir"] = sys.argv[1]
    sw = sweep_from_config(cfg)
    df = sw.run()
    print(df[["sleeve_thickness", "interference", "status", "svm_rotor_max", "sh_sleeve_max"]].to_string())
    for f in plot_results(df, sw.dir / "figures", "mech", cfg.get("plot")):
        print("figure:", f)
