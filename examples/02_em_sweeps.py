"""Small EM sweeps for the two sleeved topologies (needs Ansys Electronics Desktop).

    python examples/02_em_sweeps.py [workdir]

Each point is one Maxwell 2D transient run over 60 electrical degrees (about a minute
on a desktop PC). Results are cached: running the script again only solves missing points.
"""

import sys
from pathlib import Path

from sipm.cli import load_config, plot_results, sweep_from_config

here = Path(__file__).parent / "configs"
for name in ("em_sleeved_bridged", "em_sleeved_bridgeless"):
    cfg = load_config(here / f"{name}.toml")
    if len(sys.argv) > 1:
        cfg["sweep"]["workdir"] = sys.argv[1]
    sw = sweep_from_config(cfg)
    df = sw.run()
    cols = [c for c in ("gamma", "sleeve_thickness", "lip_depth", "status", "torque_avg") if c in df]
    print(df[cols].to_string())
    for f in plot_results(df, sw.dir / "figures", "em"):
        print("figure:", f)
