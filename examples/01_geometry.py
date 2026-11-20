"""Draw one pole of each topology with the default dimensions (saved states of the author's
Maxwell models) and print the derived magnet lengths. Licence free.

    python examples/01_geometry.py [output_dir]
"""

import sys
from pathlib import Path

from sipm import TOPOLOGIES, build_rotor, default_params
from sipm.post import geometry_plot

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("figures")
for topo in TOPOLOGIES:
    p = default_params(topo)
    m = build_rotor(p).metrics
    print(f"{topo:20s} lamination OD {2 * p.r_lam:.2f} mm, magnet lengths (derived) "
          f"{m['pm1_length']:.3f} / {m['pm2_length']:.3f} mm, q-axis web {m['q_web']:.2f} mm")
    geometry_plot(p, out=out / f"geometry_{topo}.png")
    print("  ->", out / f"geometry_{topo}.png")
