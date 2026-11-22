"""Validation runs reported in docs/validation.md (needs MAPDL and AEDT).

    python examples/05_validation.py [workdir]
"""

import json
import sys

from sipm.config import Config
from sipm.maxwell import MaxwellSession
from sipm.params import default_params
from sipm.validation import RING_CASES, mesh_study, open_circuit_check, ring_check

cfg = Config(workdir=sys.argv[1]) if len(sys.argv) > 1 else Config()

for case in RING_CASES:  # isotropic, orthotropic, thermal (orthotropic, rotor 100 K above fitting)
    r = ring_check(case, cfg=cfg)
    for row in r["rows"]:
        print(f"ring {case:11s} {row['quantity']:<40s} FE {row['fe']:10.4f}  closed form {row['closed_form']:10.4f}"
              f"  error {100 * row['rel_error']:+.3f} %")

for topo in ("sleeved_bridged", "sleeved_bridgeless"):
    for row in mesh_study(default_params(topo), cfg=cfg):
        print("mesh", topo, row)

p = default_params("sleeved_bridged")
with MaxwellSession(cfg) as s:
    print(json.dumps(open_circuit_check(s.open_circuit(p), p), indent=1, default=float))
    for g in (90, 100, 110, 120, 125, 130, 135, 140, 150, 160, 170):
        r = s.solve(p.replace(gamma=float(g)))
        print(f"gamma {g:5.1f} deg  torque {r['torque_avg']:8.2f} N m")
