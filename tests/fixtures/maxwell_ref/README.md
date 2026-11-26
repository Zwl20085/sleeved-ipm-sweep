# Maxwell reference fixtures

Exact outlines of the rotor regions of the author's original Maxwell 2D models of the two
sleeved rotors, and of their stator, used by `tests/test_maxwell_reference.py` to check that
the package geometry reproduces them.

| File | Maxwell design | State |
|---|---|---|
| `stator_saved.json` | stator of the models | one stator pole pitch and its six coil regions |
| `sleeved_bridged_saved.json`, `_v1` .. `_v7` | `01_ICEM_01_simplerotor` | saved state (Wcf = 2 mm) and 7 variable sets |
| `sleeved_bridgeless_saved.json`, `_v1` .. `_v6` | `02_ICEM_01_teslamotor1` | saved state (Wcf = 0.5 mm) and 6 variable sets |

The variable sets span the original sweep ranges (`Wcf` 0.5-2, `RD_WO_lay*` 0.5-1.5,
`RD_FB_lay*` 0.5-1.5, `w_tip` 0.5-2 mm); `index.json` lists them with the Maxwell
overrides and the derived magnet lengths.

**Origin.** The `.aedt` project files were read as text; every rotor object's construction
history (polylines, rectangles, equation-based curves, unite/subtract, mirror/rotate/move,
object coordinate systems, MoveEdges, fillets) was replayed with the variable expressions
evaluated for each variable set (`tools/maxwell_ref`). At the saved state the replay agrees
with every vertex and edge position AEDT caches in the file (max. deviation 1.5e-13 mm).
No Ansys software was run. No result of the Maxwell models is stored, only geometry.

**Content of one file.** `topology`; `design`; `maxwell_overrides` (Maxwell expressions);
`params` (the package parameters that reproduce the state on top of the topology defaults,
except the two deliberate bridgeless deviations of docs/methods.md 1.4: rib cut and layer-1
OD openings);
`derived` (`RD_Lpm_lay1/2`, mm); `regions`: `lamination` (faces: `area`, `outer`, `holes`),
`magnets` (four faces), `sleeve`. Pole frame: shaft axis at the origin, d-axis on +x, one pole
from -22.5 to +22.5 deg, mm. Segments:

* `["L", x0, y0, x1, y1]` line;
* `["A", x0, y0, x1, y1, cx, cy, r, sweep_deg]` circular arc (counter-clockwise positive);
* `["C", x0, y0, x1, y1, [[x, y] x 9]]` equation-based curve with nine points at uniform curve
  parameter from start to end.

Outer loops run counter-clockwise, holes clockwise. Coordinates are rounded to 1e-10 mm.

The stator file holds `topology` (a sleeved topology, for the package defaults), `design`,
`params` (empty), `stator` (one face) and `coils` (six faces). The stator is the same in both
sleeved designs, where it is stored as a body without construction history; its outline was
taken from a Maxwell design that keeps the stator's history, replayed the same way. That
design is not part of this package, so **the stator file cannot be regenerated**:
`tools/make_maxwell_fixtures.py` writes only the rotor files and keeps `stator_saved.json` as
it is. Treat it as a frozen reference.

**Regenerate** the rotor files (needs the `.aedt` project and `shapely`, no licence):

    python tools/make_maxwell_fixtures.py --icem PATH/TO/ICEM_project.aedt --out tests/fixtures/maxwell_ref

The run takes about 30 minutes (one topology per process can run in parallel with `--only`).
