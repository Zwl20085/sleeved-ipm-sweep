Test fixtures

* `mapdl_results_bridged.txt` - results file written by the generated APDL deck for the
  default `sleeved_bridged` design (MAPDL 2024 R1). Used to test the parser only.
* `mapdl_sleeve_nodes.txt` - six rows in the layout of the sleeve node file the deck writes
  (r, theta, SX, SY, SZ, SXY in the cylindrical system, Fortran E21.12 as MAPDL 2024 R1 writes
  it); values hand-made (some extremes copied from a real default-point run). Used to test
  `parse_sleeve_nodes` only.
* `mapdl_error_excerpt.out` - excerpt in the layout of a MAPDL output file with one ERROR
  block (text of a real rigid-body error message), used to test error extraction.
* `aedt_report.csv` - first rows of a real AEDT 2024 R1 report export (time in us, torque in
  NewtonMeter), used to test unit handling of the report parser.
