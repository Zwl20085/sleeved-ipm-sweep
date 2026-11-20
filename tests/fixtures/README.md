Test fixtures

* `mapdl_results_bridged.txt` - results file written by the generated APDL deck for the
  default `sleeved_bridged` design (MAPDL 2024 R1). Used to test the parser only.
* `mapdl_error_excerpt.out` - excerpt in the layout of a MAPDL output file with one ERROR
  block (text of a real rigid-body error message), used to test error extraction.
* `aedt_report.csv` - first rows of a real AEDT 2024 R1 report export (time in us, torque in
  NewtonMeter), used to test unit handling of the report parser.
