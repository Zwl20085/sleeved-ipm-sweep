"""Regenerate the README figures in docs/img (no Ansys needed).

    python docs/make_figures.py            # geometry figures
    python docs/make_figures.py --results  # also the result figures from examples/data

* geometry_topologies.png - default parameters of the two sleeved topologies;
* geometry_parametric.png - the sleeved bridge-retained rotor with one parameter group changed;
* torque_vs_gamma.png     - examples/data/sample_em_results.csv (sleeved bridge-retained,
                            1.0 and 2.0 mm sleeve) and, as a reference curve only,
                            sample_em_conventional.csv: a baseline solved with the author's
                            conventional rotor, whose geometry is not part of the package;
* stress_map.png          - examples/data/sample_mech_results.csv, allowable stress 400 MPa
                            (example value, a user input) and interference <= 0.2 x thickness.

The two result figures are drawn from the sample CSV files (outputs of this toolchain with the
default geometry, see examples/data/README.md); they are only redrawn with ``--results``.

Operating conditions in the notes are read from the CSV files, not typed in.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from matplotlib.figure import Figure

from sipm import default_params, post

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "data"
IMG = ROOT / "docs" / "img"

SIGMA_ALLOW = 400.0           # MPa, example lamination strength (user input)
INTERFERENCE_FRACTION = 0.2   # interference <= 0.2 x sleeve thickness

# geometry_parametric.png: panel heading -> parameters changed from the sleeved bridge-retained defaults
PARAMETRIC_PANELS: list[tuple[str, dict[str, float]]] = [
    ("Defaults: 1.5 mm sleeve", {}),
    ("sleeve_thickness = 3 mm: magnets shorten", {"sleeve_thickness": 3.0}),
    ("Outer bridges = 2 mm", {"pm1_bridge_thickness": 2.0, "pm2_bridge_thickness": 2.0}),
    ("pm1_rib_half_width = 1.5, pm2_web_width = 1.2 mm", {"pm1_rib_half_width": 1.5, "pm2_web_width": 1.2}),
]


def torque_figure(out: Path) -> Figure:
    sleeved = pd.read_csv(DATA / "sample_em_results.csv")
    conv = pd.read_csv(DATA / "sample_em_conventional.csv")  # baseline curve (topology label only)
    if set(conv["current_rms"].unique()) != set(sleeved["current_rms"].unique()):
        raise ValueError("the two EM tables were solved at different currents")
    df = pd.concat([conv, sleeved], ignore_index=True)
    # two topologies -> the note keeps only the shared current and speed
    return post.mtpa_lines(df, ["topology", "sleeve_thickness"], note=post.conditions_note(df, "em"), out=out)


def stress_figure(out: Path) -> Figure:
    mech = pd.read_csv(DATA / "sample_mech_results.csv")
    return post.stress_map(mech, SIGMA_ALLOW, INTERFERENCE_FRACTION, note=post.conditions_note(mech, "mech"),
                           out=out)


def parametric_figure(out: Path) -> Figure:
    """The same sleeved bridge-retained rotor with one group of parameters changed per panel."""
    variants = [(heading, default_params("sleeved_bridged", **changes)) for heading, changes in PARAMETRIC_PANELS]
    return post.topology_figure(out=out, params=[p for _, p in variants], headings=[h for h, _ in variants],
                                callouts=False)


def main(argv: list[str] | None = None) -> None:
    import sys

    argv = sys.argv[1:] if argv is None else argv
    IMG.mkdir(parents=True, exist_ok=True)
    post.topology_figure(out=IMG / "geometry_topologies.png")
    parametric_figure(IMG / "geometry_parametric.png")
    from make_construction_figures import main as construction_figures

    construction_figures()
    if "--results" in argv:
        torque_figure(IMG / "torque_vs_gamma.png")
        stress_figure(IMG / "stress_map.png")


if __name__ == "__main__":
    main()
