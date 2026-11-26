"""Composite sleeve constants against fibre volume fraction and winding angle (no Ansys needed).

    python examples/06_composite_sleeve.py [out_dir]

Prints the sleeve constants that sipm derives from a composite spec (carbon_sm fibre class in an
epoxy matrix, Chamis micromechanics, balanced +/- winding) and, with ``out_dir``, saves a figure
of them. The presets are illustrative class-typical numbers, not the data of a commercial grade.
"""

import sys
from pathlib import Path

import numpy as np

from sipm.params import default_params

SPEC = dict(sleeve_fibre="carbon_sm", sleeve_matrix="epoxy", sleeve_micromechanics="chamis")
COLS = ("sleeve_E_hoop", "sleeve_E_axial", "sleeve_E_radial", "sleeve_G_hoop_axial", "sleeve_nu_hoop_axial",
        "sleeve_cte_hoop", "sleeve_cte_axial", "sleeve_density")


def constants(vf: float, angle: float) -> dict:
    p = default_params("sleeved_bridged", sleeve_vf=vf, sleeve_wind_angle=angle, **SPEC)
    return {k: getattr(p, k) for k in COLS}


def table(title: str, xname: str, xs, rows) -> None:
    print(title)
    print(f"{xname:>8s}  E_hoop  E_axial E_radial  G_h-a  nu_h-a  cte_hoop  cte_axial  density")
    print(f"{'':>8s}  (GPa)   (GPa)   (GPa)    (GPa)  (-)     (1e-6/K)  (1e-6/K)   (kg/m3)")
    for x, c in zip(xs, rows):
        print(f"{x:8.3g}  {c['sleeve_E_hoop']:6.1f}  {c['sleeve_E_axial']:6.2f}  {c['sleeve_E_radial']:6.2f}   "
              f"{c['sleeve_G_hoop_axial']:5.2f}  {c['sleeve_nu_hoop_axial']:6.3f}  {1e6 * c['sleeve_cte_hoop']:8.2f}  "
              f"{1e6 * c['sleeve_cte_axial']:9.2f}  {c['sleeve_density']:7.0f}")
    print()


vfs = np.round(np.arange(0.40, 0.701, 0.05), 3)
by_vf = [constants(v, 0.0) for v in vfs]
table("Hoop winding (0 deg), against fibre volume fraction:", "vf", vfs, by_vf)
print("Table II of the paper: E_hoop 125, E_radial 8.8 GPa, nu 0.28, 1800 kg/m3; this preset reaches\n"
      f"E_hoop = 125 GPa at vf = 0.536 ({constants(0.536, 0.0)['sleeve_E_hoop']:.1f} GPa).\n")

angles = np.arange(0.0, 90.1, 7.5)
by_angle = [constants(0.55, a) for a in angles]
table("vf = 0.55, against the winding angle from hoop (balanced +/- winding):", "deg", angles, by_angle)

if len(sys.argv) > 1:
    import matplotlib.pyplot as plt

    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6), constrained_layout=True)
    for k, lab in (("sleeve_E_hoop", "E hoop"), ("sleeve_E_axial", "E axial"), ("sleeve_G_hoop_axial", "G hoop-axial")):
        ax1.plot(angles, [c[k] for c in by_angle], marker="o", ms=3, label=lab)
    ax1.set_xlabel("Winding angle from hoop (deg)")
    ax1.set_ylabel("Modulus (GPa)")
    ax1.legend(frameon=False)
    for k, lab in (("sleeve_cte_hoop", "hoop"), ("sleeve_cte_axial", "axial")):
        ax2.plot(angles, [1e6 * c[k] for c in by_angle], marker="o", ms=3, label=lab)
    ax2.axhline(12.0, color="0.5", ls="--", lw=0.8, label="steel (12e-6 /K)")
    ax2.set_xlabel("Winding angle from hoop (deg)")
    ax2.set_ylabel("Thermal expansion (1e-6 /K)")
    ax2.legend(frameon=False)
    f = out / "composite_constants.png"
    fig.savefig(f, dpi=150)
    print(f"figure: {f}")
