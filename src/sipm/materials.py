"""Electromagnetic material data.

No vendor data are shipped. The default lamination curve is a *generic*
non-oriented silicon-steel BH curve built from a two-term tanh fit
(J = 1.30 T tanh(H/180) + 0.65 T tanh(H/15000), B = mu0 H + J): about 1.41 T at
2.5 kA/m, 1.69 T at 10 kA/m and 1.97 T at 30 kA/m. It is NOT the 27SW1400 grade
named in the paper. Supply your own curve with ``lam_bh = "path/to/bh.csv"``
(two columns H [A/m], B [T]) or use an AEDT library material with
``lam_bh = "aedt:<name>"``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

MU0 = 4e-7 * math.pi


def generic_bh(n: int = 40, h_max: float = 2e5) -> np.ndarray:
    """Generic NO-steel BH table, columns H [A/m], B [T], starting at (0, 0)."""
    H = np.concatenate([[0.0], np.geomspace(20.0, h_max, n - 1)])
    J = 1.30 * np.tanh(H / 180.0) + 0.65 * np.tanh(H / 15000.0)
    return np.column_stack([H, MU0 * H + J])


def load_bh(spec: str) -> np.ndarray | None:
    """BH table for ``Params.lam_bh``; None when an AEDT library material is requested."""
    if spec == "generic":
        return generic_bh()
    if spec.startswith("aedt:"):
        return None
    path = Path(spec)
    if not path.exists():
        raise FileNotFoundError(f"BH curve file not found: {spec}")
    data = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    if data.shape[1] != 2:
        raise ValueError("BH file must have two columns: H [A/m], B [T]")
    if not (np.all(np.diff(data[:, 0]) > 0) and np.all(np.diff(data[:, 1]) > 0)):
        raise ValueError("BH curve must be strictly increasing in H and B")
    if data[0, 0] != 0 or data[0, 1] != 0:
        data = np.vstack([[0.0, 0.0], data])
    return data


def magnet_coercivity(br: float, mur: float) -> float:
    """Coercivity [A/m] of a linear magnet: Hc = Br / (mu0 mur)."""
    return br / (MU0 * mur)
