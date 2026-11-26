"""Construction points of one magnet layer (the ``RD_*`` construction of the Maxwell models).

Pole frame: d-axis on +x, the upper magnet of the V (y > 0) is constructed, the lower one
is its mirror image. For layer k (1 = outer, 2 = inner) the independent quantities are

    xpm    d-axis crossing of the magnet centre line         (Maxwell RD_xpm_layk)
    theta  angle of the magnet axis to the d-axis = V/2       (RD_thetapm_layk)
    W      slot width = nominal magnet thickness              (RD_Wpm_layk)
    dpm    gap along the pole-side wall, slot corner -> magnet (RD_dpm_layk)
    K      magnet length / free slot length                   (RD_K_Lpm_layk)
    FB     rib / web measure (layer 1: half centre rib)        (RD_FB_layk)
    INNO   extra height of the slot's inner end                (RD_INNO_layk)
    tINNO  tilt of the slot's inner end                        (RD_thetaINNO_layk)
    WO     outer bridge thickness below the lamination OD      (RD_WO_layk)

The two slot walls are the lines 1 (pole side, d-axis intercept x1) and 2 (q side, x2),
parallel to the magnet axis and W apart. Point names follow the Maxwell variables
``RD_x<n>_layk`` / ``RD_y<n>_layk``; every expression below is the Maxwell expression in mm.
The magnet length is derived: Lpm = K (|P3 - P5| - dpm), P3 = line 1 on the circle r_lam - WO.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .primitives import GeometryError, Pt


@dataclass(frozen=True)
class LayerInput:
    k: int
    xpm: float
    theta: float        # rad
    W: float
    dpm: float
    K: float
    FB: float
    INNO: float
    tINNO: float        # rad
    WO: float


@dataclass(frozen=True)
class LayerPoints:
    inp: LayerInput
    r_lam: float
    P: dict[str, Pt] = field(default_factory=dict)
    scalars: dict[str, float] = field(default_factory=dict)

    @property
    def length(self) -> float:
        """Derived magnet length Lpm (mm)."""
        return self.scalars["Lpm"]

    @property
    def u(self) -> Pt:
        t = self.inp.theta
        return (math.cos(t), math.sin(t))

    def p(self, name: str) -> Pt:
        return self.P[name]


def _line_circle_x(x0: float, m: float, R: float) -> float:
    """Larger root x of (x, m (x - x0)) on the circle of radius R (Maxwell RD_x3 / RD_x4)."""
    A = 1 + m * m
    B = -2 * m * m * x0
    C = -R * R + m * m * x0 * x0
    disc = B * B - 4 * A * C
    if disc < 0:
        raise GeometryError(f"slot wall (d-axis intercept {x0:.3f} mm) misses the bridge circle r = {R:.3f} mm")
    return (-B + math.sqrt(disc)) / (2 * A)


def layer_points(inp: LayerInput, r_lam: float) -> LayerPoints:
    """Evaluate the Maxwell construction of one layer for a lamination radius ``r_lam``."""
    k = inp.k
    th, W, dpm, FB, INNO, ti = inp.theta, inp.W, inp.dpm, inp.FB, inp.INNO, inp.tINNO
    c, s = math.cos(th), math.sin(th)
    if not 0 < th < math.pi / 2:
        raise GeometryError(f"layer {k}: magnet angle must be in (0, 180) deg (included V angle)")
    R = r_lam - inp.WO
    if R <= 0:
        raise GeometryError(f"layer {k}: outer bridge thicker than the lamination radius")
    m = math.tan(th)
    x1 = inp.xpm + W / 2 / s
    x2 = inp.xpm - W / 2 / s
    x3 = _line_circle_x(x1, m, R)
    y3 = m * (x3 - x1)
    x4 = _line_circle_x(x2, m, R)
    y4 = m * (x4 - x2)
    dFB = FB + INNO
    x5, y5 = x1 + dFB / math.tan(th), dFB
    d2 = W / s
    cons1 = (math.tan(th) * math.tan(ti)) / (math.tan(th) + math.tan(ti))
    dY = d2 * cons1
    dX = dY / math.tan(th)
    x6, y6 = x2 + dFB / math.tan(th) + dX, dFB + dY
    x8, y8 = x5 + dpm * c, y5 + dpm * s
    x9, y9 = x8 - W * s, y8 + W * c
    distance3 = math.hypot(x3 - x5, y3 - y5) - dpm
    free = (x3 - x5) * c + (y3 - y5) * s - dpm   # signed: P3 must lie beyond the magnet's inner end
    Lpm = inp.K * distance3
    x10, y10 = x9 + Lpm * c, y9 + Lpm * s
    x11, y11 = x8 + Lpm * c, y8 + Lpm * s
    a_mid = math.atan(y3 / x3) / 2 + math.atan(y4 / x4) / 2
    mid = (R * math.cos(a_mid), R * math.sin(a_mid))
    P: dict[str, Pt] = {
        "P1": (x1, 0.0), "P2": (x2, 0.0), "P3": (x3, y3), "P4": (x4, y4), "P5": (x5, y5), "P6": (x6, y6),
        "P8": (x8, y8), "P9": (x9, y9), "P10": (x10, y10), "P11": (x11, y11), "mid": mid,
    }
    # ---- central cavity control points (layer 2) ----------------------------------------
    P["FB"] = (x5 - FB * math.sin(ti), y5 - FB * math.cos(ti))
    P["FBd"] = (x6 + 1.0, y6 - 1.0 * math.tan(ti))
    scalars = dict(x1=x1, x2=x2, dFB=dFB, d2=d2, distance3=distance3, free_length=free, Lpm=Lpm, R_bridge=R,
                   angle_middle=a_mid)
    # ---- polar angle where the magnet centre line meets the lamination OD -----------------
    # (Maxwell angle_temp_layk, used by the bridgeless rotor to place the OD openings). The
    # centre line (xpm, 0) + t (cos th, sin th) meets r = r_lam at t = -xpm cos th + sqrt(root).
    # The model's expression for layer 1 uses the layer-2 d-axis crossing in the x term (a
    # typo: 13.616 deg instead of 12.394 deg at the saved state); the package uses the true
    # angle for both layers (see docs/methods.md, 1.4).
    root = r_lam ** 2 - (inp.xpm * s) ** 2
    if root < 0:
        raise GeometryError(f"layer {k}: magnet centre line misses the lamination OD")
    t_od = -inp.xpm * c + math.sqrt(root)
    scalars["opening_reference"] = math.atan2(t_od * s, inp.xpm + t_od * c)
    return LayerPoints(inp=inp, r_lam=r_lam, P=P, scalars=scalars)
