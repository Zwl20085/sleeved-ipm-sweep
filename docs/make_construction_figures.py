"""Draw the construction of one rotor pole per topology (docs/img/construction_<topology>.png).

    python docs/make_construction_figures.py

Each figure shows one pole, zoomed to x 35..75 mm, y -6..28 mm (lamination grey, holes white,
magnets red, sleeve teal) and marks the independent variables of the Maxwell construction as
dimension arrows, labelled with the package parameter name and, in parentheses, the Maxwell
variable. Layer 1 is drawn in colour C0, layer 2 in C2. The construction points are those of
``rp.layers[k - 1].P`` (see docs/methods.md, section 1).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Arc as MplArc  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

from sipm import default_params  # noqa: E402
from sipm.geometry import RotorPole, build_rotor  # noqa: E402
from sipm.geometry.construction import LayerPoints  # noqa: E402
from sipm.geometry.primitives import chain_polyline  # noqa: E402
from sipm.geometry.rotor import opening_angles  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
IMG = ROOT / "docs" / "img"

XLIM = (35.0, 75.0)
YLIM = (-6.0, 28.0)
LAYER_COLOUR = {1: "C0", 2: "C2"}
FS = 6.6  # label font size
ROW0 = 0.285  # figure fraction: top of the legend table
ROW_DY = 0.0205
Pt = tuple[float, float]

def _fill_face(ax: Axes, outer, holes, face: str) -> None:
    ax.add_patch(Polygon(chain_polyline(outer, 0.5, 40), closed=True, fc=face, ec="0.35", lw=0.5, zorder=1))
    for h in holes:
        ax.add_patch(Polygon(chain_polyline(h, 0.5, 40), closed=True, fc="white", ec="0.35", lw=0.5, zorder=2))


def _draw_bodies(ax: Axes, rp: RotorPole) -> None:
    for f in rp.lamination:
        _fill_face(ax, f.outer, f.holes, "0.78")
    for m in rp.magnets:
        ax.add_patch(Polygon(chain_polyline(m.chain, 0.5, 40), closed=True, fc="#d9534f", ec="0.2", lw=0.5,
                             alpha=0.85, zorder=3))
    if rp.sleeve is not None:
        ro, ri = rp.r_sleeve_out, rp.r_lam
        th = np.linspace(-rp.half_pitch, rp.half_pitch, 200)
        xy = np.vstack([np.c_[ro * np.cos(th), ro * np.sin(th)], np.c_[ri * np.cos(th[::-1]), ri * np.sin(th[::-1])]])
        ax.add_patch(Polygon(xy, closed=True, fc="#2a9d8f", ec="0.2", lw=0.5, alpha=0.8, zorder=3))
    ax.plot(XLIM, (0, 0), ls="-.", color="0.3", lw=0.6, zorder=4)
    ax.text(XLIM[1] - 0.4, 0.35, "d-axis", fontsize=FS, ha="right", va="bottom", color="0.3", zorder=6)
    hp = rp.half_pitch
    r = np.array([rp.r_shaft, rp.r_lam + 3.0])
    ax.plot(r * math.cos(hp), r * math.sin(hp), ls="-.", color="0.3", lw=0.6, zorder=4)
    ax.text(39.0, 27.3, "q-axis", fontsize=FS, ha="left", va="top", color="0.3", rotation=22.5, zorder=6,
            rotation_mode="anchor")


class Legend:
    """Numbered tags drawn next to the features; the legend text is collected for the table below."""

    def __init__(self, ax: Axes) -> None:
        self.ax = ax
        self.rows: dict[int, list[tuple[str, str, str]]] = {1: [], 2: [], 0: []}

    def tag(self, group: int, text: str, at: Pt, colour: str) -> None:
        rows = self.rows[group]
        name = f"{group if group else 'S'}{'abcdefghij'[len(rows)]}"
        rows.append((name, text, colour))
        self.ax.text(at[0], at[1], name, fontsize=FS + 0.6, color=colour, ha="center", va="center", zorder=9,
                     fontweight="bold", bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=colour, lw=0.5))


def _dim(ax: Axes, a: Pt, b: Pt, colour: str, style: str = "<->", ls: str = "-", lw: float = 0.9) -> None:
    ax.annotate("", xy=b, xytext=a, zorder=7,
                arrowprops=dict(arrowstyle=style, color=colour, lw=lw, ls=ls, shrinkA=0, shrinkB=0,
                                mutation_scale=7))


def _mid(a: Pt, b: Pt) -> Pt:
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)


def _mark_layer(lg: Legend, rp: RotorPole, L: LayerPoints) -> None:
    ax, p, k = lg.ax, rp.p, L.inp.k
    c = LAYER_COLOUR[k]
    P = L.P
    th = L.inp.theta
    u = np.array(L.u)
    n_pole = np.array([math.sin(th), -math.cos(th)])  # normal to the axis, towards the pole side
    x0 = (L.inp.xpm, 0.0)
    ax.plot(*x0, "o", color=c, ms=3.5, zorder=8)
    ax.plot([x0[0], x0[0] + 22 * u[0]], [0, 22 * u[1]], ls=":", color=c, lw=0.8, zorder=5)
    lg.tag(k, f"pm{k}_d_crossing (RD_xpm_lay{k}) = {L.inp.xpm:.2f} mm", (x0[0], -1.7), c)
    rad = 7.0 if k == 1 else 6.0
    ax.add_patch(MplArc(x0, 2 * rad, 2 * rad, theta1=0, theta2=math.degrees(th), color=c, lw=0.9, zorder=7))
    mid_a = th / 2
    lg.tag(k, f"pm{k}_angle / 2 (RD_thetapm_lay{k}) = {math.degrees(th):.2f} deg",
           (x0[0] + (rad + 1.6) * math.cos(mid_a), (rad + 1.6) * math.sin(mid_a)), c)
    _dim(ax, P["P8"], P["P9"], c)
    lg.tag(k, f"pm{k}_width (RD_Wpm_lay{k}) = {L.inp.W:.2f} mm",
           tuple(np.array(_mid(P["P8"], P["P9"])) - 1.9 * u), c)
    _dim(ax, P["P5"], P["P8"], "k", lw=0.8)
    lg.tag(k, f"pm{k}_end_gap (RD_dpm_lay{k}) = {L.inp.dpm:.2f} mm",
           tuple(np.array(_mid(P["P5"], P["P8"])) + 1.5 * n_pole), c)
    off = 1.6 * n_pole
    a, b = np.array(P["P8"]) + off, np.array(P["P11"]) + off
    _dim(ax, tuple(a), tuple(b), c, ls="--")
    lg.tag(k, f"Lpm = {L.length:.2f} mm, derived (K_Lpm x free slot length), dashed",
           tuple((a + b) / 2 + 1.5 * n_pole), c)
    P3 = np.array(P["P3"])
    out = P3 / float(np.hypot(*P3)) * rp.r_lam
    _dim(ax, tuple(P3), tuple(out), "k", style="<|-|>", lw=0.6)
    sub = ", sets Lpm only" if p.is_bridgeless else ""
    lg.tag(k, f"pm{k}_bridge_thickness (RD_WO_lay{k}) = {L.inp.WO:.2f} mm{sub}",
           tuple(out / np.hypot(*out) * (rp.r_sleeve_out + 1.7)), c)
    P5 = P["P5"]
    if k == 1:
        _dim(ax, (P5[0], 0.0), P5, "k", style="<|-|>", lw=0.6)
        lg.tag(k, f"pm1_rib_half_width (RD_FB_lay1) = {L.inp.FB:.2f} mm", (P5[0] + 0.3, -1.7), c)
    else:
        FBp = P["FB"]
        _dim(ax, P5, FBp, "k", style="<|-|>", lw=0.6)
        what = "rib half width here" if p.is_bridgeless else "web beside the cavity"
        lg.tag(k, f"pm2_web_width (RD_FB_lay2) = {L.inp.FB:.2f} mm ({what})", (FBp[0] + 0.3, -1.7), c)
    s0 = getattr(p, f"pm{k}_slot_start")
    end = (x0[0] + s0 * u[0], s0 * u[1])
    _dim(ax, x0, end, "k", style="<|-|>", lw=0.9)
    lg.tag(k, f"pm{k}_slot_start = {s0:g} mm (Maxwell strip start, along the axis)",
           tuple(np.array(_mid(x0, end)) - 1.6 * n_pole), c)


def _mark_sleeve_and_lips(lg: Legend, rp: RotorPole) -> None:
    ax, p = lg.ax, rp.p
    if rp.sleeve is not None:
        a = math.radians(4.0)
        t = np.array((math.cos(a), math.sin(a)))
        _dim(ax, tuple(rp.r_lam * t), tuple(rp.r_sleeve_out * t), "k", style="<|-|>", lw=0.6)
        lg.tag(0, f"sleeve_thickness (Wcf) = {p.sleeve_thickness:g} mm, radial", tuple((rp.r_sleeve_out + 1.7) * t),
               "#1f6f66")
    if p.is_bridgeless:
        a = math.radians(15.5)
        t = np.array((math.cos(a), math.sin(a)))
        _dim(ax, tuple((rp.r_lam - p.lip_depth) * t), tuple(rp.r_lam * t), "k", style="<|-|>", lw=0.6)
        lg.tag(0, f"lip_depth (w_tip) = {p.lip_depth:g} mm, radial", tuple((rp.r_sleeve_out + 1.7) * t), "0.2")
        for L in rp.layers:
            k = L.inp.k
            lo, hi = opening_angles(p, L)
            ro = rp.r_sleeve_out + 0.6
            ang = np.linspace(lo, hi, 20)
            ax.plot(ro * np.cos(ang), ro * np.sin(ang), color=LAYER_COLOUR[k], lw=1.6, zorder=8)
            for a_ in (lo, hi):
                ax.plot([rp.r_lam * math.cos(a_), (ro + 0.4) * math.cos(a_)],
                        [rp.r_lam * math.sin(a_), (ro + 0.4) * math.sin(a_)], color=LAYER_COLOUR[k], lw=0.6,
                        zorder=8)
            mid = 0.5 * (lo + hi)
            lg.tag(0, f"pm{k}_opening_angle ({'rib2_open' if k == 1 else 'rib1_open'}) = "
                      f"{getattr(p, f'pm{k}_opening_angle'):g} deg, opening "
                      f"{math.degrees(lo):.2f}..{math.degrees(hi):.2f} deg (q-side edge = reference - "
                      f"pm{k}_opening_shift {getattr(p, f'pm{k}_opening_shift'):g} deg)",
                   ((ro + 1.5) * math.cos(mid), (ro + 1.5) * math.sin(mid)), LAYER_COLOUR[k])


def _legend_table(fig: Figure, lg: Legend) -> None:
    """Write the tag texts below the axes: layer 1 left, layer 2 right, others at the bottom."""
    for group, x in ((1, 0.01), (2, 0.52)):
        y = ROW0
        for name, text, colour in lg.rows[group]:
            fig.text(x, y, f"{name}  {text}", fontsize=FS, color=colour, ha="left", va="top")
            y -= ROW_DY
    y = ROW0 - ROW_DY * max(len(lg.rows[1]), len(lg.rows[2])) - 0.006
    for name, text, colour in lg.rows[0]:
        fig.text(0.01, y, f"{name}  {text}", fontsize=FS, color=colour, ha="left", va="top")
        y -= ROW_DY


def construction_figure(topology: str, out: Path) -> Figure:
    kw = {"sleeve_thickness": 2.0} if topology == "sleeved_bridged" else {}
    p = default_params(topology, **kw)
    rp = build_rotor(p)
    fig = plt.figure(figsize=(7.0, 8.4))
    ax = fig.add_axes((0.08, 0.32, 0.90, 0.65))
    _draw_bodies(ax, rp)
    lg = Legend(ax)
    for L in rp.layers:
        _mark_layer(lg, rp, L)
    _mark_sleeve_and_lips(lg, rp)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_aspect("equal", anchor="N")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"{topology}, sleeve {p.sleeve_thickness:g} mm: independent variables of the construction", fontsize=9)
    _legend_table(fig, lg)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return fig


def main() -> None:
    IMG.mkdir(parents=True, exist_ok=True)
    for topology in ("sleeved_bridged", "sleeved_bridgeless"):
        construction_figure(topology, IMG / f"construction_{topology}.png")
        print(f"figure: {IMG / f'construction_{topology}.png'}")


if __name__ == "__main__":
    main()
