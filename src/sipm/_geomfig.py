"""Cross-section figures of one rotor + stator pole (used by :mod:`sipm.post`).

Every colour has a legend entry, and the features the toolchain sweeps are called out with
leader lines whose end points are computed from the geometry objects (``build_rotor`` /
``build_stator``) so they follow any parameter change. A callout whose feature does not
exist in a topology (e.g. the top bridge of the bridgeless rotor) is skipped.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from . import _style as st

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.legend import Legend

    from .geometry import RotorPole, StatorPole
    from .params import Params

GEOMETRY_COLORS = {
    "stator": "#d9c3a5",
    "rotor": "#9fb7cc",
    "magnet": "#3060c0",
    "conductor": "#d08030",
    "sleeve": "#4d4d4d",
    "air": "#ffffff",
}

PANEL_HEADINGS = {
    "sleeved_bridged": "Sleeved, bridge-retained",
    "sleeved_bridgeless": "Sleeved, bridgeless",
}

CALLOUT_TEXT = {
    "top_bridge": "Outer bridge",
    "centre_rib": "Centre rib",
    "sleeve": "Carbon-fibre sleeve",
    "slot_opening": "OD opening",
    "lip": "Pole-tip lip",
    "air_gap": "Air gap",
    "magnet_length": "Magnet length\n(derived)",
    "cavity": "Central cavity",
    "lightening_hole": "Lightening hole",
    "d_axis": "d-axis",
}
# Corner of the panel that holds the callout text. The two left corners lie outside the
# 360/poles-degree pole sector, so they hold no geometry, for 6 poles and more; at 4 poles
# a long label can reach over the stator back (it keeps its white box and stays readable).
_CORNER = {"top_bridge": "tl", "air_gap": "tl", "slot_opening": "tl", "lip": "tl", "magnet_length": "tl",
           "lightening_hole": "bl", "centre_rib": "bl", "sleeve": "bl", "cavity": "bl"}
# Callouts per panel of topology_figure: each feature once, on the panel where it is typical.
PANEL_CALLOUTS = {
    # top-left labels stack in the order of their anchors (outer bridge, air gap, magnet) so that
    # the leader lines do not cross
    "sleeved_bridged": ("top_bridge", "air_gap", "magnet_length", "centre_rib", "sleeve", "cavity", "d_axis"),
    "sleeved_bridgeless": ("slot_opening", "lip", "lightening_hole", "d_axis"),
}
MAGNET_ARROW_MM = 2.0
CALLOUT_ROW_STEP = 0.1        # axes fraction between stacked callout labels in one corner
CALLOUT_MARGIN = 0.03         # axes fraction from the panel edge to a callout label
LIMIT_PAD_MM = 1.5            # margin around the drawn pole


def _polar(r: float, th: float) -> np.ndarray:
    return r * np.array([math.cos(th), math.sin(th)])


def callout_anchors(rp: RotorPole, sp: StatorPole) -> dict[str, np.ndarray]:
    """Data-coordinate (mm) end point of each callout leader that exists for this geometry.

    * top_bridge      - middle of the layer-2 outer bridge, above the slot's outermost point
                        ("mid" on the bridge circle; bridge-retained rotors only);
    * centre_rib      - on the d-axis under the layer-1 slot's inner pole-side corner (not bridgeless);
    * sleeve          - mid-thickness of the sleeve at 0.55 x half pole pitch below the d-axis;
    * slot_opening    - centre of the layer-2 OD opening (bridgeless rotor only);
    * lip             - middle of the lip between the two OD openings (bridgeless rotor only);
    * air_gap         - mid air gap on the tooth centre line nearest 0.55 x half pole pitch;
    * magnet_length   - middle of the upper layer-2 magnet's pole-side face;
    * cavity          - centre of the central cavity (where the rotor has one);
    * lightening_hole - inside the lower q-axis lightening hole (its label sits bottom left);
    * d_axis          - start of the d-axis label, 1 mm outside the shaft radius.
    """
    p = rp.p
    hp = rp.half_pitch
    L1, L2 = rp.layers
    out: dict[str, np.ndarray] = {}
    if not p.is_bridgeless:
        mid = np.array(L2.p("mid"))
        out["top_bridge"] = mid * (1 + 0.5 * p.pm2_bridge_thickness / np.linalg.norm(mid))
        out["centre_rib"] = np.array([L1.p("P5")[0] - 0.5, 0.0])
    out["sleeve"] = _polar(0.5 * (rp.r_lam + rp.r_sleeve_out), -0.55 * hp)
    if p.is_bridgeless:
        from .geometry.rotor import opening_angles

        lo2, hi2 = opening_angles(p, L2)
        lo1, hi1 = opening_angles(p, L1)
        out["slot_opening"] = _polar(rp.r_lam, 0.5 * (lo2 + hi2))
        out["lip"] = _polar(rp.r_lam - 0.5 * p.lip_depth, 0.5 * (hi1 + lo2))
    tau = 2 * math.pi / p.slots
    teeth = np.concatenate([sp.slot_angles - tau / 2, sp.slot_angles[-1:] + tau / 2])
    teeth = teeth[np.abs(teeth) < hp - 1e-9]
    if len(teeth):
        th = float(teeth[np.argmin(np.abs(teeth - 0.55 * hp))])
        out["air_gap"] = _polar(0.5 * (p.r_sleeve_out + p.r_stator_in), th)
    out["magnet_length"] = 0.5 * (np.array(L2.p("P8")) + np.array(L2.p("P11")))
    if "cavity" in rp.features:
        out["cavity"] = np.array(rp.features["cavity"])
    out["lightening_hole"] = _polar(0.5 * (p.q_hole_inner_radius + p.q_hole_apex_radius), -(hp - 0.02))
    out["d_axis"] = np.array([rp.r_shaft + 1.0, 0.0])
    return out


def _legend_handles() -> tuple[list, list[str]]:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    c = GEOMETRY_COLORS
    pm = Patch(fc=c["magnet"], ec="k", lw=0.4)
    arrow = Line2D([], [], ls="none", marker=r"$\rightarrow$", ms=9, color="k")
    h = [(Patch(fc=c["stator"], ec="k", lw=0.4), "Stator lamination"),
         (Patch(fc=c["rotor"], ec="k", lw=0.4), "Rotor lamination"),
         ((pm, arrow), "Permanent magnet (arrow: magnetisation direction)"),
         (Patch(fc=c["conductor"], ec="k", lw=0.4), "Conductors (coil region)"),
         (Patch(fc=c["sleeve"], ec="k", lw=0.4), "Carbon-fibre sleeve")]
    h.append((Patch(fc=c["air"], ec="k", lw=0.4), "Air (magnet pocket / flux barrier)"))
    return [a for a, _ in h], [b for _, b in h]


def _draw_pole(ax: Axes, rp: RotorPole, sp: StatorPole) -> None:
    from .geometry.primitives import chain_polyline

    c = GEOMETRY_COLORS
    regions = [(sp.loop, c["stator"]), (rp.sleeve, c["sleeve"])]
    for loop, col in regions:
        P = loop.discretize(0.5)
        ax.fill(P[:, 0], P[:, 1], color=col, ec="k", lw=0.4)
    for face in rp.lamination:
        P = chain_polyline(face.outer, 0.5, 32)
        ax.fill(P[:, 0], P[:, 1], color=c["rotor"], ec="k", lw=0.4)
        for h in face.holes:
            P = chain_polyline(h, 0.5, 32)
            ax.fill(P[:, 0], P[:, 1], color=c["air"], ec="k", lw=0.4)
    for mg in rp.magnets:
        P = mg.loop.discretize(10.0)
        ax.fill(P[:, 0], P[:, 1], color=c["magnet"], ec="k", lw=0.4)
        ctr = mg.centre
        ax.annotate("", xy=ctr + MAGNET_ARROW_MM * mg.magnetisation, xytext=ctr,
                    arrowprops=dict(arrowstyle="->", lw=0.8, color="k", shrinkA=0, shrinkB=0))
    for cd in sp.conductors:
        P = cd.loop.pts
        ax.fill(P[:, 0], P[:, 1], color=c["conductor"], ec="k", lw=0.2)


def _draw_callouts(ax: Axes, rp: RotorPole, sp: StatorPole, names: Sequence[str]) -> None:
    anchors = callout_anchors(rp, sp)
    box = dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.9)
    arrow = dict(arrowstyle="-|>", lw=0.6, color="k", shrinkA=1, shrinkB=0, mutation_scale=6)
    if "d_axis" in names:
        x1 = rp.p.stator_od / 2 + 1.0
        ax.plot([rp.r_shaft - 1.0, x1], [0.0, 0.0], color="k", lw=0.6, ls=st.LINE_STYLES[2], zorder=3)
        ax.text(*anchors["d_axis"] + np.array([0.0, 0.6]), CALLOUT_TEXT["d_axis"], ha="left", va="bottom",
                fontsize=st.NOTE_FONTSIZE, style="italic", zorder=4, bbox=box)
    stack = {"tl": 0, "bl": 0}
    for name in names:
        if name == "d_axis" or name not in anchors:
            continue
        corner = _CORNER[name]
        k = stack[corner]
        stack[corner] += 1
        y = 1 - CALLOUT_MARGIN - CALLOUT_ROW_STEP * k if corner == "tl" else CALLOUT_MARGIN + CALLOUT_ROW_STEP * k
        ax.annotate(CALLOUT_TEXT[name], xy=anchors[name], xycoords="data", xytext=(CALLOUT_MARGIN, y),
                    textcoords="axes fraction", ha="left", va="top" if corner == "tl" else "bottom",
                    fontsize=st.NOTE_FONTSIZE, bbox=box, arrowprops=arrow, zorder=6)


def _limits(ax: Axes, rp: RotorPole, sp: StatorPole) -> None:
    P = np.vstack([sp.loop.discretize(1.0), rp.sector.discretize(1.0)])
    pad = LIMIT_PAD_MM
    ax.set_xlim(P[:, 0].min() - pad, P[:, 0].max() + pad)
    ym = np.abs(P[:, 1]).max() + pad
    ax.set_ylim(-ym, ym)


def _fit_headings(fig: Figure, axes: Sequence[Axes]) -> None:
    """Break a panel heading that is wider than its panel over two lines (at ': ' if it has
    one, else at the space nearest the middle)."""
    fig.draw_without_rendering()
    for ax in axes:
        t = ax.title
        if t.get_window_extent().width <= ax.get_window_extent().width or "\n" in t.get_text():
            continue
        s = t.get_text()
        if ": " in s:
            s = s.replace(": ", ":\n", 1)
        else:
            mid = len(s) // 2
            spaces = [i for i, ch in enumerate(s) if ch == " "]
            if spaces:
                i = min(spaces, key=lambda j: abs(j - mid))
                s = s[:i] + "\n" + s[i + 1:]
        t.set_text(s)


def _note_text(p: Params) -> str:
    return (f"{p.slots} slots / {p.poles} poles, one pole ({360 / p.poles:g}°) shown; "
            "dimensions of the author's Maxwell models")


def _fig_legend_and_note(fig: Figure, axes: Sequence[Axes], ncol: int, note: str | None) -> Legend:
    """Legend and note under the panels, placed after measuring the panels' extent."""
    fig.draw_without_rendering()
    inv = fig.transFigure.inverted()
    y0 = min(inv.transform(ax.get_tightbbox().corners())[:, 1].min() for ax in axes)
    handles, labels = _legend_handles()
    from matplotlib.legend_handler import HandlerTuple

    leg = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y0 - 0.01), ncol=ncol,
                     frameon=False, handlelength=1.8, handleheight=1.0, columnspacing=1.2, handletextpad=0.5,
                     handler_map={tuple: HandlerTuple(ndivide=1)})
    if note:
        fig.draw_without_rendering()
        yl = inv.transform(leg.get_window_extent().corners())[:, 1].min()
        fig.text(0.5, yl - 0.01, note, ha="center", va="top", fontsize=st.NOTE_FONTSIZE, style="italic")
    return leg


def _callout_names(callouts: bool | Sequence[str]) -> list[str]:
    if callouts is True:
        return list(CALLOUT_TEXT)
    if isinstance(callouts, str):
        callouts = [callouts]
    names = list(callouts or [])
    unknown = [c for c in names if c not in CALLOUT_TEXT]
    if unknown:
        raise ValueError(f"unknown callout(s) {unknown}; choose from {list(CALLOUT_TEXT)}")
    return names


@st.styled
def geometry_plot(p: Params, out: str | Path | None = None, ax: Axes | None = None, heading: str | None = None,
                  callouts: bool | Sequence[str] = True, legend: bool | None = None,
                  note: str | bool | None = True) -> Figure:
    """Draw one pole of the rotor and stator for a Params object.

    ``heading`` defaults to the topology in words; ``callouts`` is True (every feature that
    exists in this topology), False, or a list of names from :data:`CALLOUT_TEXT`. With
    ``ax=None`` a stand-alone figure is made and ``legend``/``note`` default to on; with an
    ``ax`` they default to off (the caller draws one legend for all panels). Unknown callout
    names raise ValueError; a known feature that this topology lacks is skipped."""
    from .geometry import build_rotor, build_stator

    names = _callout_names(callouts)
    plt = st.pyplot()
    rp, sp = build_rotor(p), build_stator(p)
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(st.COLUMN_W_IN, st.COLUMN_W_IN))
    else:
        fig = ax.figure
    _draw_pole(ax, rp, sp)
    _draw_callouts(ax, rp, sp, names)
    _limits(ax, rp, sp)
    ax.set_aspect("equal")
    ax.set_title(heading if heading is not None else PANEL_HEADINGS.get(p.topology, p.topology),
                 fontsize=st.LABEL_FONTSIZE - 0.5)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    st.style_axes(ax)
    if standalone:
        _fit_headings(fig, [ax])
    if legend if legend is not None else standalone:
        txt = _note_text(p) if note is True else (note or None)
        _fig_legend_and_note(fig, [ax], ncol=2, note=txt)
    return st.save(fig, out)


@st.styled
def topology_figure(out: str | Path | None = None, params: Sequence[Params] | None = None,
                    note: str | bool | None = True, width: float = st.FULL_W_IN,
                    headings: Sequence[str] | None = None,
                    callouts: bool | Sequence[str] | None = None) -> Figure:
    """Side-by-side cross-sections of the rotor topologies (default parameters unless
    ``params`` is given), panel headings (a), (b), ..., one shared legend and a note on what is shown.

    ``headings`` replaces the topology names after the (a), (b), ... prefix, one per panel
    (e.g. the parameter that was changed); ``callouts`` overrides the per-topology default."""
    from .params import TOPOLOGIES, default_params

    plt = st.pyplot()
    params = list(params) if params is not None else [default_params(t) for t in TOPOLOGIES]
    n = len(params)
    if headings is not None and len(headings) != n:
        raise ValueError(f"headings has {len(headings)} entries for {n} panels")
    fig, axs = plt.subplots(1, n, figsize=(width, width / n * 1.12), sharey=True, squeeze=False)
    axs = axs[0]
    for i, (ax, p) in enumerate(zip(axs, params)):
        text = headings[i] if headings is not None else PANEL_HEADINGS.get(p.topology, p.topology)
        heading = f"({chr(ord('a') + i)}) {text}"
        panel_callouts = callouts if callouts is not None else PANEL_CALLOUTS.get(p.topology, True)
        geometry_plot(p, ax=ax, heading=heading, callouts=panel_callouts, legend=False)
        if i:
            ax.set_ylabel("")
    fig.subplots_adjust(wspace=0.06)
    _fit_headings(fig, axs)
    txt = _note_text(params[0]) if note is True else (note or None)
    _fig_legend_and_note(fig, list(axs), ncol=3, note=txt)
    return st.save(fig, out)
