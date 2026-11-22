"""Plots of sweep results (matplotlib, Agg backend safe).

* :func:`torque_scatter`  - average torque against swept parameters, coloured by torque
                            (the layout of the paper's Fig. 3);
* :func:`mtpa_lines`      - torque against current advance angle, one line per design;
* :func:`stress_lines`    - peak stress against interference, one line per sleeve thickness
                            (layout of Figs. 4 and 6);
* :func:`stress_map`      - stress interpolated over (sleeve thickness, interference) with the
                            strength contour, the interference <= fraction x thickness limit
                            and the feasible region (layout of Figs. 5 and 7);
* :func:`torque_waveforms`, :func:`kac_lines`;
* :func:`geometry_plot`, :func:`topology_figure` - rotor/stator cross-sections.

Every function takes a results DataFrame (``Sweep.results()`` or the CSV it writes) and
returns the matplotlib Figure; pass ``out=`` to save it. The figure stays open: close it
(``matplotlib.pyplot.close(fig)``) when drawing many. Figures follow the house style in
:mod:`sipm._style` (serif fonts, bold axis titles, series told apart by line style, open
markers, legends that do not cover data, no title unless one is passed). The style is applied
with ``matplotlib.rc_context``; the caller's global ``rcParams`` are left untouched. Input
tables, label dicts and parameter objects are never modified.
"""

from __future__ import annotations

import math
import textwrap
import warnings
from collections.abc import Hashable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from . import _style as st
from ._geomfig import GEOMETRY_COLORS, PANEL_HEADINGS, geometry_plot, topology_figure

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

__all__ = ["LABELS", "TOPOLOGY_NAMES", "PANEL_HEADINGS", "GEOMETRY_COLORS", "torque_scatter", "mtpa_lines",
           "stress_lines", "stress_map", "feasible_points", "torque_waveforms", "kac_lines", "geometry_plot",
           "topology_figure", "series_label", "conditions_note"]

NM = r"N$\cdot$m"

# geometry parameters a sweep typically varies (new names first; the old names identify the
# columns of results written before the geometry rebase)
GEOMETRY_SWEEP_COLUMNS = ("pm1_bridge_thickness", "pm2_bridge_thickness", "pm1_rib_half_width", "pm2_web_width",
                          "lip_depth", "pm1_opening_angle", "pm2_opening_angle",
                          "pm1_bridge_top", "pm1_bridge_bottom", "pm2_bridge_top", "pm2_bridge_bottom",
                          "rotor_slot_opening")

LABELS = {
    "gamma": "Current advance angle from +d axis (elec. deg)",
    "sleeve_thickness": "Sleeve thickness (mm)",
    "interference": "Radial interference (mm)",
    "pm1_bridge_thickness": "Layer-1 outer bridge (mm)",
    "pm2_bridge_thickness": "Layer-2 outer bridge (mm)",
    "pm1_rib_half_width": "Layer-1 half centre rib (mm)",
    "pm2_web_width": "Layer-2 web beside cavity (mm)",
    "lip_depth": "Lip depth (mm)",
    "pm1_opening_angle": "Layer-1 OD opening (deg)",
    "pm2_opening_angle": "Layer-2 OD opening (deg)",
    # columns of results written before the geometry rebase (examples/data/sample_*.csv)
    "rotor_slot_opening": "Rotor slot opening (mm)",
    "pm1_bridge_top": "Layer-1 top bridge (mm)",
    "pm1_bridge_bottom": "Layer-1 bottom bridge (mm)",
    "pm2_bridge_top": "Layer-2 top bridge (mm)",
    "pm2_bridge_bottom": "Layer-2 bottom bridge (mm)",
    "current_rms": "Phase current (A rms)",
    "em_rpm": "Speed (rpm)",
    "mech_rpm": "Speed (rpm)",
    "torque_avg": f"Average torque ({NM})",
    "svm_rotor_max": "Max von Mises at named lamination sections (MPa)",
    "svm_lam_max": "Max von Mises at named lamination sections (MPa)",
    "svm_lam_max_all": "Max von Mises, whole lamination, mesh-dependent (MPa)",
    "svm_pm_max": "Max von Mises, PM (MPa)",
    "sh_sleeve_max": "Max hoop stress, sleeve (MPa)",
    "sr_sleeve_min": "Min radial stress, sleeve (MPa)",
    "sr_sleeve_max": "Max radial tension, sleeve (MPa)",
    "tau_sleeve_max": "Max radial-hoop shear, sleeve (MPa)",
    "mech_temperature": "Rotor temperature (°C)",
    "sleeve_wind_angle": "Sleeve winding angle from hoop (deg)",
    "sleeve_vf": "Sleeve fibre volume fraction (-)",
    "sleeve_fi_max_stress": "Max-stress failure index, sleeve (-)",
    "sleeve_fi_tsai_wu": "Tsai-Wu failure index (1/strength ratio), sleeve (-)",
    "sleeve_fi_hashin_fibre": "Hashin fibre failure index, sleeve (-)",
    "sleeve_fi_hashin_matrix": "Hashin matrix failure index, sleeve (-)",
    "sleeve_sf": "Tsai-Wu strength ratio (safety factor), sleeve (-)",
}

# Words used for a topology in legends (series names).
TOPOLOGY_NAMES = {
    "conventional_baseline": "Conventional baseline (author's rotor, not in the package)",
    "sleeved_bridged": "Sleeved bridge-retained",
    "sleeved_bridgeless": "Sleeved bridgeless",
}

AXIS_LABEL_WRAP = 38          # characters per line of an axis or colour-bar title
FEASIBLE_TOL = 1e-12          # tolerance of the interference limit test (mm)
THOUSANDS_FROM = 10000        # speeds from this value on are grouped "20 000"
THIN_SPACE = " "         # thousands separator in speeds
STAR_KW = dict(ls="none", marker="*", ms=9, mfc="none", mec="k", mew=0.9)
WASH_COLOR = "0.92"           # infeasible-interference wash of the stress map

_plt = st.pyplot
_save = st.save

OutPath = str | Path | None  # where to save a figure; None: do not save


# ----------------------------------------------------------------------------- text helpers
def _ok(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Rows with status 'ok' (if there is a status column) and a finite ``col``."""
    if col not in df.columns:
        raise KeyError(f"results table has no column {col!r}; columns: {list(df.columns)}")
    d = df[df["status"] == "ok"] if "status" in df.columns else df
    return d[np.isfinite(pd.to_numeric(d[col], errors="coerce"))]


def _wrap(label: str, width: int = AXIS_LABEL_WRAP) -> str:
    """Break a long axis title over lines (mathtext segments are never split)."""
    if len(label) <= width or "$" in label:
        return label
    return "\n".join(textwrap.wrap(label, width))


def _num(v: Any) -> str:
    v = float(v)
    return f"{v:.1f}" if v.is_integer() else f"{v:g}"


def _is_missing(v: Any) -> bool:
    """None, NaN, NaT or pd.NA (scalars only)."""
    return v is None or (np.ndim(v) == 0 and bool(pd.isna(v)))


def _split_label(col: str) -> tuple[str, str]:
    """'Sleeve thickness (mm)' -> ('Sleeve thickness', 'mm'); unknown columns -> (col, '')."""
    lab = LABELS.get(col)
    if not lab or not lab.endswith(")") or "(" not in lab:
        return col.replace("_", " "), ""
    name, unit = lab[: lab.rindex("(")].strip(), lab[lab.rindex("(") + 1: -1]
    return name, unit


def _capitalise(text: str) -> str:
    return text[:1].upper() + text[1:]


def series_label(group: Sequence[str], keys: Any) -> str:
    """Legend text in words for one group of a ``groupby(group)``.

    ``topology`` becomes its name ("Sleeved bridge-retained"), ``sleeve_thickness`` becomes
    "1.0 mm sleeve" ("no sleeve" for 0), any other column "<name> <value> <unit>"."""
    keys = keys if isinstance(keys, tuple) else (keys,)
    vals = dict(zip(group, keys))
    parts = []
    topo = vals.get("topology")
    if not _is_missing(topo):
        parts.append(TOPOLOGY_NAMES.get(str(topo), str(topo)))
    for col, v in vals.items():
        if col == "topology" or _is_missing(v):
            continue
        if col == "sleeve_thickness":
            parts.append("no sleeve" if float(v) == 0 else f"{_num(v)} mm sleeve")
        elif isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_)):
            name, unit = _split_label(col)
            parts.append(f"{name.lower() if parts else name} {_num(v)} {unit}".strip())
        else:
            parts.append(str(v))
    return _capitalise(", ".join(parts) if parts else "all")


def _label_candidates(keys: Any) -> list:
    keys_t = keys if isinstance(keys, tuple) else (keys,)
    return [keys_t, keys_t[0]] if len(keys_t) == 1 else [keys_t]


def _series_labels(labels: Mapping | None, group: Sequence[str], all_keys: Sequence) -> list[str]:
    """Legend text per group key: ``labels[key]`` (scalar or tuple key) or :func:`series_label`.
    Warns about ``labels`` keys that match no group, which are usually typos or float noise."""
    used: set = set()
    texts = []
    for keys in all_keys:
        text = None
        for k in _label_candidates(keys) if labels else ():
            try:
                hit = k in labels
            except TypeError:  # unhashable key
                continue
            if hit:
                text = labels[k]
                used.add(k)
                break
        texts.append(text if text is not None else series_label(group, keys))
    unused = [k for k in (labels or {}) if k not in used]
    if unused:
        warnings.warn(f"labels keys {unused} match no group of {list(group)}; their default text is used",
                      UserWarning, stacklevel=3)
    return texts


def _speed(v: float) -> str:
    """'3000', '20 000' (thin-space U+2009 thousands grouping from five digits on)."""
    return f"{int(round(v)):,}".replace(",", THIN_SPACE) if v >= THOUSANDS_FROM else f"{v:g}"


def _single_value(df: pd.DataFrame, col: str) -> Any:
    if col not in df.columns:
        return None
    vals = df[col].dropna().unique()
    return vals[0] if len(vals) == 1 else None


def conditions_note(df: pd.DataFrame, physics: str = "em") -> str | None:
    """Operating conditions shared by every row of a results table, in words, for the
    ``note`` argument: rotor topology, phase current (EM) and speed. A quantity that varies
    or is missing is left out; None if nothing is common."""
    parts = []
    topo = _single_value(df, "topology")
    if topo is not None:
        parts.append(f"{TOPOLOGY_NAMES.get(str(topo), str(topo))} rotor")
    if physics == "em":
        cur = _single_value(df, "current_rms")
        if cur is not None:
            parts.append(f"phase current {float(cur):g} A rms")
    rpm = _single_value(df, "mech_rpm" if physics == "mech" else "em_rpm")
    if rpm is not None:
        parts.append(f"{_speed(float(rpm))} rpm")
    return _capitalise(", ".join(parts)) if parts else None


# ----------------------------------------------------------------------------- layout helpers
def _headroom(ax: Axes, ys: Sequence[float], top: float = 0.16, bottom: float = 0.06, zero_floor: bool = False) -> None:
    ys = np.asarray(ys, float)
    ys = ys[np.isfinite(ys)]
    if not len(ys):
        return
    lo, hi = ys.min(), ys.max()
    span = hi - lo if hi > lo else max(abs(hi), 1.0)
    lo2 = lo - bottom * span
    if zero_floor and lo >= 0:
        lo2 = max(0.0, lo2)
    ax.set_ylim(lo2, hi + top * span)


def _reference_handle(label: str, ls: Any) -> Any:
    from matplotlib.lines import Line2D

    return Line2D([], [], color=st.REFERENCE_COLOR, lw=1.0, ls=ls, label=label)


def _finish(ax: Axes, out: OutPath, title: str | None, note: str | None, pts: Sequence,
            handles: Sequence | None = None, ncol: int = 2, avoid: Sequence = ()) -> Figure:
    """Common tail of the single-axes line plots: grid, ticks, legend above the axes, the
    note inside the axes where it covers no data, the optional title, save."""
    ax.grid(**st.GRID_KW)
    st.style_axes(ax)
    handles = list(handles) if handles is not None else ax.get_legend_handles_labels()[0]
    leg = st.legend_above(ax, ncol=ncol, handles=handles) if handles else None  # nothing drawn: no legend
    if note:
        st.place_note(ax, note, np.asarray(pts, float).reshape(-1, 2), avoid=avoid)
    st.title_above_legend(ax, title, leg)
    return _save(ax.figure, out)


# ----------------------------------------------------------------------------- torque plots
def _scatter_panel(ax: Axes, d: pd.DataFrame, x: str, norm: Any, cmap: Any, reference: float | None) -> bool:
    """One panel of :func:`torque_scatter`; True if it drew the gamma envelope."""
    ax.scatter(d[x], d["torque_avg"], s=14, facecolors="none", edgecolors=cmap(norm(d["torque_avg"])),
               linewidths=0.8)
    envelope = x == "gamma" and len(d) > 0
    if envelope:
        env = d.groupby("gamma")["torque_avg"].max()
        ax.plot(env.index, env.values, color="k", lw=1.1, ls=st.LINE_STYLES[1])
    if reference is not None:
        ax.axhline(reference, color=st.REFERENCE_COLOR, lw=1.0, ls=st.LINE_STYLES[2])
    ax.set_xlabel(_wrap(LABELS.get(x, x), 30))
    _headroom(ax, list(d["torque_avg"]) + ([reference] if reference is not None else []), top=0.06)
    st.style_axes(ax)
    return envelope


@st.styled
def torque_scatter(df: pd.DataFrame, xcols: list[str], reference: float | None = None,
                   reference_label: str = "reference design peak", out: OutPath = None,
                   title: str | None = None) -> Figure:
    """Average torque against each column in ``xcols`` (first panel usually 'gamma').

    In the gamma panel the dashed line is the upper envelope (MTPA peak per angle)."""
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    if not xcols:
        raise ValueError("torque_scatter needs at least one column in xcols")
    plt = _plt()
    d = _ok(df, "torque_avg")
    n, ncol = len(xcols), min(3, len(xcols))
    nrow = math.ceil(n / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.6 * ncol + 0.9, 2.4 * nrow), squeeze=False, sharey=True)
    norm = Normalize(vmin=d["torque_avg"].min(), vmax=d["torque_avg"].max()) if len(d) else Normalize(0.0, 1.0)
    cmap = plt.get_cmap("coolwarm")
    handles = [Line2D([], [], ls="none", marker="o", mfc="none", mec="0.3", ms=4, label="Solved design points")]
    envelope = False
    for j, (ax, x) in enumerate(zip(axs.flat, xcols)):
        envelope |= _scatter_panel(ax, d, x, norm, cmap, reference)
        if j % ncol == 0:
            ax.set_ylabel(LABELS["torque_avg"])
    if envelope:
        handles.append(Line2D([], [], color="k", lw=1.1, ls=st.LINE_STYLES[1],
                              label="Upper envelope (max over designs per angle)"))
    if reference is not None:
        handles.append(_reference_handle(reference_label, st.LINE_STYLES[2]))
    for ax in list(axs.flat)[n:]:
        ax.axis("off")
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.85, pad=0.02)
    cb.set_label(f"Average torque ({NM}), point colour", fontweight="bold", fontsize=st.LABEL_FONTSIZE)
    cb.ax.tick_params(labelsize=st.TICK_FONTSIZE, direction="in")
    leg = st.legend_above(axs[0, 0], ncol=1, handles=handles, left=True)
    st.title_above_legend(axs[0, 0], title, leg)
    return _save(fig, out)


def _mark_peak(ax: Axes, g: pd.DataFrame) -> Any:
    """Open star and angle at the highest torque of one curve (positional, so a repeated
    index label cannot select several rows)."""
    k = int(np.argmax(g["torque_avg"].to_numpy(float)))
    gx, gy = float(g["gamma"].iloc[k]), float(g["torque_avg"].iloc[k])
    ax.plot(gx, gy, zorder=4, **STAR_KW)
    return ax.annotate(f"{gx:g}°", (gx, gy), xytext=(0, 5), textcoords="offset points", ha="center",
                       va="bottom", fontsize=st.NOTE_FONTSIZE - 0.5)


@st.styled
def mtpa_lines(df: pd.DataFrame, group: list[str] | str = "sleeve_thickness", out: OutPath = None,
               title: str | None = None, labels: Mapping[Hashable, str] | None = None, note: str | None = None,
               figsize: tuple[float, float] = (st.COLUMN_W_IN, st.DEFAULT_H_IN)) -> Figure:
    """Average torque against current advance angle, one line per design.

    ``labels`` maps a group key (the value, or the tuple of values for several group columns)
    to its legend text; missing keys get a text built by :func:`series_label`, and keys that
    match no group raise a warning. ``note`` (e.g. "420 A rms, 3000 rpm") is written inside
    the axes where it covers no data. The highest solved torque of each curve is marked with
    an open star and its angle."""
    from matplotlib.lines import Line2D

    plt = _plt()
    d = _ok(df, "torque_avg")
    group = [group] if isinstance(group, str) else list(group)
    fig, ax = plt.subplots(figsize=figsize)
    groups = list(d.groupby(group, dropna=False, sort=True))
    texts = _series_labels(labels, group, [k for k, _ in groups])
    handles, pts, annotations = [], [], []
    for i, ((_, g), text) in enumerate(zip(groups, texts)):
        g = g.sort_values("gamma")
        (h,) = ax.plot(g["gamma"], g["torque_avg"], ms=st.SERIES_MARKER_SIZE, mew=st.SERIES_MARKER_EDGE,
                       label=text, **st.series_style(i))
        handles.append(h)
        annotations.append(_mark_peak(ax, g))
        pts.extend(zip(g["gamma"], g["torque_avg"]))
    handles.append(Line2D([], [], label="Maximum torque per curve (angle in deg)", **STAR_KW))
    ax.set_xlabel(_wrap(LABELS["gamma"], 46))
    ax.set_ylabel(LABELS["torque_avg"])
    _headroom(ax, d["torque_avg"], top=0.14, bottom=0.08)
    if len(d):
        gmin, gmax = d["gamma"].min(), d["gamma"].max()
        pad = 0.04 * (gmax - gmin if gmax > gmin else 1.0)
        ax.set_xlim(gmin - pad, gmax + pad)
    return _finish(ax, out, title, note, pts, handles=handles, ncol=1 if len(handles) > 4 else 2,
                   avoid=annotations)


# ----------------------------------------------------------------------------- stress plots
@st.styled
def stress_lines(df: pd.DataFrame, quantity: str = "svm_rotor_max", reference: float | None = None,
                 reference_label: str = "reference design", out: OutPath = None, title: str | None = None,
                 note: str | None = None) -> Figure:
    """Peak ``quantity`` against radial interference, one line per sleeve thickness."""
    plt = _plt()
    d = _ok(df, quantity)
    fig, ax = plt.subplots(figsize=(st.COLUMN_W_IN, st.DEFAULT_H_IN))
    handles, pts = [], []
    groups = list(d.groupby("sleeve_thickness"))
    for j, (t, g) in enumerate(groups):
        g = g.sort_values("interference")
        (h,) = ax.plot(g["interference"], g[quantity], ms=st.SERIES_MARKER_SIZE, mew=st.SERIES_MARKER_EDGE,
                       label=series_label(["sleeve_thickness"], t), **st.series_style(j))
        handles.append(h)
        pts.extend(zip(g["interference"], g[quantity]))
    ys = list(d[quantity])
    if reference is not None:
        ls = st.LINE_STYLES[len(groups) % len(st.LINE_STYLES)]
        ax.axhline(reference, color=st.REFERENCE_COLOR, lw=1.0, ls=ls)
        handles.append(_reference_handle(f"{reference_label} ({reference:g} MPa)", ls))
        ys.append(reference)
    ax.set_xlabel(LABELS["interference"])
    ax.set_ylabel(_wrap(LABELS.get(quantity, quantity), 30))
    _headroom(ax, ys, top=0.08, bottom=0.08)
    return _finish(ax, out, title, note, pts, handles=handles, ncol=min(3, len(handles)))


def _interpolate(x: np.ndarray, y: np.ndarray, z: np.ndarray, n: int, method: str) -> tuple[np.ndarray, ...]:
    """``griddata`` of z over an n x n grid spanning the solved points; clear errors for
    tables that cannot span an area."""
    from scipy.interpolate import griddata

    try:
        from scipy.spatial import QhullError
    except ImportError:  # pragma: no cover - scipy < 1.10 keeps it in the private module
        from scipy.spatial.qhull import QhullError

    if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        raise ValueError("stress_map needs at least two sleeve thicknesses and two interference values")
    X, Y = np.meshgrid(np.linspace(x.min(), x.max(), n), np.linspace(y.min(), y.max(), n))
    try:
        Z = griddata((x, y), z, (X, Y), method=method)
    except QhullError as e:
        raise ValueError("stress_map cannot interpolate: the solved (sleeve thickness, interference) points "
                         "are collinear, so they span no area; solve a grid of points") from e
    return X, Y, Z


def _draw_limits(ax: Axes, X: np.ndarray, Y: np.ndarray, Z: np.ndarray, sigma_allow: float,
                 fraction: float) -> list:
    """Strength contour, interference limit, its wash and the feasible hatch; their legend handles."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    handles = []
    stress_ls, limit_ls = st.LINE_STYLES[2], st.LINE_STYLES[1]
    if np.isfinite(Z).any() and np.nanmin(Z) <= sigma_allow <= np.nanmax(Z):
        ax.contour(X, Y, Z, levels=[sigma_allow], colors="k", linestyles=[stress_ls], linewidths=1.5)
        handles.append(Line2D([], [], color="k", lw=1.5, ls=stress_ls,
                              label=f"Stress = {sigma_allow:g} MPa (allowable, user input)"))
    x_row, limit_row = X[0], fraction * X[0]
    ax.plot(x_row, limit_row, color="k", lw=1.2, ls=limit_ls)
    handles.append(Line2D([], [], color="k", lw=1.2, ls=limit_ls,
                          label=f"Interference = {fraction:g} × sleeve thickness (limit)"))
    if (limit_row < Y.max()).any():
        ax.fill_between(x_row, np.minimum(limit_row, Y.max()), Y.max(), lw=0, color=WASH_COLOR, alpha=0.75)
        handles.append(Patch(facecolor=WASH_COLOR, alpha=0.9, edgecolor="0.5", lw=0.5,
                             label="Interference above limit (washed out)"))
    feasible = ~((Y > fraction * X + FEASIBLE_TOL) | ~(Z <= sigma_allow))
    if feasible.any():
        ax.contourf(X, Y, np.where(feasible, 1.0, np.nan), levels=[0.5, 1.5], colors="none", hatches=["///"])
        ax.contour(X, Y, feasible.astype(float), levels=[0.5], colors="k", linewidths=0.8)
        handles.append(Patch(facecolor="none", edgecolor="k", hatch="///", lw=0.8,
                             label="Feasible region (both limits met)"))
    return handles


@st.styled
def stress_map(df: pd.DataFrame, sigma_allow: float, interference_fraction: float = 0.2,
               quantity: str = "svm_rotor_max", method: str = "linear", n: int = 121, out: OutPath = None,
               title: str | None = None, note: str | None = None,
               figsize: tuple[float, float] = (st.COLUMN_W_IN, 3.0)) -> Figure:
    """Interpolated stress over (sleeve thickness, interference).

    Feasible region: stress <= ``sigma_allow`` and interference <= ``interference_fraction``
    x sleeve thickness. Both thresholds are user inputs; no default strength is assumed.
    Interpolation is ``scipy.interpolate.griddata`` (no extrapolation outside the hull).
    Every drawn element has a legend entry; the legend (with ``note`` as its heading, e.g.
    the rotor and speed) sits above the axes."""
    from matplotlib.lines import Line2D

    plt = _plt()
    d = _ok(df, quantity)
    x, y, z = d["sleeve_thickness"].to_numpy(float), d["interference"].to_numpy(float), d[quantity].to_numpy(float)
    X, Y, Z = _interpolate(x, y, z, n, method)
    fig, ax = plt.subplots(figsize=figsize)
    cf = ax.contourf(X, Y, Z, levels=14, cmap="coolwarm")
    cb = fig.colorbar(cf, ax=ax, pad=0.03)
    cb.set_label(_wrap(LABELS.get(quantity, quantity), 34), fontweight="bold", fontsize=st.LABEL_FONTSIZE - 1)
    cb.ax.tick_params(labelsize=st.TICK_FONTSIZE, direction="in")
    handles = _draw_limits(ax, X, Y, Z, sigma_allow, interference_fraction)
    # open circles (a '+' on a spine reads like a tick mark); clip_on=False keeps edge points visible
    ax.scatter(x, y, s=16, facecolors="none", edgecolors="k", marker="o", linewidths=0.8, zorder=4, clip_on=False)
    handles.append(Line2D([], [], ls="none", marker="o", mfc="none", mec="k", ms=4, mew=0.8, label="Solved points"))
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_xlabel(LABELS["sleeve_thickness"])
    ax.set_ylabel(LABELS["interference"])
    st.style_axes(ax)
    leg = st.legend_above(ax, ncol=1, handles=handles, title=note, anchor_y=1.03, left=True)
    st.title_above_legend(ax, title, leg)
    return _save(fig, out)


def feasible_points(df: pd.DataFrame, sigma_allow: float, interference_fraction: float = 0.2,
                    quantity: str = "svm_rotor_max") -> pd.DataFrame:
    """Solved rows that meet both limits of :func:`stress_map` (a filtered copy of ``df``)."""
    d = _ok(df, quantity)
    limit = interference_fraction * d["sleeve_thickness"] + FEASIBLE_TOL
    return d[(d[quantity] <= sigma_allow) & (d["interference"] <= limit)]


# ----------------------------------------------------------------------------- other plots
@st.styled
def torque_waveforms(waves: Mapping[str, Mapping[str, Sequence[float]]], p_poles: int, rpm: float,
                     out: OutPath = None, title: str | None = None, note: str | None = None) -> Figure:
    """``waves``: label -> {'Time': [...], 'Moving1.Torque': [...]} (SI)."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(st.COLUMN_W_IN, st.DEFAULT_H_IN))
    fe = rpm * p_poles / 2 / 60.0
    ys, pts = [], []
    for i, (lab, w) in enumerate(waves.items()):
        t = np.asarray(w["Time"], float)
        tq = np.asarray(w["Moving1.Torque"], float)
        sty = st.series_style(i)
        ax.plot(360 * fe * t, tq, ls=sty["linestyle"], color=sty["color"], lw=1.1, label=lab)
        ys.extend(tq)
        pts.extend(zip(360 * fe * t, tq))
    ax.set_xlabel("Rotor position (elec. deg)")
    ax.set_ylabel(f"Torque ({NM})")
    _headroom(ax, ys, top=0.08, bottom=0.08)
    return _finish(ax, out, title, note, pts, ncol=min(3, len(waves)))


@st.styled
def kac_lines(df: pd.DataFrame, group: str = "topology", out: OutPath = None,
              labels: Mapping[Hashable, str] | None = None, note: str | None = None,
              title: str | None = None) -> Figure:
    """AC loss factor against speed, one line per ``group`` value."""
    plt = _plt()
    d = _ok(df, "kac")
    fig, ax = plt.subplots(figsize=(st.COLUMN_W_IN, st.DEFAULT_H_IN))
    groups = list(d.groupby(group))
    texts = _series_labels(labels, [group], [k for k, _ in groups])
    pts = []
    for i, ((_, g), text) in enumerate(zip(groups, texts)):
        g = g.sort_values("em_rpm")
        ax.plot(g["em_rpm"], g["kac"], ms=st.SERIES_MARKER_SIZE, mew=st.SERIES_MARKER_EDGE, label=text,
                **st.series_style(i))
        pts.extend(zip(g["em_rpm"], g["kac"]))
    ax.set_xlabel("Rotor speed (rpm)")
    ax.set_ylabel("AC loss factor $p_\\mathrm{ac}/p_\\mathrm{dc}$, active length (-)")
    _headroom(ax, d["kac"], top=0.08, bottom=0.08)
    return _finish(ax, out, title, note, pts, ncol=min(3, len(groups)))
