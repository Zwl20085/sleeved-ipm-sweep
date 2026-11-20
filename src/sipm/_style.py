"""House figure style for :mod:`sipm.post` (applied locally, never globally).

Rules: serif (Times New Roman) text with STIX mathtext; axis titles bold and 2 pt larger
than the tick numbers; series told apart by line style first and colour second; open
markers for sparse data only; no plot title; legends that do not cover data; one-column
size 3.5 x 2.7 in by default.

The style is applied through :func:`style`, a ``matplotlib.rc_context``, so importing
``sipm`` or calling a plot function never changes the caller's global ``rcParams``.
"""

from __future__ import annotations

import functools
import os
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.legend import Legend
    from matplotlib.text import Text
    from matplotlib.transforms import Bbox

TICK_FONTSIZE = 8.0
LABEL_FONTSIZE = TICK_FONTSIZE + 2.0
LEGEND_FONTSIZE = TICK_FONTSIZE
NOTE_FONTSIZE = TICK_FONTSIZE
COLUMN_W_IN = 3.5
DEFAULT_H_IN = 2.7
FULL_W_IN = 7.2
DPI = 300
SAVE_PAD_IN = 0.03

# Serif fallbacks in order; DejaVu Serif ships with matplotlib, so the list never runs dry
# (a machine without Times New Roman silently falls back, no warning is logged).
SERIF_FONTS = ["Times New Roman", "Times", "Nimbus Roman", "Nimbus Roman No9 L", "DejaVu Serif"]

# Separable line styles: one plain dash only, the rest carry dots.
LINE_STYLES = [
    "-",                         # solid
    (0, (6, 2)),                 # dashed (the one plain dash)
    (0, (6, 2, 1, 2)),           # dash-dot
    (0, (1, 1.5)),               # dotted
    (0, (6, 2, 1, 2, 1, 2)),     # dash-dot-dot
    (0, (10, 2, 1, 2)),          # long-dash-dot
    (0, (10, 2, 1, 2, 1, 2)),    # long-dash-dot-dot
]
# Open markers, used for sparse series only.
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
# Colour is a secondary cue (colour-blind safe, distinct in grey scale together with the line style).
COLORS = ["#1f4e9c", "#c0392b", "#2e8b57", "#7b3fa0", "#d4800f", "#4d4d4d", "#008b8b"]
REFERENCE_COLOR = "0.35"     # horizontal reference lines
SERIES_MARKER_SIZE = 4.0
SERIES_MARKER_EDGE = 0.9
GRID_KW = dict(alpha=0.3, ls=":", lw=0.6)
NOTE_BBOX = dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=0.9)
LEGEND_HANDLE_LENGTH = 2.6

# Backends that never open a window; any other one is replaced by Agg unless SIPM_INTERACTIVE is set.
NON_INTERACTIVE_BACKENDS = ("agg", "pdf", "svg", "ps", "pgf", "cairo", "template")

RC = {
    "font.family": ["serif"],
    "font.serif": SERIF_FONTS,
    "mathtext.fontset": "stix",
    "axes.labelsize": LABEL_FONTSIZE,
    "axes.labelweight": "bold",
    "axes.titlesize": LABEL_FONTSIZE,
    "axes.linewidth": 0.8,
    "xtick.labelsize": TICK_FONTSIZE,
    "ytick.labelsize": TICK_FONTSIZE,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": LEGEND_FONTSIZE,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.framealpha": 1.0,
    "legend.fancybox": False,
    "legend.handlelength": LEGEND_HANDLE_LENGTH,
    "lines.linewidth": 1.3,
    "lines.markersize": 4.5,
    "hatch.linewidth": 0.8,
    "savefig.dpi": DPI,
}

# Note positions tried in order (axes fraction x, y, ha, va).
NOTE_CANDIDATES = (
    (0.5, 0.03, "center", "bottom"), (0.03, 0.03, "left", "bottom"), (0.97, 0.03, "right", "bottom"),
    (0.03, 0.97, "left", "top"), (0.97, 0.97, "right", "top"), (0.5, 0.97, "center", "top"),
    (0.03, 0.5, "left", "center"), (0.97, 0.5, "right", "center"),
)
NOTE_CLEARANCE = (1.05, 1.15)  # note box grown by these factors (x, y) before the overlap test


def _configured_backend() -> str | None:
    """The backend name set by the user, or None while matplotlib would still auto-select.

    ``matplotlib.get_backend()`` resolves the automatic choice by importing Qt, Tk, ... in
    turn; a broken Qt install then crashes inside the import (Windows: 'fatal exception
    0xc0000139' under faulthandler). Reading the raw rcParams entry avoids that probe."""
    import matplotlib

    raw = dict.__getitem__(matplotlib.rcParams, "backend")
    return raw.lower() if isinstance(raw, str) else None


def pyplot() -> Any:
    """``matplotlib.pyplot``, switching to the Agg backend unless SIPM_INTERACTIVE is set or a
    non-interactive backend is already selected. The GUI backends are never probed."""
    import matplotlib

    if _configured_backend() not in NON_INTERACTIVE_BACKENDS and not os.environ.get("SIPM_INTERACTIVE"):
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def save(fig: Figure, out: str | Path | None) -> Figure:
    """Write ``fig`` to ``out`` (if given) at :data:`DPI` with a tight bounding box."""
    if out:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=SAVE_PAD_IN)
    return fig


@contextmanager
def style() -> Iterator[None]:
    """Context in which the house style is active; global rcParams are restored on exit."""
    import matplotlib

    with matplotlib.rc_context(RC):
        yield


def styled(func: Callable[..., Figure]) -> Callable[..., Figure]:
    """Decorator: run a plot function inside :func:`style` and lay the figure out before
    leaving it, so tick labels are created with the house fonts."""

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Figure:
        with style():
            fig = func(*args, **kwargs)
            if fig is not None:
                fig.draw_without_rendering()
            return fig

    return wrapper


def style_axes(ax: Axes) -> None:
    """Inward ticks on all four sides, thin spines; baked into the axes so later redraws
    outside the style context keep them."""
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=TICK_FONTSIZE)
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)


def series_style(i: int) -> dict:
    """Line style, colour and open marker for series ``i``."""
    return dict(linestyle=LINE_STYLES[i % len(LINE_STYLES)], color=COLORS[i % len(COLORS)],
                marker=MARKERS[i % len(MARKERS)], markerfacecolor="none")


def legend_above(ax: Axes, ncol: int = 2, title: str | None = None, handles: Sequence | None = None,
                 labels: Sequence[str] | None = None, anchor_y: float = 1.02, left: bool = False) -> Legend:
    """Legend parked above the axes: it can never cover data. ``left`` aligns it (and its
    title) with the left edge of the axes instead of centring it. ``ncol`` below 1 is read as 1."""
    loc, x = ("lower left", 0.0) if left else ("lower center", 0.5)
    kw = dict(loc=loc, bbox_to_anchor=(x, anchor_y), ncol=max(1, int(ncol)), handlelength=LEGEND_HANDLE_LENGTH,
              handletextpad=0.5, columnspacing=1.0, borderaxespad=0.0, frameon=False,
              title=title, title_fontsize=NOTE_FONTSIZE, alignment="left" if left else "center")
    if handles is not None:
        return ax.legend(handles, labels if labels is not None else [h.get_label() for h in handles], **kw)
    return ax.legend(**kw)


def title_above_legend(ax: Axes, title: str | None, leg: Legend | None = None) -> None:
    """Optional axes title, raised clear of a legend parked above the axes."""
    if not title:
        return
    pad = 6.0
    if leg is not None:
        fig = ax.figure
        fig.draw_without_rendering()
        h = leg.get_window_extent().height / fig.dpi * 72.0
        pad += h + 4.0
    ax.set_title(title, pad=pad, fontsize=LABEL_FONTSIZE)


def place_note(ax: Axes, text: str, points: np.ndarray | Sequence | None,
               candidates: Sequence[tuple[float, float, str, str]] | None = None, avoid: Sequence = ()) -> Text:
    """Write ``text`` inside ``ax`` at the first candidate spot whose box covers none of
    ``points`` (data coordinates, N x 2), no line of the axes and none of the ``avoid``
    artists. Tries each candidate once, so it always terminates. If every spot is taken the
    note becomes the title of the axes legend (when it has none), else it goes above the axes.
    Axis limits must be final before the call. Returns the Text that holds the note."""
    fig = ax.figure
    fig.draw_without_rendering()  # fixes the transforms once; candidates are then measured, not drawn
    obstacles = _obstacle_points(ax, points)
    avoid_boxes = [a.get_window_extent() for a in avoid if a is not None]
    for x, y, ha, va in candidates or NOTE_CANDIDATES:
        t = ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=NOTE_FONTSIZE, bbox=NOTE_BBOX,
                    zorder=5)
        box = t.get_window_extent().expanded(*NOTE_CLEARANCE)
        if not _any_inside(box, obstacles) and not any(box.overlaps(b) for b in avoid_boxes):
            return t
        t.remove()
    leg = ax.get_legend()
    if leg is not None and not leg.get_title().get_text():
        leg.set_title(text, prop=dict(size=NOTE_FONTSIZE))
        return leg.get_title()
    return ax.text(0.0, 1.02, text, transform=ax.transAxes, ha="left", va="bottom", fontsize=NOTE_FONTSIZE)


def _obstacle_points(ax: Axes, points: np.ndarray | Sequence | None) -> np.ndarray:
    """Display coordinates of the data points and of every visible line of ``ax`` (densified)."""
    pts = np.asarray(points, float).reshape(-1, 2) if points is not None and len(points) else np.empty((0, 2))
    parts = [ax.transData.transform(pts)] if len(pts) else []
    for ln in ax.get_lines():
        xy = ln.get_xydata()
        if ln.get_visible() and len(xy):
            parts.append(_densify(ln.get_transform().transform(xy)))
    return np.vstack(parts) if parts else np.empty((0, 2))


def _any_inside(box: Bbox, P: np.ndarray) -> bool:
    P = np.asarray(P, float).reshape(-1, 2)
    inside = (P[:, 0] >= box.x0) & (P[:, 0] <= box.x1) & (P[:, 1] >= box.y0) & (P[:, 1] <= box.y1)
    return bool(inside.any())


def _densify(P: np.ndarray, step: float = 2.0) -> np.ndarray:
    """Points every ``step`` display pixels along a polyline (so a box between vertices is caught)."""
    P = np.asarray(P, float)
    P = P[np.isfinite(P).all(axis=1)]
    if len(P) < 2:
        return P
    out = [P[:1]]
    for a, b in zip(P[:-1], P[1:]):
        n = max(1, int(np.hypot(*(b - a)) / step))
        s = np.linspace(0, 1, n + 1)[1:, None]
        out.append(a + s * (b - a))
    return np.vstack(out)
