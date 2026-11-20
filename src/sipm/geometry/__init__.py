"""Licence-free 2D geometry of one rotor pole and one stator pole pitch.

Frame ("pole frame"): origin on the shaft axis, d-axis of the modelled pole along +x, the
pole spans [-pi/p, +pi/p] (p = pole pairs). Lengths in mm.

The rotor is constructed the way the author's Maxwell models construct it (see
docs/methods.md, section 1): each magnet layer is a straight slot fixed by the d-axis
crossing of its magnet centre line, its angle and width; the magnet length is derived. Module
map:

* :mod:`.primitives`   exact segments (line, arc, cubic Hermite), chains, fillets, biarcs and
                       the solver-facing :class:`Loop`;
* :mod:`.csg`          exact Boolean evaluation of line/arc regions;
* :mod:`.construction` the per-layer ``RD_*`` construction points and the derived magnet length;
* :mod:`.rotor`        the two sleeved rotor topologies;
* :mod:`.checks`       feasibility checks and derived metrics;
* :mod:`.locations`    named load-carrying sections where the headline stress is read;
* :mod:`.stator`       stator pole pitch, conductors and winding layout;
* :mod:`.polygons`     sampled-polygon helpers.
"""

from .construction import LayerInput, LayerPoints, layer_points
from .polygons import fillet, points_in_polygon, polyline_distance, seg_seg_distance, self_intersects
from .primitives import (
    BIARC_TOL,
    Arc,
    GeometryError,
    Hermite,
    Line,
    Loop,
    chain_area,
    chain_to_loop,
    polygon,
)
from .rotor import LaminationFace, Magnet, RotorPole, build_rotor, construction_points
from .stator import Conductor, StatorPole, build_stator, winding_layout


def pole_model(p):
    """Rotor + stator for one pole (pole frame)."""
    return build_rotor(p), build_stator(p)


__all__ = [
    "BIARC_TOL", "Arc", "Conductor", "GeometryError", "Hermite", "LaminationFace", "LayerInput", "LayerPoints",
    "Line", "Loop", "Magnet", "RotorPole", "StatorPole", "build_rotor", "build_stator", "chain_area",
    "chain_to_loop", "construction_points", "fillet", "layer_points", "points_in_polygon", "pole_model",
    "polygon", "polyline_distance", "seg_seg_distance", "self_intersects", "winding_layout",
]
