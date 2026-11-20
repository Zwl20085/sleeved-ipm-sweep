"""Load a Maxwell 2D design from an ``.aedt`` text file and replay its geometry history.

``load_design(path, design)`` parses the project once (cached), ``build(design_ctx,
overrides, naming)`` evaluates the variables (with optional overrides such as
``{"Wcf": "1mm"}``) and returns a :class:`~.aedt_geom.GeometryModel` whose
``final_body(pid)`` gives the exact outline of every object.

Topological naming: AEDT refers to edges/vertices by ID. In the saved (nominal) state the
replay checks every ID against the positions AEDT cached in the file; for another variable
set the IDs are mapped through the signatures recorded in the nominal replay
(:func:`nominal_naming`). No path is hard-coded: the caller passes the ``.aedt`` file.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .aedt_expr import Evaluator
from .aedt_geom import GeometryModel
from .aedt_parse import AedtFile, Block
from .aedt_vars import design_variables, make_evaluator, project_datasets, project_variables


@dataclass(frozen=True)
class DesignContext:
    path: str
    design: str
    block: Block
    pvars: tuple
    datasets: dict
    dvars: tuple


@lru_cache(maxsize=4)
def _project(path: str) -> tuple:
    f = AedtFile.load(path)
    return (f, project_variables(f.project_block("GlobalVariables")),
            project_datasets(f.project_block("ProjectDatasets")), f.design_blocks())


def load_design(path: str, design: str) -> DesignContext:
    f, pv, ds, designs = _project(path)
    if design not in designs:
        raise KeyError(f"design {design!r} not in {path}; available: {sorted(designs)}")
    blk = designs[design]
    return DesignContext(path, design, blk, tuple(pv), ds, tuple(design_variables(blk)))


def evaluator(ctx: DesignContext, overrides: dict | None = None) -> Evaluator:
    return make_evaluator(list(ctx.dvars), list(ctx.pvars), ctx.datasets, overrides)


def build(ctx: DesignContext, overrides: dict | None = None, naming: dict | None = None) -> GeometryModel:
    return GeometryModel(ctx.block, evaluator(ctx, overrides), naming=naming)


def nominal_naming(ctx: DesignContext, parts: list[str]) -> tuple[GeometryModel, dict]:
    """Replay the saved state, resolve mirror-clone numbering and return the ID signatures
    needed to replay other variable sets."""
    gm = build(ctx)
    gm.resolve_mirror_rules()
    for name in parts:
        gm.final_body(gm.part_by_name(name).pid)
    return gm, gm.naming


def mm(ev: Evaluator, name: str) -> float:
    """Value of a length variable in mm."""
    return ev.var(name).v * 1e3
