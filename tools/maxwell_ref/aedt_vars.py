"""Variable tables (design + project) and datasets from parsed AEDT blocks."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .aedt_expr import Dataset, Evaluator, ExprError, dependency_order, dim_str
from .aedt_parse import Block


@dataclass(frozen=True)
class VarDef:
    name: str
    flags: str
    expr: str
    line: int
    extra: tuple  # remaining VariableProp args (ranges etc.) as raw values


def variable_props(block: Block) -> list[VarDef]:
    out = []
    for f in block.fields("VariableProp"):
        args = f.value()
        if len(args) < 4 or not isinstance(args[3], str):
            raise ValueError(f"line {f.line}: unexpected VariableProp {f.raw[:120]!r}")
        out.append(VarDef(name=args[0], flags=args[1], expr=args[3], line=f.line, extra=tuple(args[4:])))
    return out


def project_variables(globals_block: Block) -> list[VarDef]:
    return variable_props(globals_block)


def project_datasets(ds_block: Optional[Block]) -> dict[str, Dataset]:
    out: dict[str, Dataset] = {}
    if ds_block is None:
        return out
    defs = ds_block.opt_block("DatasetDefinitions")
    if defs is None:
        return out
    for b in defs.blocks():
        pts = b.block("Coordinates")["Points"]
        if len(pts) % 2:
            raise ValueError(f"dataset {b.name} has odd number of values")
        xs = tuple(float(v) for v in pts[0::2])
        ys = tuple(float(v) for v in pts[1::2])
        out[b.name] = Dataset(b.name, xs, ys)
    return out


def design_variables(design: Block) -> list[VarDef]:
    return variable_props(design.block("ModelSetup").block("Properties"))


def make_evaluator(dvars: list[VarDef], pvars: list[VarDef], datasets: dict[str, Dataset],
                   overrides: Optional[dict[str, str]] = None) -> Evaluator:
    return Evaluator(design_vars={v.name: v.expr for v in dvars},
                     project_vars={v.name: v.expr for v in pvars},
                     datasets=datasets, overrides=dict(overrides or {}))


def variable_table(ev: Evaluator, dvars: list[VarDef], pvars: list[VarDef]) -> tuple[list[dict], list[dict]]:
    """Evaluate every variable; returns (rows, failures). Never raises for a single variable."""
    rows, fails = [], []
    table = {v.name: v.expr for v in pvars}
    table.update({v.name: v.expr for v in dvars})
    try:
        order = dependency_order(table, set())
    except ExprError as e:
        fails.append({"name": "<dependency order>", "error": str(e)})
        order = list(table)
    pos = {n: i for i, n in enumerate(order)}
    for scope, vs in (("project", pvars), ("design", dvars)):
        for v in vs:
            row: dict[str, Any] = {"name": v.name, "scope": scope, "flags": v.flags, "expression": v.expr,
                                   "source_line": v.line, "dependency_rank": pos.get(v.name)}
            if v.name in ev.overrides:
                row["override_expression"] = ev.overrides[v.name]
            try:
                q = ev.var(v.name)
                row["si_value"] = q.v
                row["si_unit"] = dim_str(q.dim)
            except ExprError as e:
                row["error"] = str(e)
                fails.append({"name": v.name, "expression": v.expr, "error": str(e)})
            rows.append(row)
    return rows, fails
