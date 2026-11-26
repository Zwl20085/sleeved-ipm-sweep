"""Compare a built rotor with a Maxwell reference fixture (tests/fixtures/maxwell_ref).

A fixture segment matches a package segment when both are of the same kind and their end
points (and, for arcs, centre and radius; for curves, nine points at uniform curve parameter)
agree within ``tol``. Regions are matched by area first.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from sipm.geometry import build_rotor, build_stator, chain_area
from sipm.params import default_params

FIXTURES = Path(__file__).parent / "fixtures" / "maxwell_ref"
KIND = {"L": "line", "A": "arc", "C": "hermite"}


@dataclass
class Report:
    max_dev: float = 0.0                       # worst vertex / centre / sample deviation (mm)
    max_area_rel: float = 0.0                  # worst relative region area difference
    problems: list[str] = field(default_factory=list)
    unmatched: list[tuple[str, list]] = field(default_factory=list)  # (where, fixture segment)


def load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


STATOR_FIXTURE = "stator_saved.json"


def fixture_names() -> list[str]:
    """Rotor fixtures (every JSON file except the index and the stator fixture)."""
    return sorted(p.name for p in FIXTURES.glob("*.json") if p.name not in ("index.json", STATOR_FIXTURE))


def _pkg_seg(s) -> dict:
    d = {"kind": s.kind, "p0": s.start(), "p1": s.end()}
    if s.kind == "arc":
        d["c"], d["r"] = s.c, s.r
    if s.kind == "hermite":
        d["samples"] = [s.point(k / 8) for k in range(9)]
    return d


def _seg_dev(f: list, m: dict) -> float:
    if KIND[f[0]] != m["kind"]:
        return math.inf
    dev = max(math.dist(f[1:3], m["p0"]), math.dist(f[3:5], m["p1"]))
    if f[0] == "A":
        dev = max(dev, math.dist(f[5:7], m["c"]), abs(f[7] - m["r"]))
    if f[0] == "C":
        dev = max(dev, max(math.dist(a, b) for a, b in zip(f[5], m["samples"])))
    return dev


def compare_loop(fl: list, ml, where: str, rep: Report, tol: float, skip=None) -> None:
    ms = [_pkg_seg(s) for s in ml]
    if len(fl) != len(ms) and skip is None:
        rep.problems.append(f"{where}: {len(fl)} fixture segments vs {len(ms)} package segments")
    for f in fl:
        if skip is not None and skip(f):
            continue
        d = min(_seg_dev(f, m) for m in ms)
        if d > tol:
            rep.unmatched.append((where, f))
        else:
            rep.max_dev = max(rep.max_dev, d)


def _loop_key(seglist) -> tuple[float, float]:
    """Mean of the segment start points of a fixture loop."""
    xs = [s[1] for s in seglist]
    ys = [s[2] for s in seglist]
    return (sum(xs) / len(xs), sum(ys) / len(ys))


def _chain_key(chain) -> tuple[float, float]:
    """Mean of the segment start points of a package chain."""
    pts = [s.start() for s in chain]
    return (sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts))


def compare_face(ff: dict, outer, holes, where: str, rep: Report, tol: float, skip=None) -> None:
    compare_loop(ff["outer"], outer, where + "/outer", rep, tol, skip)
    if len(ff["holes"]) != len(holes):
        rep.problems.append(f"{where}: {len(ff['holes'])} fixture holes vs {len(holes)} package holes")
    for i, fh in enumerate(ff["holes"]):
        c = _loop_key(fh)
        best = min(holes, key=lambda h: math.dist(c, _chain_key(h)))
        compare_loop(fh, best, f"{where}/hole{i}", rep, tol, skip)


def compare_fixture(doc: dict, tol: float = 1e-6, skip=None, rp=None) -> tuple[Report, object]:
    """Build the package rotor for the fixture's parameters and compare every region."""
    rp = rp or build_rotor(default_params(doc["topology"], **doc["params"]))
    rep = Report()
    reg = doc["regions"]
    faces = sorted(rp.lamination, key=lambda f: -f.area())
    ffaces = sorted(reg["lamination"], key=lambda f: -f["area"])
    if len(faces) != len(ffaces):
        rep.problems.append(f"lamination: {len(ffaces)} fixture faces vs {len(faces)} package faces")
    for i, (ff, mf) in enumerate(zip(ffaces, faces)):
        rep.max_area_rel = max(rep.max_area_rel, abs(mf.area() - ff["area"]) / ff["area"])
        compare_face(ff, mf.outer, list(mf.holes), f"lamination{i}", rep, tol, skip)
    for i, fm in enumerate(reg["magnets"]):
        c = _loop_key(fm["outer"])
        mg = min(rp.magnets, key=lambda m: math.dist(c, _chain_key(m.chain)))
        rep.max_area_rel = max(rep.max_area_rel, abs(chain_area(mg.chain) - fm["area"]) / fm["area"])
        compare_face(fm, mg.chain, [], f"magnet{i}", rep, tol)
    if "sleeve" in reg:
        from sipm.geometry.primitives import annular_sector

        hp = rp.half_pitch
        ch = annular_sector(rp.r_lam, rp.r_sleeve_out, -hp, hp)
        rep.max_area_rel = max(rep.max_area_rel, abs(chain_area(ch) - reg["sleeve"]["area"]) / reg["sleeve"]["area"])
        compare_face(reg["sleeve"], ch, [], "sleeve", rep, tol)
    return rep, rp


def compare_stator(doc: dict, tol: float = 1e-6) -> Report:
    p = default_params(doc["topology"], **doc["params"])
    sp = build_stator(p)
    rep = Report()
    rep.max_area_rel = abs(chain_area(sp.chain) - doc["stator"]["area"]) / doc["stator"]["area"]
    compare_face(doc["stator"], sp.chain, [], "stator", rep, tol)
    for i, fc in enumerate(doc["coils"]):
        c = _loop_key(fc["outer"])
        cd = min(sp.conductors, key=lambda q: math.dist(c, _chain_key(q.chain)))
        compare_face(fc, cd.chain, [], f"coil{i}", rep, tol)
    return rep
