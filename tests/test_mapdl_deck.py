import math
import re
from pathlib import Path

import numpy as np
import pytest

from sipm.geometry import Loop, build_rotor
from sipm.mapdl import (
    RESULTS_FILE,
    errors_from_output,
    mech_model_from_rotor,
    mech_model_ring,
    parse_results,
    write_deck,
)
from sipm.params import TOPOLOGIES, default_params

FIX = Path(__file__).parent / "fixtures"


def _deck(tmp_path, topo, **kw):
    p = default_params(topo, **kw)
    m = mech_model_from_rotor(build_rotor(p))
    return p, write_deck(p, m, tmp_path / "x.inp").read_text()


@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_deck_structure(tmp_path, topo):
    p, d = _deck(tmp_path, topo)
    lines = d.splitlines()
    assert lines[0] == "/BATCH"
    for cmd in ("/PREP7", "ET,1,PLANE183,,,0", "AMESH,ALL", "/SOLU", "OMEGA,0,0,OMG", "SOLVE", "/POST1",
                "*CFOPEN", "('done = 1')", "/EXIT,NOSAVE"):
        assert cmd in d, cmd
    omega = float(re.search(r"\*SET,OMG,(\S+)", d).group(1))
    assert math.isclose(omega, 20000 * 2 * math.pi / 60, rel_tol=1e-9)
    # symmetry on both cut lines, half pole = 22.5 deg
    assert d.count("D,ALL,UY,0") == 2 and "22.5-" in d
    # magnets as separate areas with their own contact pair, bonded by default (MPC, the
    # designed clearance bridged: initial gap ignored, absolute pinball wider than the gap)
    assert "ET,12,CONTA172" in d and p.magnet_contact == "bonded"
    for cmd in ("KEYOPT,12,2,2", "KEYOPT,12,9,1", "KEYOPT,12,12,5"):
        assert cmd in d, cmd
    pinb = float(re.search(r"R,20,,,1.0,0.1,0,-(\S+)", d).group(1))
    assert pinb > max(p.pm_side_clearance, p.pm_end_clearance)
    assert d.count("NUMSTR,AREA,20") == 2  # two half-pole magnets: areas 201, 202
    # lamination nodes dumped for the named stress sections; global peak kept
    assert "*VWRITE,NXA(1),NYA(1),NSA(1)" in d and "*GET,SVM_LAM,SORT,,MAX" in d
    assert "AATT,3,1,1,11" in d  # sleeve uses the cylindrical element system
    assert "RMODIF,10,10,DELTA" in d
    assert math.isclose(float(re.search(r"\*SET,DELTA,(\S+)", d).group(1)), p.interference)


def test_deck_frictional_magnets(tmp_path):
    """Frictional option: open-clearance contact with stabilisation (bridge-retained rotor);
    the bridgeless rotor refuses it with a clear error."""
    p, d = _deck(tmp_path, "sleeved_bridged", magnet_contact="frictional")
    assert "KEYOPT,12,9,0" in d and "KEYOPT,12,15,2" in d and "KEYOPT,12,12,5" not in d
    assert "R,20,,,1.0,0.1,-0.005" in d
    with pytest.raises(ValueError, match="not supported for sleeved_bridgeless"):
        _deck(tmp_path, "sleeved_bridgeless", magnet_contact="frictional")


def test_deck_orthotropic_mapping(tmp_path):
    p, d = _deck(tmp_path, "sleeved_bridged")
    get = lambda k: float(re.search(rf"MP,{k},3,(\S+)", d).group(1))  # noqa: E731
    assert math.isclose(get("EX"), 8.8e3)          # radial
    assert math.isclose(get("EY"), 125e3)          # hoop (fibre)
    assert math.isclose(get("PRXY"), 0.28 * 8.8 / 125)  # minor ratio in MAPDL convention
    assert math.isclose(get("DENS"), 1800e-12)
    assert math.isclose(float(re.search(r"MP,EX,1,(\S+)", d).group(1)), 200e3)


def test_deck_plane_strain_and_mesh(tmp_path):
    p, d = _deck(tmp_path, "sleeved_bridged", plane="strain", mech_mesh=0.4, mech_mesh_fine=0.1)
    assert "ET,1,PLANE183,,,2" in d
    assert "ESIZE,0.4" in d and "LESIZE,ALL,0.1" in d
    assert "AESIZE,ALL,0.4" in d  # sleeve 1.5 mm / 3 = 0.5 > 0.4


def test_deck_arcs_and_loops(tmp_path):
    _, d = _deck(tmp_path, "sleeved_bridgeless")
    for m in re.finditer(r"LARC,(\d+),(\d+),(\d+),(\S+)", d):
        assert float(m.group(4)) > 0
    # bridgeless half pole: three separate lamination pieces, the pole cap carries the d-axis hole
    assert [n for n in (1, 2, 3, 4) if f"NUMSTR,AREA,{n}\n" in d] == [1, 2, 3]
    assert d.count("ASBA,") == 0  # the d-axis hole is cut by the symmetry line: part of the cap outline


def _deck_areas(deck: str) -> dict[int, tuple[float, int | None]]:
    """Interpret the geometry commands of a deck: area number -> (area mm^2, material)."""
    kp: dict[int, np.ndarray] = {}
    lines: dict[int, tuple] = {}
    areas: dict[int, list] = {}
    comps: dict[str, set] = {}
    nline = narea = 0
    selected: set = set()
    for cmd in deck.splitlines():
        f = cmd.split(",")
        if f[0] == "K":
            kp[int(f[1])] = np.array([float(f[2]), float(f[3])])
        elif f[0] == "NUMSTR":
            if f[1] == "LINE":
                nline = int(f[2])
            else:
                narea = int(f[2])
        elif f[0] in ("L", "LARC"):
            lines[nline] = (int(f[1]), int(f[2]), int(f[3]) if f[0] == "LARC" else None)
            nline += 1
        elif f[0] == "LSEL" and f[1] == "S" and f[2] == "LINE":
            lsel = range(int(f[4]), int(f[5]) + 1)
        elif f[0] == "AL":
            segs = [lines[i] for i in lsel]
            loop = Loop(np.array([kp[a] for a, _, _ in segs]), [None if c is None else kp[c] for _, _, c in segs])
            areas[narea] = [loop.area(), None]
            narea += 1
        elif f[0] in ("ALLSEL",) or (f[0] == "ASEL" and f[1] == "ALL"):
            selected = set(areas)
        elif f[0] == "ASEL" and f[1] == "S":
            lo = int(f[4])
            hi = int(f[5]) if len(f) > 5 and f[5] else lo
            selected = {a for a in areas if lo <= a <= hi}
        elif f[0] == "CM" and f[2] == "AREA":
            comps[f[1]] = set(selected)
        elif f[0] == "ASBA":
            base, cut = int(f[1]), comps[f[2]]
            new = max(areas) + 1000
            areas[new] = [areas.pop(base)[0] - sum(areas.pop(c)[0] for c in cut), None]
        elif f[0] == "AATT":
            for a in selected:
                if areas[a][1] is None:
                    areas[a][1] = int(f[1])
    return {k: (v[0], v[1]) for k, v in areas.items()}


@pytest.mark.parametrize("topo", TOPOLOGIES)
def test_deck_areas_have_materials_and_match_geometry(tmp_path, topo):
    from sipm.mapdl import MAT_LAM, MAT_PM, MAT_SLV

    p = default_params(topo)
    rp = build_rotor(p)
    d = write_deck(p, mech_model_from_rotor(rp), tmp_path / "x.inp").read_text()
    areas = _deck_areas(d)
    assert all(mat is not None for _, mat in areas.values()), areas
    by = {}
    for a, mat in areas.values():
        by[mat] = by.get(mat, 0.0) + a
    lam = sum(f.area() for f in rp.half_lamination())
    assert math.isclose(by[MAT_LAM], lam, rel_tol=1e-9)
    assert math.isclose(by[MAT_LAM], rp.metrics["lamination_area"] / 2, rel_tol=1e-5)  # biarc curves
    assert math.isclose(by[MAT_PM], sum(m.loop.area() for m in rp.half_magnets()), rel_tol=1e-9)
    assert math.isclose(by[MAT_SLV], 0.5 * (rp.r_sleeve_out ** 2 - rp.r_lam ** 2) * rp.half_pitch, rel_tol=1e-7)


def test_ring_model_has_origin_constraint(tmp_path):
    p = default_params("sleeved_bridged")
    d = write_deck(p, mech_model_ring(71.1, 1.5), tmp_path / "r.inp").read_text()
    assert "D,ALL,UX,0" in d and "ET,12" not in d


def test_no_absolute_paths_in_deck(tmp_path):
    _, d = _deck(tmp_path, "sleeved_bridged")
    assert str(tmp_path) not in d
    assert RESULTS_FILE.split(".")[0] in d


def test_parse_results_fixture():
    r = parse_results(FIX / "mapdl_results_bridged.txt")
    assert r["done"] == 1 and r["cnvg"] == 1
    assert math.isclose(r["svm_lam"], 779.136731221, rel_tol=1e-9)
    assert math.isclose(r["sh_slv"], 301.701215353, rel_tol=1e-9)
    assert r["nnode"] == 19014


def test_errors_from_output_fixture():
    errs = errors_from_output(FIX / "mapdl_error_excerpt.out")
    assert len(errs) == 1
    assert "rigid body motion" in errs[0]
    assert errors_from_output(FIX / "does_not_exist.out") == []
