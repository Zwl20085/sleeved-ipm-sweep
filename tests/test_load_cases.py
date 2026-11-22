"""Load-case set (mech_cases): parameters, deck, result columns, envelope, lift-off speed."""

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sipm import loadcases as lc
from sipm.analytic import coupled_pressure, disc_ring, isotropic, liftoff_speed, orthotropic
from sipm.config import Config
from sipm.geometry import build_rotor
from sipm.mapdl import (
    CONTACT_COLUMNS,
    CONTACT_FILE,
    mech_model_from_rotor,
    mech_model_ring,
    parse_contact_elements,
    write_deck,
)
from sipm.maxwell import MECH_ONLY, em_model_key
from sipm.params import CASE_SEPARATOR, MECH_CASES, coerce, default_params, parse_mech_cases
from sipm.sweep import KEY_OMITTED_AT_DEFAULT, Sweep, point_key, resolve

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

NEW_FIELDS = ("mech_cases", "overspeed_factor", "mech_temperature_cold")
ALL = ",".join(MECH_CASES)


# ---- parameters -----------------------------------------------------------------------------
def test_defaults_single_case():
    p = default_params("sleeved_bridged")
    assert (p.mech_cases, p.overspeed_factor, p.mech_temperature_cold) == ("", 1.2, 20.0)
    assert p.mech_case_list == ()


def test_parse_cases():
    assert parse_mech_cases("") == () and parse_mech_cases("  ") == ()
    assert parse_mech_cases(" overspeed_hot , standstill_cold") == ("overspeed_hot", "standstill_cold")
    assert parse_mech_cases(ALL) == MECH_CASES


@pytest.mark.parametrize("bad", ["rated", "rated_hot,rated_warm", "rated_hot,,standstill_cold", "rated_hot,rated_hot",
                                 "rated_hot;standstill_cold"])
def test_parse_cases_rejects(bad):
    with pytest.raises(ValueError, match="mech_cases"):
        parse_mech_cases(bad)
    with pytest.raises(ValueError, match="mech_cases"):
        default_params("sleeved_bridged", mech_cases=bad)


@pytest.mark.parametrize("bad", [dict(overspeed_factor=0.9), dict(overspeed_factor=float("nan")),
                                 dict(mech_temperature_cold=-300.0), dict(mech_temperature_cold=501.0),
                                 dict(mech_cases=["rated_hot"])])
def test_invalid_rejected(bad):
    with pytest.raises(ValueError):
        default_params("sleeved_bridged", **bad)


def test_case_speed_temperature_and_for_case():
    p = default_params("sleeved_bridged", mech_rpm=18000.0, overspeed_factor=1.25, mech_temperature=140.0,
                       mech_temperature_cold=-20.0, mech_cases=ALL)
    expect = {"standstill_cold": (0.0, -20.0), "standstill_hot": (0.0, 140.0), "rated_cold": (18000.0, -20.0),
              "rated_hot": (18000.0, 140.0), "overspeed_cold": (22500.0, -20.0), "overspeed_hot": (22500.0, 140.0)}
    for case, (rpm, t) in expect.items():
        assert p.case_speed_temperature(case) == (rpm, t)
        q = p.for_case(case)
        assert (q.mech_rpm, q.mech_temperature, q.mech_cases) == (rpm, t, "")
        assert q.mech_ref_temperature == p.mech_ref_temperature and q.interference == p.interference
    with pytest.raises(ValueError):
        p.case_speed_temperature("idle_hot")


def test_coerce_and_toml_array():
    assert coerce("mech_cases", ["rated_hot", "overspeed_hot"]) == "rated_hot,overspeed_hot"
    assert coerce("overspeed_factor", 1) == 1.0 and coerce("mech_temperature_cold", "-10") == -10.0
    p = resolve({}, base={"mech_cases": ["standstill_cold", "overspeed_hot"]})
    assert p.mech_case_list == ("standstill_cold", "overspeed_hot")


def test_example_toml_resolves():
    from sipm.cli import sweep_from_config

    cfg = tomllib.loads((Path(__file__).parents[1] / "examples/configs/mech_load_cases.toml").read_text())
    plan = sweep_from_config(cfg).plan()
    assert len(plan) == 9 and all(err is None for *_, err in plan)
    p = plan[0][2]
    assert p.mech_case_list == MECH_CASES and p.mech_temperature == 120.0 and p.mech_temperature_cold == 20.0


# ---- keys -----------------------------------------------------------------------------------
def test_new_fields_in_key_lists():
    assert set(NEW_FIELDS) <= set(KEY_OMITTED_AT_DEFAULT) and set(NEW_FIELDS) <= MECH_ONLY
    p = default_params("sleeved_bridged")
    q = p.replace(mech_cases=ALL, overspeed_factor=1.3, mech_temperature_cold=-40.0)
    assert em_model_key(p) == em_model_key(q)
    base = point_key(p, "mech")
    for k, v in dict(mech_cases="rated_hot", overspeed_factor=1.5, mech_temperature_cold=0.0).items():
        assert point_key(p.replace(**{k: v}), "mech") != base
    # back at the default, the key is the pinned one again (tests/test_composite_sleeve.py)
    assert point_key(q.replace(mech_cases="", overspeed_factor=1.2, mech_temperature_cold=20.0), "mech") == base


# ---- deck -----------------------------------------------------------------------------------
# sha256 of the deck text ("\n" newlines) written by the code before mech_cases existed (HEAD of the
# composite-sleeve work, generated from `git archive HEAD src`). A deliberate change of the deck or of
# the geometry has to update these.
PINNED_DECKS = {
    ("sleeved_bridged", "iso"): "0dbe06c6d82f923ce5b554ae71bc9c0f1be707e23066a045eb49c26e912ea9e5",
    ("sleeved_bridged", "hot"): "85ed5b91f36a52d230608b593d6de4832c6c953a0a1e5854dc7e7e91ba9ff90d",
    # bridgeless: the thin-sleeve bore arc sized like the lamination OD (one added LESIZE line; the decks
    # before it hashed ca6115e7... (iso) and 06241907... (hot))
    ("sleeved_bridgeless", "iso"): "16dbc235f8b14b3290b362ed21be512c61a4577929b007764ceab5ad78d2f5b6",
    ("sleeved_bridgeless", "hot"): "5ebb8ed61adbefe579a4fb75b8e2f5d700d78cc1ff5a4395e51189d121921a40",
    ("ring", "iso"): "45d2736261b67693dbdeca57bc502f22e1c80fd3f2886115fcefa6cb9a8da92e",
    ("ring", "hot"): "47c4ea4657ec82fea1c2baf0ac120dddc24897e1c55dbcb837bd7d765f90b87b",
}


def _deck_text(tmp_path, topo, p, **kw):
    model = mech_model_ring(71.1, 1.5) if topo == "ring" else mech_model_from_rotor(build_rotor(p))
    return write_deck(p, model, tmp_path / "d.inp", **kw).read_text(encoding="ascii")


@pytest.mark.parametrize("key", list(PINNED_DECKS))
def test_deck_unchanged_without_cases(tmp_path, key):
    topo, tag = key
    kw = dict(mech_temperature=120.0) if tag == "hot" else {}
    p = default_params("sleeved_bridged" if topo == "ring" else topo, **kw)
    text = _deck_text(tmp_path, topo, p)
    assert hashlib.sha256(text.encode()).hexdigest() == PINNED_DECKS[key]
    assert "CEMX" not in text and CONTACT_FILE.split(".")[0] not in text
    # mech_cases does not enter the deck itself (one deck per case, see Params.for_case)
    assert _deck_text(tmp_path, topo, p.replace(mech_cases=ALL, overspeed_factor=1.5)) == text


@pytest.mark.parametrize("topo, kw, sized", [
    ("sleeved_bridgeless", {}, True),                                  # 0.5 mm sleeve: area size 0.167 < 0.5
    ("sleeved_bridgeless", dict(sleeve_thickness=1.5), False),         # area size = mech_mesh: deck as before
    ("sleeved_bridgeless", dict(sleeve_thickness=0.9, mech_mesh=0.8), True),
    ("sleeved_bridged", {}, False),                                    # 1.5 mm: equal sizes, default deck pinned
    ("sleeved_bridged", dict(sleeve_thickness=1.0), True),             # thin bridge-retained sleeve: same rule
    ("sleeved_bridged", dict(mech_mesh=0.8, mech_mesh_fine=0.3), True),   # coarse mesh: 0.5 < 0.8
    ("sleeved_bridged", dict(sleeve_thickness=2.0), False),
    ("sleeved_bridged", dict(mech_mesh=0.3, mech_mesh_fine=0.08), False),
    ("ring", {}, False),                                               # 1.5 mm ring of the closed-form check
])
def test_bore_arc_sized_like_lamination(tmp_path, topo, kw, sized):
    from sipm.mapdl import bore_arc_sized

    p = default_params("sleeved_bridged" if topo == "ring" else topo, **kw)
    lines = _deck_text(tmp_path, topo, p).splitlines()
    bore = [ln for ln in lines if ln.startswith("LESIZE,30003")]
    assert bore == ([f"LESIZE,30003,{p.mech_mesh:.10g},,,,1"] if sized else [])
    model = mech_model_ring(71.1, 1.5) if topo == "ring" else mech_model_from_rotor(build_rotor(p))
    assert bore_arc_sized(p, model) is sized
    if sized:   # after the area size of the sleeve (which it overrides on that line) and before meshing
        aesize = "AESIZE,ALL,%.10g" % min(p.mech_mesh, p.sleeve_thickness / 3)
        assert lines.index(aesize) < lines.index(bore[0]) < lines.index("AMESH,ALL")


def test_bore_arc_rule_adds_one_line_only(tmp_path):
    """On a thin bridge-retained sleeve the rule adds the one LESIZE line and changes nothing else."""
    p = default_params("sleeved_bridged", sleeve_thickness=1.0)
    lines = _deck_text(tmp_path, "sleeved_bridged", p).splitlines()
    i = lines.index("LESIZE,30003,0.5,,,,1")
    assert lines[i - 1] == "ALLSEL" and lines[i - 2] == "AESIZE,ALL,%.10g" % (1.0 / 3)
    # layered sleeves size their bore arc themselves and are not affected
    from sipm.mapdl import bore_arc_sized

    q = default_params("sleeved_bridged", sleeve_thickness=1.0, sleeve_fibre="carbon_sm", sleeve_layers="0:1,30:1")
    assert not bore_arc_sized(q, mech_model_from_rotor(build_rotor(q)))


@pytest.mark.parametrize("topo", ["sleeved_bridged", "sleeved_bridgeless", "ring"])
def test_deck_contact_element_export(tmp_path, topo):
    p = default_params("sleeved_bridged" if topo == "ring" else topo)
    plain = _deck_text(tmp_path, topo, p).splitlines()
    lines = _deck_text(tmp_path, topo, p, contact_elements=True).splitlines()
    i = lines.index("*VWRITE,CEID(1),CEX(1),CEY(1),CEP(1),CES(1)")
    assert lines[i + 1] == "(5E21.12)" and lines[i + 2] == "*CFCLOS" and lines[i - 1] == "*VMASK,CMSK(1)"
    stem, ext = CONTACT_FILE.split(".")
    assert lines[i - 2] == f"*CFOPEN,{stem},{ext}"
    # written while the sleeve contact elements (type 10) are selected, after the pressure and
    # status tables exist and before the open elements are reselected
    j = lines.index("*GET,CEMX,ELEM,0,NUM,MAXD")
    assert lines.index("ESEL,S,TYPE,,10") < lines.index("ETABLE,CSTA,CONT,STAT") < j < i
    assert i < lines.index("ESEL,R,ETAB,CSTA,-1,1.5")
    for item in ("CENT,X", "CENT,Y", "ETAB,CPRS", "ETAB,CSTA"):
        assert any(ln.startswith("*VGET,") and ln.endswith(f",ELEM,1,{item}") for ln in lines), item
    # nothing else changes
    block = set(range(j, i + 3))
    assert [ln for k, ln in enumerate(lines) if k not in block] == plain


def test_parse_contact_elements(tmp_path):
    f = tmp_path / CONTACT_FILE
    f.write_text("   0.961200000000E+04   0.71E+02   0.25E+00   0.3148E+01   0.2E+01\n"
                 "   0.961100000000E+04   0.65E+02   0.27E+02   0.3150E+01   0.3E+01\n")
    e = parse_contact_elements(f)
    assert tuple(e) == CONTACT_COLUMNS
    assert e["elem"].tolist() == [9611, 9612] and e["elem"].dtype.kind == "i"   # sorted by element
    assert e["pres"].tolist() == [3.150, 3.148] and e["status"].tolist() == [3.0, 2.0]
    f.write_text("1 2 3\n")
    with pytest.raises(ValueError, match="columns"):
        parse_contact_elements(f)


# ---- columns and envelope -------------------------------------------------------------------
def _fake(case_index, **over):
    """Synthetic single-case result: values that differ per case."""
    i = case_index
    r = dict(svm_lam_max=200.0 + 10 * i, svm_rotor_max=200.0 + 10 * i, svm_lam_where="qhole",
             sh_sleeve_max=150 + 20 * i, sr_sleeve_max=0.1 * i, sleeve_fi_max_stress=0.1 + 0.01 * i,
             sleeve_fi_tsai_wu=0.2 + 0.01 * i,
             sleeve_fi_hashin_fibre=0.05 * i, sleeve_fi_hashin_matrix=0.3 - 0.01 * i, sleeve_sf=5.0 - 0.1 * i,
             sleeve_open_fraction=0.0, contact_p_min=4.0 - 0.5 * i, n_nodes=100, mech_wall_s=1.0,
             _raw={"x": i}, _sleeve_nodes={"r": np.zeros(2)})
    r.update(over)
    return r


def test_column_names():
    assert CASE_SEPARATOR == "__"
    assert lc.case_column("sh_sleeve_max", "overspeed_hot") == "sh_sleeve_max__overspeed_hot"
    assert lc.worst_column("sleeve_sf") == "sleeve_sf__worst"
    assert lc.worst_case_column("sleeve_sf") == "sleeve_sf__worst_case"
    # every per-case column is a valid identifier (MATLAB readtable keeps it as is)
    for out in ("sleeve_fi_hashin_matrix", "svm_lam_max_all"):
        for case in MECH_CASES:
            c = lc.case_column(out, case)
            assert c.isidentifier() and len(c) <= 63


def test_base_case():
    assert lc.base_case(("standstill_cold", "rated_hot", "overspeed_hot")) == "rated_hot"
    assert lc.base_case(("overspeed_hot", "standstill_cold")) == "overspeed_hot"
    with pytest.raises(ValueError):
        lc.base_case(())


def test_envelope_worst_and_governing_case():
    res = {c: _fake(i) for i, c in enumerate(MECH_CASES)}
    env = lc.envelope(res)
    assert env["svm_lam_max__worst"] == 250.0 and env["svm_lam_max__worst_case"] == "overspeed_hot"
    assert env["sleeve_sf__worst"] == pytest.approx(4.5) and env["sleeve_sf__worst_case"] == "overspeed_hot"
    assert env["contact_p_min__worst"] == 1.5 and env["contact_p_min__worst_case"] == "overspeed_hot"
    assert env["sleeve_fi_hashin_matrix__worst"] == 0.3
    assert env["sleeve_fi_hashin_matrix__worst_case"] == "standstill_cold"
    # all equal -> the first case in canonical order governs
    assert env["sleeve_open_fraction__worst"] == 0.0 and env["sleeve_open_fraction__worst_case"] == "standstill_cold"
    # every output of a single-layer point; the interlaminar shear exists only for a layered sleeve
    # (tests/test_multilayer_sleeve.py) and is left out like any missing output
    assert set(env) == {f"{n}__worst{s}" for n, _ in lc.ENVELOPE for s in ("", "_case")
                        if n != "tau_sleeve_interface_max"}
    layered = lc.envelope({c: _fake(i, tau_sleeve_interface_max=0.1 * (3 - i)) for i, c in enumerate(MECH_CASES)})
    assert layered["tau_sleeve_interface_max__worst"] == pytest.approx(0.3)
    assert layered["tau_sleeve_interface_max__worst_case"] == "standstill_cold"


def test_envelope_ties_nan_and_missing():
    res = {"overspeed_hot": _fake(0, svm_lam_max=300.0, contact_p_min=float("nan")),
           "standstill_cold": _fake(1, svm_lam_max=300.0, contact_p_min=2.0)}
    del res["overspeed_hot"]["sr_sleeve_max"], res["standstill_cold"]["sr_sleeve_max"]
    env = lc.envelope(res)
    assert env["svm_lam_max__worst_case"] == "standstill_cold"     # tie: canonical order, not dict order
    assert env["contact_p_min__worst"] == 2.0                       # NaN skipped
    assert "sr_sleeve_max__worst" not in env                         # absent everywhere: no column


def test_combine_order_independent_and_base():
    res = {c: _fake(i) for i, c in enumerate(MECH_CASES)}
    p1 = default_params("sleeved_bridged", mech_cases=ALL)
    p2 = p1.replace(mech_cases=",".join(reversed(MECH_CASES)))
    a, b = lc.combine(p1, res, wall_s=7.0), lc.combine(p2, dict(reversed(list(res.items()))), wall_s=7.0)
    pub = lambda d: {k: v for k, v in d.items() if not k.startswith("_")}  # noqa: E731
    assert pub(a) == pub(b) and list(pub(a)) == list(pub(b))   # same values and the same column order
    assert a["mech_case"] == "rated_hot" and a["svm_lam_max"] == res["rated_hot"]["svm_lam_max"]
    assert a["mech_wall_s"] == 7.0 and a["_raw"] == res["rated_hot"]["_raw"] and set(a["_cases"]) == set(MECH_CASES)
    for case in MECH_CASES:
        assert a[f"sh_sleeve_max__{case}"] == res[case]["sh_sleeve_max"]
        assert a[f"svm_lam_where__{case}"] == "qhole"
        assert f"_raw__{case}" not in a
    # without rated_hot the first listed case fills the un-suffixed columns
    p3 = p1.replace(mech_cases="overspeed_cold,standstill_hot")
    c = lc.combine(p3, {k: res[k] for k in ("overspeed_cold", "standstill_hot")})
    assert c["mech_case"] == "overspeed_cold" and c["contact_p_min"] == res["overspeed_cold"]["contact_p_min"]
    assert "sh_sleeve_max__rated_hot" not in c and "liftoff_rpm_cold" not in c   # no contact arrays given
    with pytest.raises(ValueError, match="missing"):
        lc.combine(p1, {"rated_hot": res["rated_hot"]})


# ---- lift-off -------------------------------------------------------------------------------
def _elems(p0, c, s, status=None, n0=500):
    p0, c = np.asarray(p0, float), np.asarray(c, float)
    n = len(p0)
    th = np.radians(np.linspace(0.2, 22.3, n))
    pres = np.maximum(p0 - c * s * s, 0.0)
    st = np.where(pres > 0, 2.0, 0.0) if status is None else np.full(n, status, float)
    return dict(elem=np.arange(n0, n0 + n), x=71.1 * np.cos(th), y=71.1 * np.sin(th), pres=pres, status=st)


# element 0 has the lowest pressure at standstill, element 2 falls fastest and opens first
P0, C = [3.0, 4.0, 5.0, 6.0], [1.0e-9, 1.5e-9, 4.0e-9, -1e-10]


def test_liftoff_from_elements_governing_element_changes():
    s, el = lc.liftoff_from_elements(0.0, _elems(P0, C, 0.0), 20000.0, _elems(P0, C, 20000.0))
    assert el == 502 and s == pytest.approx(math.sqrt(5.0 / 4.0e-9), rel=1e-12)
    # the minimum of the pressures at each speed would extrapolate to the wrong speed
    pmin0, pmin1 = min(P0), min(np.array(P0) - np.array(C) * 20000.0 ** 2)
    naive = math.sqrt(pmin0 / ((pmin0 - pmin1) / 20000.0 ** 2))
    assert abs(naive - s) / s > 0.05
    # any two speeds of the linear range give the same answer, in either order
    s2, _ = lc.liftoff_from_elements(25000.0, _elems(P0, C, 25000.0), 10000.0, _elems(P0, C, 10000.0))
    assert s2 == pytest.approx(s, rel=1e-9)


def test_liftoff_from_elements_edge_cases():
    rising = [-1e-9, -2e-9]
    a, b = _elems([1, 2], rising, 0.0), _elems([1, 2], rising, 1e4)
    assert lc.liftoff_from_elements(0.0, a, 1e4, b) == (math.inf, -1)
    with pytest.raises(ValueError, match="different speeds"):
        lc.liftoff_from_elements(1e4, _elems(P0, C, 1e4), 1e4, _elems(P0, C, 1e4))
    a, b = _elems(P0, C, 0.0), _elems(P0, C, 1e4)
    b["x"] = b["x"] + 1e-3
    with pytest.raises(ValueError, match="centroids"):
        lc.liftoff_from_elements(0.0, a, 1e4, b)
    with pytest.raises(ValueError, match="different contact elements"):
        lc.liftoff_from_elements(0.0, a, 1e4, _elems(P0, C, 1e4, n0=1))


def test_liftoff_at_temperature_rules():
    exact = math.sqrt(5.0 / 4.0e-9)   # 35355 rpm
    ent = [("rated_cold", 20000.0, _elems(P0, C, 20000.0)), ("standstill_cold", 0.0, _elems(P0, C, 0.0)),
           ("overspeed_cold", 24000.0, _elems(P0, C, 24000.0))]
    rpm, basis = lc.liftoff_at_temperature(ent)
    assert rpm == pytest.approx(exact, rel=1e-9) and basis == "standstill_cold+overspeed_cold"
    # overspeed beyond lift-off: its contact is partly open, the closed pair is used
    big = [("standstill_cold", 0.0, _elems(P0, C, 0.0)), ("rated_cold", 20000.0, _elems(P0, C, 20000.0)),
           ("overspeed_cold", 40000.0, _elems(P0, C, 40000.0))]
    rpm, basis = lc.liftoff_at_temperature(big)
    assert rpm == pytest.approx(exact, rel=1e-9) and basis == "standstill_cold+rated_cold"
    # open at standstill -> 0; open at the lowest (non-zero) speed -> NaN
    shut = [("standstill_hot", 0.0, _elems([0.0, 1.0], [1e-9, 1e-9], 0.0)),
            ("rated_hot", 2e4, _elems([1.0, 1.0], [1e-9, 1e-9], 2e4))]
    assert lc.liftoff_at_temperature(shut) == (0.0, "open at standstill_hot")
    rpm, basis = lc.liftoff_at_temperature([("rated_hot", 40000.0, _elems(P0, C, 40000.0)),
                                            ("overspeed_hot", 48000.0, _elems(P0, C, 48000.0))])
    assert math.isnan(rpm) and "lowest" in basis
    # one closed case only, or none at all
    rpm, basis = lc.liftoff_at_temperature([("standstill_cold", 0.0, _elems(P0, C, 0.0)),
                                            ("overspeed_cold", 40000.0, _elems(P0, C, 40000.0))])
    assert math.isnan(rpm) and "two" in basis
    rpm, basis = lc.liftoff_at_temperature([])
    assert math.isnan(rpm) and "no case" in basis
    # a case open below the extrapolated speed is reported, not hidden
    odd = [("standstill_cold", 0.0, _elems(P0, C, 0.0)), ("rated_cold", 20000.0, _elems(P0, C, 20000.0, status=0.0)),
           ("overspeed_cold", 24000.0, _elems(P0, C, 24000.0))]
    rpm, basis = lc.liftoff_at_temperature(odd)
    assert "inconsistent: rated_cold" in basis
    # the same speed twice (hot = cold temperature) counts once
    dup = ent + [("standstill_hot", 0.0, _elems(P0, C, 0.0))]
    assert lc.liftoff_at_temperature(dup)[0] == pytest.approx(exact, rel=1e-9)


@pytest.mark.parametrize("bump, flagged", [(0.001, False), (0.009, False), (0.011, True), (0.087, True)])
def test_liftoff_three_point_linearity_flag(bump, flagged):
    """A middle closed case off the two-point line by more than LINEARITY_TOL (1 % of the element's
    standstill pressure) is flagged; the value is still the two-point extrapolation."""
    assert lc.LINEARITY_TOL == 1e-2
    mid = _elems(P0, C, 20000.0)
    mid["pres"] = mid["pres"] + bump * np.array(P0) * (np.arange(len(P0)) == 1)   # one element off its line
    ent = [("standstill_cold", 0.0, _elems(P0, C, 0.0)), ("rated_cold", 20000.0, mid),
           ("overspeed_cold", 24000.0, _elems(P0, C, 24000.0))]
    assert lc.linearity_deviation(0.0, ent[0][2], 24000.0, ent[2][2], 20000.0, mid) == pytest.approx(bump, rel=1e-6)
    rpm, basis = lc.liftoff_at_temperature(ent)
    assert rpm == pytest.approx(math.sqrt(5.0 / 4.0e-9), rel=1e-9)
    assert ("not linear: rated_cold deviates" in basis) is flagged
    assert basis.startswith("standstill_cold+overspeed_cold")


def test_liftoff_bridgeless_not_computed():
    p = default_params("sleeved_bridgeless", mech_cases="standstill_cold,rated_cold")
    res = {"standstill_cold": {"_contact_elems": _elems(P0, C, 0.0)},
           "rated_cold": {"_contact_elems": _elems(P0, C, 20000.0)}}
    out = lc.liftoff_outputs(p, res)
    for label in ("cold", "hot"):
        assert math.isnan(out[f"liftoff_rpm_{label}"])
        assert "sleeved_bridgeless" in out[f"liftoff_basis_{label}"]
        assert "docs/validation.md section 11" in out[f"liftoff_basis_{label}"]
    assert lc.nonlinear_reason(default_params("sleeved_bridged")) == ""


def test_liftoff_outputs_temperatures_and_linearity():
    p = default_params("sleeved_bridged", mech_rpm=20000.0, mech_temperature=120.0, mech_cases=ALL)
    hot_p0 = [x + 2.0 for x in P0]
    res = {}
    for case in MECH_CASES:
        rpm, t = p.case_speed_temperature(case)
        res[case] = {"_contact_elems": _elems(hot_p0 if t == 120.0 else P0, C, rpm)}
    out = lc.liftoff_outputs(p, res)
    assert out["liftoff_rpm_cold"] == pytest.approx(math.sqrt(5.0 / 4.0e-9), rel=1e-9)
    assert out["liftoff_rpm_hot"] == pytest.approx(math.sqrt(7.0 / 4.0e-9), rel=1e-9)
    assert out["liftoff_basis_hot"] == "standstill_hot+overspeed_hot"
    # equal temperatures: both labels use all cases at that temperature
    q = p.replace(mech_temperature=20.0, mech_cases="standstill_cold,rated_hot")
    res2 = {"standstill_cold": {"_contact_elems": _elems(P0, C, 0.0)},
            "rated_hot": {"_contact_elems": _elems(P0, C, 20000.0)}}
    out2 = lc.liftoff_outputs(q, res2)
    assert out2["liftoff_rpm_cold"] == out2["liftoff_rpm_hot"] == pytest.approx(math.sqrt(5.0 / 4.0e-9), rel=1e-9)
    # friction or loose magnets: not linear, no value
    for kw in (dict(mu_sleeve=0.1), dict(magnet_contact="frictional")):
        o = lc.liftoff_outputs(p.replace(**kw), res)
        assert math.isnan(o["liftoff_rpm_cold"]) and "not linear" in o["liftoff_basis_cold"]


# ---- lift-off by search (contacts not linear in speed) --------------------------------------
class Contact:
    """Synthetic contact: element pressures p0_e - c_e s^2 - d_e s^4 (d = 0: the linear model), an element
    open (status 0) where its pressure is <= 0. ``solve`` counts the calls like a MAPDL run would."""

    def __init__(self, p0=P0, c=C, d=(0.0, 0.0, 0.0, 0.0), f=(0.0, 0.0, 0.0, 0.0)):
        self.p0, self.c, self.d, self.f = (np.asarray(v, float) for v in (p0, c, d, f))
        self.calls = []

    def pres(self, s):
        return self.p0 - self.c * s * s - self.d * s ** 4 - self.f * s ** 6

    def elems(self, s):
        pres = self.pres(s)
        el = _elems(np.maximum(pres, 0.0), np.zeros(len(pres)), 0.0)
        el["pres"], el["status"] = np.maximum(pres, 0.0), np.where(pres > 0, 2.0, 0.0)
        return el

    def solve(self, s):
        self.calls.append(s)
        return self.elems(s)

    def first_opening(self):
        """Lowest speed at which an element's pressure reaches 0 (bisection to 1e-9 rpm)."""
        lo, hi = 0.0, 1e6
        for _ in range(200):
            m = 0.5 * (lo + hi)
            lo, hi = (m, hi) if np.all(self.pres(m) > 0) else (lo, m)
        return hi

    def cases(self, speeds=(("standstill_cold", 0.0), ("rated_cold", 20000.0), ("overspeed_cold", 24000.0))):
        return [(n, s, self.elems(s)) for n, s in speeds]


QUARTIC = (0.0, 0.0, 2.0e-18, 0.0)   # element 2 opens earlier than its two-point line predicts
SEXTIC = (0.0, 0.0, 1.0e-27, 0.0)    # ... and is not quadratic in s^2 either


def test_quadratic_zero_exact_for_quadratic_in_speed_squared():
    k = Contact(d=QUARTIC)
    pts = [(s, k.elems(s)) for s in (0.0, 20000.0, 24000.0)]
    assert lc._quadratic_zero(pts) == pytest.approx(k.first_opening(), rel=1e-9)
    lin = Contact()
    assert lc._quadratic_zero([(s, lin.elems(s)) for s in (0.0, 20000.0, 24000.0)]) == pytest.approx(
        math.sqrt(5.0 / 4.0e-9), rel=1e-9)
    never = Contact(c=[-1e-9] * 4)
    assert lc._quadratic_zero([(s, never.elems(s)) for s in (0.0, 1e4, 2e4)]) == math.inf


@pytest.mark.parametrize("tol", [100.0, 10.0])
def test_liftoff_search_brackets_above_the_cases(tol):
    """All cases closed: first trial just beyond the extrapolation through the three fastest closed
    cases (quadratic in s^2 per element), then further estimates; the result lies within tol / 2 of the
    true lift-off speed, and the basis states the iteration and its solve count."""
    k = Contact(d=QUARTIC, f=SEXTIC)
    exact = k.first_opening()
    rpm, basis, trials = lc.liftoff_search(k.cases(), k.solve, tol)
    est = lc._quadratic_zero([(s, k.elems(s)) for s in (0.0, 20000.0, 24000.0)])
    line, _ = lc.liftoff_from_elements(20000.0, k.elems(20000.0), 24000.0, k.elems(24000.0))
    assert k.calls[0] == pytest.approx(est + 0.45 * tol, rel=1e-12) and line > est > exact
    assert abs(rpm - exact) <= tol / 2 and len(trials) == len(k.calls) <= 8
    assert [t["rpm"] for t in trials] == k.calls
    assert basis.startswith("found by iteration: pressure positive at ")
    assert f"{len(k.calls)} extra solve(s)" in basis and f"tolerance {tol:g} rpm" in basis
    assert "standstill_cold+rated_cold+overspeed_cold" in basis
    lo = max([24000.0] + [t["rpm"] for t in trials if t["closed"]])
    hi = min(t["rpm"] for t in trials if not t["closed"])
    assert hi - lo <= tol and rpm == pytest.approx(0.5 * (lo + hi))


def test_liftoff_search_reproduces_the_linear_extrapolation():
    """On a linear contact the search finds the two-point lift-off speed within its tolerance."""
    k = Contact()
    line, basis_line = lc.liftoff_at_temperature(k.cases())
    rpm, basis, _ = lc.liftoff_search(k.cases(), k.solve, 50.0)
    assert line == pytest.approx(math.sqrt(5.0 / 4.0e-9), rel=1e-9)
    assert abs(rpm - line) <= 50.0 and "found by iteration" in basis


def test_liftoff_search_bracket_from_the_cases():
    """An open overspeed case brackets the search: no trial outside (rated, overspeed)."""
    k = Contact(d=QUARTIC)
    speeds = (("standstill_cold", 0.0), ("rated_cold", 20000.0), ("overspeed_cold", 40000.0))
    rpm, basis, trials = lc.liftoff_search(k.cases(speeds), k.solve, 100.0)
    assert all(20000.0 < s < 40000.0 for s in k.calls)
    # fewer solves than plain bisection of the 20000 rpm bracket to 100 rpm (8)
    assert len(k.calls) < math.ceil(math.log2(20000.0 / 100.0)) and abs(rpm - k.first_opening()) <= 50.0
    assert "standstill_cold+rated_cold+overspeed_cold" in basis


def test_liftoff_search_definition_as_linear_path():
    """Open at standstill -> 0, open at the lowest non-zero speed -> NaN, as on the linear path, with no
    extra solve; no case -> NaN; only a closed standstill -> NaN (no speed scale)."""
    k = Contact(p0=OPEN_AT_STANDSTILL)
    shut = [("standstill_hot", 0.0, k.elems(0.0)), ("rated_hot", 20000.0, k.elems(20000.0))]
    assert lc.liftoff_search(shut, k.solve, 100.0)[:2] == lc.liftoff_at_temperature(shut)
    assert lc.liftoff_at_temperature(shut) == (0.0, "open at standstill_hot")
    k2 = Contact()
    high = [("rated_hot", 40000.0, k2.elems(40000.0)), ("overspeed_hot", 48000.0, k2.elems(48000.0))]
    rpm, basis, _ = lc.liftoff_search(high, k2.solve, 100.0)
    assert math.isnan(rpm) and basis == lc.liftoff_at_temperature(high)[1]
    assert math.isnan(lc.liftoff_search([], k2.solve, 100.0)[0])
    rpm, basis, _ = lc.liftoff_search([("standstill_cold", 0.0, k2.elems(0.0))], k2.solve, 100.0)
    assert math.isnan(rpm) and "non-zero speed" in basis
    assert k.calls == [] and k2.calls == []
    with pytest.raises(ValueError, match="tolerance"):
        lc.liftoff_search(k2.cases(), k2.solve, 0.0)


def test_liftoff_search_bounded_and_failures():
    k = Contact(d=QUARTIC)
    k.f = np.array(SEXTIC)
    rpm, basis, trials = lc.liftoff_search(k.cases(), k.solve, 1e-3, max_solves=2)
    assert math.isnan(rpm) and len(k.calls) == 2 and "tolerance 0.001 rpm not reached" in basis
    # a contact that never opens: trial speeds grow by SEARCH_GROWTH until the cap
    never = Contact(c=[-1e-9] * 4)
    speeds = (("standstill_cold", 0.0), ("rated_cold", 20000.0))
    rpm, basis, trials = lc.liftoff_search(never.cases(speeds), never.solve, 100.0, max_solves=3)
    assert math.isnan(rpm) and "still closed" in basis
    assert never.calls == pytest.approx([30000.0, 45000.0, 67500.0])   # no finite extrapolation: 1.5 x
    # a failing trial run: NaN with the reason, nothing raised

    def boom(s):
        raise RuntimeError("no convergence")

    rpm, basis, trials = lc.liftoff_search(Contact(d=QUARTIC).cases(), boom, 100.0)
    assert math.isnan(rpm) and "failed (RuntimeError: no convergence)" in basis and trials == []


def test_search_applies():
    b = default_params("sleeved_bridgeless", mech_cases=ALL)
    assert lc.search_applies(b) and lc.search_applies(default_params("sleeved_bridgeless"))
    assert not lc.search_applies(b.replace(liftoff_search_tol=0.0))
    assert not lc.search_applies(b.replace(mu_sleeve=0.2))
    assert not lc.search_applies(default_params("sleeved_bridged", mech_cases=ALL))
    with pytest.raises(ValueError, match="liftoff_search_tol"):
        default_params("sleeved_bridgeless", liftoff_search_tol=-1.0)


OPEN_AT_STANDSTILL = [0.0, 1.0, 5.0, 6.0]


class FakeSingle:
    """Stand-in for mapdl._solve_single: records the runs (folder, rpm, degC) and returns the synthetic
    contact of the run's temperature (120 degC = hot)."""

    def __init__(self, cold, hot):
        self.runs, self.cold, self.hot = [], cold, hot

    def __call__(self, q, cfg, wd, model, on_start=None, contact_elements=False, t0=None):
        assert q.mech_cases == "" and contact_elements
        self.runs.append((Path(wd).name, q.mech_rpm, q.mech_temperature))
        k = self.hot if q.mech_temperature == 120.0 else self.cold
        return _fake(0, _contact_elems=k.elems(q.mech_rpm))


def _cases_point(monkeypatch, topo, cold=None, hot=None, force=None, **kw):
    """solve_mech_cases of the six cases (hot = 120 degC) with the fake single-case solver."""
    import sipm.mapdl as mp

    fake = FakeSingle(cold or Contact(d=QUARTIC), hot or Contact(p0=OPEN_AT_STANDSTILL, d=QUARTIC))
    monkeypatch.setattr(mp, "_solve_single", fake)
    p = default_params(topo, **{"mech_temperature": 120.0, "mech_cases": ALL, **kw})
    model = mp.mech_model_from_rotor(build_rotor(p))
    return fake, p, mp.solve_mech_cases(p, Config(min_free_ram_gb=0), Path("unused"), model, liftoff_search=force)


def test_solve_mech_cases_bridgeless_searches(monkeypatch, tmp_path):
    fake, p, r = _cases_point(monkeypatch, "sleeved_bridgeless", mech_rpm=26000.0)
    case_runs = [x for x in fake.runs if not x[0].startswith("liftoff_")]
    extra = [x for x in fake.runs if x[0].startswith("liftoff_")]
    assert [x[0] for x in case_runs] == list(MECH_CASES)
    # cold: first opening near 29.5k rpm, bracketed by rated (closed, 26000) and overspeed (open, 31200);
    # hot: open at standstill, no search
    assert 26000 < fake.cold.first_opening() < 31200
    assert extra and all(x[0].startswith("liftoff_cold_") and x[2] == 20.0 and 26000 < x[1] < 31200 for x in extra)
    assert r["liftoff_rpm_hot"] == 0.0 and r["liftoff_basis_hot"] == "open at standstill_hot"
    assert abs(r["liftoff_rpm_cold"] - fake.cold.first_opening()) <= p.liftoff_search_tol / 2
    assert f"{len(extra)} extra solve(s)" in r["liftoff_basis_cold"]
    assert [t["rpm"] for t in r["_liftoff_trials"]["cold"]] == [x[1] for x in extra]
    assert r["_liftoff_trials"]["hot"] == []
    # the per-case columns are those of the cases alone
    assert r["contact_p_min__overspeed_cold"] == _fake(0)["contact_p_min"]
    # equal temperatures: one search serves both labels
    fake2, _, r2 = _cases_point(monkeypatch, "sleeved_bridgeless", hot=Contact(d=QUARTIC), mech_temperature_cold=120.0)
    assert r2["liftoff_rpm_cold"] == r2["liftoff_rpm_hot"] and all(x[0].startswith("liftoff_cold_")
                                                                    for x in fake2.runs if x[0].startswith("liftoff"))
    # switched off: NaN with the not-linear reason, no extra run
    fake3, _, r3 = _cases_point(monkeypatch, "sleeved_bridgeless", liftoff_search_tol=0.0)
    assert len(fake3.runs) == 6 and math.isnan(r3["liftoff_rpm_cold"]) and "not linear" in r3["liftoff_basis_cold"]


def test_solve_mech_cases_bridged_untouched(monkeypatch):
    """The bridge-retained rotor: the six case runs only and the two-point values; forced on, the search
    runs and agrees with them within its tolerance where the contact is linear."""
    fake, p, r = _cases_point(monkeypatch, "sleeved_bridged", hot=Contact(p0=[x + 2 for x in P0], d=QUARTIC))
    assert [x[0] for x in fake.runs] == list(MECH_CASES) and "_liftoff_trials" not in r
    assert r["liftoff_basis_cold"].startswith("standstill_cold+overspeed_cold")   # two-point line, flagged
    expect = lc.liftoff_outputs(p, {c: {"_contact_elems": (fake.hot if t == 120.0 else fake.cold).elems(s)}
                                    for c in MECH_CASES for s, t in [p.case_speed_temperature(c)]})
    assert {k: r[k] for k in expect} == expect
    fake4, p4, r4 = _cases_point(monkeypatch, "sleeved_bridged", cold=Contact(), hot=Contact(p0=[x + 2 for x in P0]),
                                 force=True)
    assert any(x[0].startswith("liftoff_") for x in fake4.runs)
    assert abs(r4["liftoff_rpm_cold"] - math.sqrt(5.0 / 4.0e-9)) <= p4.liftoff_search_tol
    assert "found by iteration" in r4["liftoff_basis_cold"]


# ---- cache: solver rules without a key change -----------------------------------------------
def test_point_rules():
    from sipm.sweep import RULE_BORE_ARC as B
    from sipm.sweep import RULE_LIFTOFF_SEARCH as S
    from sipm.sweep import point_rules

    bd, bl = default_params("sleeved_bridged"), default_params("sleeved_bridgeless")
    for phys in ("em", "acloss"):
        assert point_rules(bd.replace(sleeve_thickness=1.0), phys) == point_rules(bl, phys) == ()
    assert point_rules(bd, "mech") == () and point_rules(bd.replace(sleeve_thickness=2.0), "mech") == ()
    assert point_rules(bd.replace(mech_cases=ALL), "mech") == ()
    assert point_rules(bd.replace(sleeve_thickness=1.0), "mech") == (B,)
    assert point_rules(bd.replace(mech_mesh=0.8), "mech") == (B,)
    assert point_rules(bl, "mech") == (B,)
    assert point_rules(bl.replace(mech_cases=ALL), "mech") == (B, S)
    assert point_rules(bl.replace(sleeve_thickness=1.5, mech_cases=ALL), "mech") == (S,)
    assert point_rules(bl.replace(sleeve_thickness=1.5, mech_cases=ALL, liftoff_search_tol=0.0), "mech") == ()
    layered = bd.replace(sleeve_thickness=1.0, sleeve_fibre="carbon_sm", sleeve_layers="0:1,30:1")
    assert point_rules(layered, "mech") == ()


def test_sample_table_keys_unchanged():
    """The published sweep's keys (examples/data, written by the first release) are the keys of today,
    for the thin-sleeve points the bore-arc rule changed as for the others."""
    df = pd.read_csv(Path(__file__).parents[1] / "examples/data/sample_mech_results.csv")
    for _, row in df.iterrows():
        p = resolve({"sleeve_thickness": row["sleeve_thickness"], "interference": row["interference"]},
                    topology="sleeved_bridged")
        assert point_key(p, "mech") == row["key"]
    assert point_key(resolve({"sleeve_thickness": 1.0, "interference": 0.0}), "mech") == "25eaaf4afeddc5da"


class CountingSolver:
    def __init__(self):
        self.n = 0

    def __call__(self, p, wd):
        self.n += 1
        return dict(svm_rotor_max=100.0 + self.n)


def test_stale_cache_records_are_resolved(tmp_path):
    """A cached record of a point a rule applies to, written without the rule (earlier code), is solved
    again under the same key; records of other points, and records that list the rule, are reused."""
    cfg = Config(workdir=tmp_path, min_free_ram_gb=0)
    rows = [{"sleeve_thickness": 1.0}, {"sleeve_thickness": 1.5}]
    sw = Sweep("st", "mech", rows, topology="sleeved_bridged", cfg=cfg)
    s1 = CountingSolver()
    sw.run(solver=s1, progress=None)
    assert s1.n == 2
    plan = sw.plan()
    thin, thick = (sw._cache_path(plan[i][3]) for i in (0, 1))
    assert json.loads(thin.read_text())["rules"] == ["bore-arc-mesh"]
    assert "rules" not in json.loads(thick.read_text())          # an unaffected record is as before
    s2 = CountingSolver()
    sw.run(solver=s2, progress=None)
    assert s2.n == 0                                               # both current
    # simulate the caches of the earlier code: no rules in the records
    for f in (thin, thick):
        rec = json.loads(f.read_text())
        rec.pop("rules", None)
        f.write_text(json.dumps(rec))
    df = sw.results()
    assert list(df["status"]) == ["pending", "ok"] and "bore-arc-mesh" in df["message"].iloc[0]
    s3 = CountingSolver()
    df = sw.run(solver=s3, progress=None)
    assert s3.n == 1 and list(df["status"]) == ["ok", "ok"]        # only the thin sleeve was re-solved
    assert [p[3] for p in sw.plan()] == [p[3] for p in plan]       # under the same keys
    assert json.loads(thin.read_text())["rules"] == ["bore-arc-mesh"]


# ---- closed-form lift-off -------------------------------------------------------------------
A, B = 71.1, 72.6
STEEL = dict(E=200e3, nu=0.3, rho=7850e-12)


def _ring(rho=1800e-12):
    return orthotropic(8.8e3, 125e3, 0.28, rho, alpha_r=30e-6, alpha_t=0.0)


@pytest.mark.parametrize("dT", [0.0, 100.0, -40.0])
def test_analytic_liftoff_zero_pressure(dT):
    D, R = isotropic(**STEEL, alpha=12e-6), _ring()
    w = liftoff_speed(A, B, D, R, 0.1, dT)
    assert 0 < w < math.inf
    p0 = coupled_pressure(A, B, D, R, 0.0, 0.1, dT)
    assert abs(coupled_pressure(A, B, D, R, w, 0.1, dT)) < 1e-9 * p0
    assert disc_ring(A, B, D, R, 0.99 * w, 0.1, dT)["in_contact"]
    assert not disc_ring(A, B, D, R, 1.01 * w, 0.1, dT)["in_contact"]
    # p(w) = p0 - c w^2 exactly
    for f in (0.3, 0.7, 1.4):
        assert coupled_pressure(A, B, D, R, f * w, 0.1, dT) == pytest.approx(p0 * (1 - f * f), rel=1e-8, abs=1e-12)


def test_analytic_liftoff_values_and_limits():
    D, R = isotropic(**STEEL, alpha=12e-6), _ring()
    w_cold, w_hot = liftoff_speed(A, B, D, R, 0.1), liftoff_speed(A, B, D, R, 0.1, 100.0)
    assert w_hot > w_cold   # heating adds interference (steel grows, the fibre direction does not)
    assert w_cold * 60 / (2 * math.pi) == pytest.approx(56193.7, abs=0.1)   # docs/validation.md section 11
    assert liftoff_speed(A, B, D, R, -0.05) == 0.0     # clearance: apart at standstill
    assert liftoff_speed(A, B, D, R, 0.0) == 0.0       # line fit, zero pressure at standstill
    # a massless ring does not grow under speed while the disc does: pressure rises, no lift-off
    assert liftoff_speed(A, B, D, _ring(rho=0.0), 0.1) == math.inf
    # coupled_pressure equals disc_ring's pressure in contact and continues it below zero
    assert coupled_pressure(A, B, D, R, 0.5 * w_cold, 0.1) == disc_ring(A, B, D, R, 0.5 * w_cold, 0.1)["p_interface"]
    assert coupled_pressure(A, B, D, R, 1.5 * w_cold, 0.1) < 0


def test_ring_liftoff_reference_uses_params():
    from sipm.validation import ring_liftoff_reference, ring_params

    p = ring_params("orthotropic")
    assert ring_liftoff_reference(p, A, 1.5) == pytest.approx(56193.7, abs=0.1)
    assert ring_liftoff_reference(p, A, 1.5, 120.0) == pytest.approx(76958.9, abs=0.1)
    assert ring_liftoff_reference(p.replace(mech_temperature=120.0), A, 1.5) == pytest.approx(76958.9, abs=0.1)


# ---- sweep, CSV and labels ------------------------------------------------------------------
class CaseSolver:
    """Stand-in for solve_mech with mech_cases: the combined dict of synthetic case results."""

    def __call__(self, p, wd):
        res = {}
        for i, case in enumerate(MECH_CASES):
            rpm, _ = p.case_speed_temperature(case)
            res[case] = _fake(i, _contact_elems=_elems(P0, C, rpm))
        return lc.combine(p, {c: res[c] for c in p.mech_case_list}, wall_s=1.0)


def test_sweep_results_round_trip(tmp_path):
    cfg = Config(workdir=tmp_path, min_free_ram_gb=0)
    sw = Sweep("lc", "mech", [{"interference": 0.1}, {"interference": 0.2}], base={"mech_cases": ALL},
               topology="sleeved_bridged", cfg=cfg)
    df = sw.run(solver=CaseSolver(), progress=None)
    assert (df["status"] == "ok").all()
    rec = json.loads(next((tmp_path / "lc" / "cache").glob("*.json")).read_text())
    assert not any(k.startswith("_") for k in rec["outputs"]) and rec["extra"] == {}   # per-case arrays not cached
    back = pd.read_csv(tmp_path / "lc" / "results.csv")
    for col in ("sh_sleeve_max__overspeed_hot", "sh_sleeve_max__worst", "sh_sleeve_max__worst_case",
                "liftoff_rpm_cold", "liftoff_basis_hot", "mech_case"):
        assert col in back.columns, col
    assert (back["mech_case"] == "rated_hot").all()
    assert back["liftoff_rpm_cold"].iloc[0] == pytest.approx(math.sqrt(5.0 / 4.0e-9), rel=1e-9)


def test_labels_for_envelope_and_cases():
    from sipm.post import LABELS

    for name, how in lc.ENVELOPE:
        lab = LABELS[lc.worst_column(name)]
        assert ("highest" if how == "max" else "lowest") in lab and lab.endswith(")")
    assert LABELS["sh_sleeve_max__overspeed_hot"] == "Max hoop stress, sleeve, overspeed hot (MPa)"
    assert LABELS["liftoff_rpm_hot"] == "Lift-off speed, hot (rpm)"
    assert "sleeve_thickness__rated_hot" not in LABELS
