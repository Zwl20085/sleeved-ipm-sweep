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
    ("sleeved_bridgeless", "iso"): "ca6115e7fb6ae9ca5df4a103e93c14d1813f33eb262c492a4e339160c0901d33",
    ("sleeved_bridgeless", "hot"): "0624190794dbd31fd1a6b94b49b273222f62cfde2e104f5a958097a333002a28",
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
    assert set(env) == {f"{n}__worst{s}" for n, _ in lc.ENVELOPE for s in ("", "_case")}


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
