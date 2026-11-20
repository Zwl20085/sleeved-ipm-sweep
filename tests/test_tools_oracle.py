"""Unit tests of the reference oracle in tools/maxwell_ref (expression evaluator, file-format
parser, exact 2-D geometry). They need no .aedt file; the geometry tests need shapely."""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from maxwell_ref.aedt_expr import Dataset, Evaluator, ExprError, Q, dependency_order, dim_str  # noqa: E402
from maxwell_ref.aedt_parse import parse_args, parse_line, parse_lines, parse_value  # noqa: E402


def ev(dv=None, pv=None, ds=None, ov=None):
    return Evaluator(design_vars=dv or {}, project_vars=pv or {}, datasets=ds or {}, overrides=ov or {})


@pytest.mark.parametrize("expr,value,unit", [
    ("73.92mm", 0.07392, "m"),
    ("1.32 mm", 0.00132, "m"),
    ("2 mm*3", 0.006, "m"),
    ("5e-15mm", 5e-18, "m"),
    ("45deg", math.pi / 4, "1"),
    ("1rad", 1.0, "1"),
    ("4500rpm", 4500 * 2 * math.pi / 60, "s^-1"),
    ("22cel", 295.15, "K"),
    ("1044.8kA_per_meter", 1044800.0, "m^-1*A"),
    ("0.5ms", 5e-4, "s"),
    ("1ohm", 1.0, "m^2*kg*s^-3*A^-2"),
    ("0nH", 0.0, "m^2*kg*s^-2*A^-2"),
])
def test_units(expr, value, unit):
    q = ev().eval(expr)
    assert q.v == pytest.approx(value, rel=1e-15, abs=1e-30)
    assert dim_str(q.dim) == unit


def test_precedence_and_power():
    e = ev()
    assert e.eval("-2^2").v == -4.0
    assert e.eval("2^-1").v == 0.5
    assert e.eval("2^3^2").v == 2 ** 9  # right associative
    assert e.eval("1-2-3").v == -4.0
    assert e.eval("8/2/2").v == 2.0
    assert e.eval("3*-cos(0)").v == -3.0
    assert e.eval("(2*(0.5)^3-3*0.5^2+1)").v == pytest.approx(0.5)


def test_nested_variables_and_dependency_order():
    dv = {"Rsi": "73.92mm", "Wcf": "2mm", "g": "1.32 mm+Wcf", "R": "Rsi - g", "x": "R*cos(45deg)"}
    e = ev(dv)
    assert e.var("g").v == pytest.approx(0.00332)
    assert e.var("x").v == pytest.approx((0.07392 - 0.00332) / math.sqrt(2), rel=1e-15)
    order = dependency_order(dv, set())
    assert order.index("Wcf") < order.index("g") < order.index("R") < order.index("x")


def test_project_variables_and_overrides():
    e = ev({"nn": "$nn_rpm", "T": "60/nn/4/60"}, {"$nn_rpm": "4500"})
    assert e.var("T").v == pytest.approx(60 / 4500 / 4 / 60)
    e2 = ev({"Wcf": "2mm", "g": "1.32mm+Wcf"}, ov={"Wcf": "0.5mm"})
    assert e2.var("g").v == pytest.approx(0.00182)


@pytest.mark.parametrize("expr,value", [
    ("sin(30deg)", 0.5), ("cos(60deg)", 0.5), ("tan(45deg)", 1.0),
    ("atan(1)", math.pi / 4), ("asin(0.5)", math.pi / 6), ("acos(0.5)", math.pi / 3),
    ("atan2(1mm, -1mm)", 3 * math.pi / 4), ("sqrt(16mm*1mm)", 0.004), ("abs(-3mm)", 0.003),
    ("if(2>1, 5, 7)", 5.0), ("if(10deg>30deg,2mm,1mm)", 0.001), ("exp(0)", 1.0), ("ln(1)", 0.0),
    ("min(1mm,2mm)", 0.001), ("max(1,3)", 3.0), ("PI", math.pi), ("pi/2", math.pi / 2),
])
def test_functions(expr, value):
    assert ev().eval(expr).v == pytest.approx(value, rel=1e-14, abs=1e-15)


def test_aedt_dimensionless_promotion_is_recorded():
    e = ev({"Rsi": "73.92mm"})
    q = e.eval("Rsi-1.32/1000")  # AEDT: unitless operand taken as SI metres
    assert q.v == pytest.approx(0.0726, rel=1e-15) and dim_str(q.dim) == "m"
    assert any("dimensionless operand" in n for n in e.notes)


def test_length_field_unitless_is_model_unit():
    e = ev()
    assert e.length_mm("5") == 5.0
    assert e.length_mm("5mm") == pytest.approx(5.0)


def test_pwl_dataset():
    ds = {"$ds1": Dataset("$ds1", (20, 100, 120), (1.0, 0.8998, 0.8695))}
    e = ev(ds=ds)
    assert e.eval("pwl($ds1, 22)").v == pytest.approx(1 - 0.1002 * 2 / 80)
    assert e.eval("pwl($ds1, 110)").v == pytest.approx((0.8998 + 0.8695) / 2)


@pytest.mark.parametrize("expr", ["1mm+1s", "sqrt(-1)", "foo+1", "sin(1mm)", "1 mmm", "2*(3", "1/0", "unknownfn(1)"])
def test_errors_are_loud_and_carry_the_expression(expr):
    with pytest.raises(ExprError) as ei:
        ev().eval(expr)
    assert expr in str(ei.value) or "unknown unit" in str(ei.value)


def test_cycle_detected():
    with pytest.raises(ExprError, match="cyclic"):
        ev({"a": "b+1", "b": "a*2"}).var("a")


def test_compiled_matches_evaluator():
    dv = {"th": "9.142deg", "x": "45.08mm"}
    e = ev(dv)
    expr = "(2*(_t)^3-3*_t^2+1)*(x+1mm) + (_t^3-_t^2)*-cos(th)/1000*2 + if(th>30deg,2mm,1mm)*_t"
    f = e.compile_si(expr, ["_t"], {"_t": Q(0.0)})
    for t in (0.0, 0.13, 0.5, 0.97, 1.0):
        assert f(t) == pytest.approx(e.eval(expr, {"_t": Q(t)}).v, rel=1e-15, abs=1e-18)


def test_parser_lines_and_values():
    f = parse_line("VariableProp('g', 'UD', '', '1.32 mm+Wcf')", 1)
    assert f.kind == "call" and f.value()[3] == "1.32 mm+Wcf"
    f = parse_line("'Conductor number'='Nac'", 2)
    assert f.key == "Conductor number" and f.value() == "Nac"
    assert parse_value("[3: 1, 2, 3]") == [1, 2, 3]
    assert parse_args("933, '22cel', 13108, '22cel'") == [933, "22cel", 13108, "22cel"]
    b = parse_lines(["$begin 'A'", "\tx=1", "\t$begin 'B'", "\t\ty='s'", "\t$end 'B'", "$end 'A'"], 0, 6)
    a = b.items[0]
    assert a["x"] == 1 and a.block("B")["y"] == "s"


# ---- exact geometry ------------------------------------------------------------------------
def _square(g, x0, y0, s, tag):
    pts = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)]
    lp = tuple(g.Line(tag=f"{tag}{i}", p0=pts[i], p1=pts[(i + 1) % 4]) for i in range(4))
    return g.Sheet((g.Face(lp),))


def test_geom2d_boolean_recovers_exact_vertices():
    pytest.importorskip("shapely")
    from maxwell_ref import geom2d as g

    a, b = _square(g, 0.0, 0.0, 2.0, "a"), _square(g, 1.0, 1.0, 2.0, "b")
    log: list = []
    u = g.boolean("union", a, [b], log)
    assert g.sheet_area(u) == pytest.approx(7.0, abs=1e-12)
    verts = sorted({tuple(round(v, 12) for v in s.start()) for s in u.segments()})
    assert (2.0, 1.0) in verts and (1.0, 2.0) in verts and len(verts) == 8
    d = g.boolean("difference", a, [b], log)
    assert g.sheet_area(d) == pytest.approx(3.0, abs=1e-12)


def test_geom2d_arc_area_and_transform():
    pytest.importorskip("shapely")
    from maxwell_ref import geom2d as g

    arc = g.Arc(tag="c", c=(0.0, 0.0), r=2.0, a0=0.0, da=2 * math.pi)
    assert g.loop_area((arc,)) == pytest.approx(4 * math.pi, rel=1e-14)
    m = arc.transform(g.Affine.mirror((0.0, 0.0), (0.0, 1.0)))
    assert m.da == pytest.approx(-2 * math.pi)
    r = arc.transform(g.Affine.rotation(math.radians(30)))
    assert math.dist(r.at(0.0), (2 * math.cos(math.radians(30)), 2 * math.sin(math.radians(30)))) < 1e-14
