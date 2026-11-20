"""Evaluator for AEDT variable expressions.

Values are carried as :class:`Q` (SI magnitude + dimension exponents over
L, M, T, I, Theta).  Angles are dimensionless (radians), as in AEDT.

AEDT rules reproduced here (each one is visible in the projects parsed):

* a number may carry a unit, with or without a space: ``'1.32 mm+Wcf'``;
* a *dimensionless* operand of ``+``/``-``/comparison mixed with a dimensioned
  one is taken as an SI value of that dimension (``'Rsi-1.32/1000'`` is
  Rsi - 1.32 mm).  Every such promotion is recorded in ``Evaluator.notes``;
* trigonometric functions take radians; ``deg`` is the factor pi/180.

Anything else that cannot be evaluated (unknown name, unit, function, a
dimension clash, a cyclic definition, a math domain error) raises
:class:`ExprError` carrying the offending expression.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, Mapping, Optional, Sequence

Dim = tuple  # (L, M, T, I, Theta) as Fractions
DIMLESS: Dim = (Fraction(0),) * 5
_DIM_NAMES = ("m", "kg", "s", "A", "K")


class ExprError(ValueError):
    def __init__(self, msg: str, expr: str = "", chain: Sequence[str] = ()):
        self.expr = expr
        self.chain = list(chain)
        where = f" in expression {expr!r}" if expr else ""
        via = f" (while evaluating {' -> '.join(chain)})" if chain else ""
        super().__init__(f"{msg}{where}{via}")


def _d(L=0, M=0, T=0, I=0, K=0) -> Dim:
    return (Fraction(L), Fraction(M), Fraction(T), Fraction(I), Fraction(K))


def dim_str(d: Dim) -> str:
    parts = []
    for name, e in zip(_DIM_NAMES, d):
        if e == 0:
            continue
        parts.append(name if e == 1 else f"{name}^{e}")
    return "*".join(parts) if parts else "1"


@dataclass(frozen=True)
class Q:
    v: float
    dim: Dim = DIMLESS

    @property
    def dimless(self) -> bool:
        return self.dim == DIMLESS


# --------------------------------------------------------------------------
# units: name -> (factor to SI, dimension, offset)
_PI = math.pi
UNITS: dict[str, tuple[float, Dim, float]] = {}


def _u(names: str, factor: float, dim: Dim, offset: float = 0.0) -> None:
    for n in names.split():
        UNITS[n] = (factor, dim, offset)


_u("m meter", 1.0, _d(L=1))
_u("mm", 1e-3, _d(L=1))
_u("cm", 1e-2, _d(L=1))
_u("um", 1e-6, _d(L=1))
_u("nm", 1e-9, _d(L=1))
_u("km", 1e3, _d(L=1))
_u("in", 0.0254, _d(L=1))
_u("mil", 2.54e-5, _d(L=1))
_u("deg degree", _PI / 180, DIMLESS)
_u("rad", 1.0, DIMLESS)
_u("rpm", 2 * _PI / 60, _d(T=-1))
_u("rad_per_sec", 1.0, _d(T=-1))
_u("deg_per_sec", _PI / 180, _d(T=-1))
_u("Hz", 1.0, _d(T=-1))
_u("kHz", 1e3, _d(T=-1))
_u("MHz", 1e6, _d(T=-1))
_u("s sec", 1.0, _d(T=1))
_u("ms", 1e-3, _d(T=1))
_u("us", 1e-6, _d(T=1))
_u("ns", 1e-9, _d(T=1))
_u("A", 1.0, _d(I=1))
_u("mA", 1e-3, _d(I=1))
_u("kA", 1e3, _d(I=1))
_u("V", 1.0, _d(L=2, M=1, T=-3, I=-1))
_u("mV", 1e-3, _d(L=2, M=1, T=-3, I=-1))
_u("kV", 1e3, _d(L=2, M=1, T=-3, I=-1))
_u("ohm", 1.0, _d(L=2, M=1, T=-3, I=-2))
_u("mOhm mohm", 1e-3, _d(L=2, M=1, T=-3, I=-2))
_u("H", 1.0, _d(L=2, M=1, T=-2, I=-2))
_u("mH", 1e-3, _d(L=2, M=1, T=-2, I=-2))
_u("uH", 1e-6, _d(L=2, M=1, T=-2, I=-2))
_u("nH", 1e-9, _d(L=2, M=1, T=-2, I=-2))
_u("Wb", 1.0, _d(L=2, M=1, T=-2, I=-1))
_u("tesla T", 1.0, _d(M=1, T=-2, I=-1))
_u("N newton", 1.0, _d(L=1, M=1, T=-2))
_u("NewtonMeter", 1.0, _d(L=2, M=1, T=-2))
_u("W", 1.0, _d(L=2, M=1, T=-3))
_u("kW", 1e3, _d(L=2, M=1, T=-3))
_u("kg", 1.0, _d(M=1))
_u("kel", 1.0, _d(K=1))
_u("cel", 1.0, _d(K=1), 273.15)
_u("A_per_meter", 1.0, _d(L=-1, I=1))
_u("kA_per_meter", 1e3, _d(L=-1, I=1))
_u("siemens_per_m siemens_per_meter", 1.0, _d(L=-3, M=-1, T=3, I=2))

CONSTANTS = {"pi": _PI, "PI": _PI, "Pi": _PI}

# --------------------------------------------------------------------------
# tokenizer
_TOK = re.compile(
    r"\s*(?:(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
    r"|(?P<name>\$?[A-Za-z_][A-Za-z_0-9]*)"
    r"|(?P<op><=|>=|==|!=|&&|\|\||[-+*/^(),<>!]))"
)


def tokenize(expr: str) -> list[tuple[str, str]]:
    out, pos = [], 0
    s = expr.strip()
    while pos < len(s):
        m = _TOK.match(s, pos)
        if not m or m.end() == pos:
            raise ExprError(f"unexpected character {s[pos]!r} at {pos}", expr)
        kind = m.lastgroup
        out.append((kind, m.group(kind)))
        pos = m.end()
        while pos < len(s) and s[pos].isspace():
            pos += 1
    return out


# --------------------------------------------------------------------------
# AST
@dataclass(frozen=True)
class Num:
    q: Q


@dataclass(frozen=True)
class Name:
    name: str


@dataclass(frozen=True)
class Unary:
    op: str
    a: object


@dataclass(frozen=True)
class Bin:
    op: str
    a: object
    b: object


@dataclass(frozen=True)
class Call:
    fn: str
    args: tuple


class _Parser:
    def __init__(self, expr: str):
        self.expr = expr
        self.t = tokenize(expr)
        self.i = 0

    def peek(self, k: int = 0) -> Optional[tuple[str, str]]:
        j = self.i + k
        return self.t[j] if j < len(self.t) else None

    def take(self) -> tuple[str, str]:
        tok = self.peek()
        if tok is None:
            raise ExprError("unexpected end of expression", self.expr)
        self.i += 1
        return tok

    def expect(self, op: str) -> None:
        tok = self.take()
        if tok != ("op", op):
            raise ExprError(f"expected {op!r}, got {tok[1]!r}", self.expr)

    def parse(self):
        node = self.or_()
        if self.peek() is not None:
            raise ExprError(f"trailing token {self.peek()[1]!r}", self.expr)
        return node

    def _binary(self, sub, ops):
        node = sub()
        while self.peek() is not None and self.peek()[0] == "op" and self.peek()[1] in ops:
            op = self.take()[1]
            node = Bin(op, node, sub())
        return node

    def or_(self):
        return self._binary(self.and_, ("||",))

    def and_(self):
        return self._binary(self.cmp, ("&&",))

    def cmp(self):
        return self._binary(self.add, ("<", ">", "<=", ">=", "==", "!="))

    def add(self):
        return self._binary(self.mul, ("+", "-"))

    def mul(self):
        return self._binary(self.unary, ("*", "/"))

    def unary(self):
        tok = self.peek()
        if tok in (("op", "-"), ("op", "+"), ("op", "!")):
            self.take()
            return Unary(tok[1], self.unary())
        return self.power()

    def power(self):
        base = self.primary()
        if self.peek() == ("op", "^"):
            self.take()
            return Bin("^", base, self.unary())  # right associative, binds tighter than unary minus on the left
        return base

    def primary(self):
        kind, val = self.take()
        if kind == "num":
            q = Q(float(val))
            nxt = self.peek()
            if nxt is not None and nxt[0] == "name" and nxt[1] in UNITS:
                self.take()
                f, dim, off = UNITS[nxt[1]]
                q = Q(q.v * f + off, dim)
            elif nxt is not None and nxt[0] == "name":
                raise ExprError(f"unknown unit {nxt[1]!r} after number {val}", self.expr)
            return Num(q)
        if kind == "name":
            if self.peek() == ("op", "("):
                self.take()
                args = []
                if self.peek() != ("op", ")"):
                    args.append(self.or_())
                    while self.peek() == ("op", ","):
                        self.take()
                        args.append(self.or_())
                self.expect(")")
                return Call(val, tuple(args))
            return Name(val)
        if (kind, val) == ("op", "("):
            node = self.or_()
            self.expect(")")
            return node
        raise ExprError(f"unexpected token {val!r}", self.expr)


_PARSE_CACHE: dict[str, object] = {}


def parse(expr: str):
    if expr not in _PARSE_CACHE:
        _PARSE_CACHE[expr] = _Parser(expr).parse()
    return _PARSE_CACHE[expr]


# --------------------------------------------------------------------------
# datasets (piecewise linear)
@dataclass(frozen=True)
class Dataset:
    name: str
    x: tuple
    y: tuple

    def pwl(self, xv: float) -> float:
        xs, ys = self.x, self.y
        if xv <= xs[0]:
            i = 0
        elif xv >= xs[-1]:
            i = len(xs) - 2
        else:
            i = max(k for k in range(len(xs) - 1) if xs[k] <= xv)
        x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
        return y0 + (y1 - y0) * (xv - x0) / (x1 - x0)  # AEDT pwl extrapolates linearly


# --------------------------------------------------------------------------
@dataclass
class Evaluator:
    """Evaluates expressions against a design + project variable table."""

    design_vars: Mapping[str, str]
    project_vars: Mapping[str, str] = field(default_factory=dict)
    datasets: Mapping[str, Dataset] = field(default_factory=dict)
    overrides: Mapping[str, str] = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def __post_init__(self) -> None:
        self._cache: dict[str, Q] = {}
        self._stack: list[str] = []

    # -- variables -----------------------------------------------------
    def raw(self, name: str) -> str:
        if name in self.overrides:
            return self.overrides[name]
        if name.startswith("$"):
            if name not in self.project_vars:
                raise ExprError(f"unknown project variable {name!r}", "", self._stack)
            return self.project_vars[name]
        if name not in self.design_vars:
            raise ExprError(f"unknown variable {name!r}", "", self._stack)
        return self.design_vars[name]

    def var(self, name: str) -> Q:
        if name in self._cache:
            return self._cache[name]
        if name in self._stack:
            raise ExprError(f"cyclic definition of {name!r}", "", self._stack + [name])
        self._stack.append(name)
        try:
            q = self.eval(self.raw(name), _ctx=name)
        finally:
            self._stack.pop()
        self._cache[name] = q
        return q

    def is_var(self, name: str) -> bool:
        return name in self.overrides or name in self.design_vars or name in self.project_vars

    # -- expressions ---------------------------------------------------
    def eval(self, expr: str, local: Optional[Mapping[str, Q]] = None, _ctx: str = "") -> Q:
        if not isinstance(expr, str):
            raise ExprError(f"expression must be a string, got {type(expr).__name__}", str(expr), self._stack)
        try:
            node = parse(expr)
        except ExprError as e:
            raise ExprError(str(e).split(" in expression")[0], expr, self._stack) from None
        return self._ev(node, expr, local or {})

    def si(self, expr: str, local: Optional[Mapping[str, Q]] = None) -> float:
        return self.eval(expr, local).v

    def length_mm(self, expr: str, local: Optional[Mapping[str, Q]] = None) -> float:
        """Evaluate a model-length field (model units mm); dimensionless => mm.

        AEDT interprets a unitless number typed in a length field in the model
        unit (here mm), unlike inside an expression with mixed terms.
        """
        q = self.eval(expr, local)
        if q.dim == _d(L=1):
            return q.v * 1e3
        if q.dimless:
            self.notes.append(f"unitless length field {expr!r} taken in model units (mm)")
            return q.v
        raise ExprError(f"expected a length, got dimension {dim_str(q.dim)}", expr, self._stack)

    def angle_rad(self, expr: str) -> float:
        q = self.eval(expr)
        if not q.dimless:
            raise ExprError(f"expected an angle, got dimension {dim_str(q.dim)}", expr, self._stack)
        return q.v

    def compile_si(self, expr: str, params: Sequence[str], probe: Mapping[str, Q]) -> Callable:
        """Fast numeric version of ``expr`` as a function of ``params`` (SI floats).

        The expression is first evaluated once with the full unit-checking
        evaluator at ``probe`` (raises on any error); the compiled function then
        performs the same arithmetic on SI magnitudes (AEDT's dimensionless
        promotion uses the SI magnitude, so numeric results are identical).
        """
        self.eval(expr, probe)
        node = parse(expr)
        consts: dict[str, float] = {}

        def py(n) -> str:
            if isinstance(n, Num):
                return repr(n.q.v)
            if isinstance(n, Name):
                if n.name in params:
                    return "_p_" + n.name
                if self.is_var(n.name):
                    key = f"_c{len(consts)}"
                    consts[key] = self.var(n.name).v
                    return key
                if n.name in CONSTANTS:
                    return repr(CONSTANTS[n.name])
                if n.name in UNITS:
                    return repr(UNITS[n.name][0])
                raise ExprError(f"unknown name {n.name!r}", expr)
            if isinstance(n, Unary):
                return {"-": f"(-{py(n.a)})", "+": py(n.a), "!": f"(0.0 if {py(n.a)} else 1.0)"}[n.op]
            if isinstance(n, Bin):
                a, b = py(n.a), py(n.b)
                if n.op == "^":
                    return f"({a}**{b})"
                if n.op in ("&&", "||"):
                    return f"(1.0 if ({a} {'and' if n.op == '&&' else 'or'} {b}) else 0.0)"
                if n.op in ("<", ">", "<=", ">=", "==", "!="):
                    return f"(1.0 if {a} {n.op} {b} else 0.0)"
                return f"({a}{n.op}{b})"
            if isinstance(n, Call):
                args = [py(a) for a in n.args]
                if n.fn == "if":
                    return f"({args[1]} if {args[0]} else {args[2]})"
                fmap = {"sin": "math.sin", "cos": "math.cos", "tan": "math.tan", "asin": "math.asin",
                        "acos": "math.acos", "atan": "math.atan", "atan2": "math.atan2", "sqrt": "math.sqrt",
                        "abs": "abs", "exp": "math.exp", "ln": "math.log", "log": "math.log", "log10": "math.log10"}
                if n.fn not in fmap:
                    raise ExprError(f"function {n.fn!r} not compilable", expr)
                return f"{fmap[n.fn]}({', '.join(args)})"
            raise ExprError(f"cannot compile {n!r}", expr)

        src = f"lambda {', '.join('_p_' + p for p in params)}: {py(node)}"
        return eval(src, {"math": math, **consts})  # noqa: S307 - source built from our own AST

    def _promote(self, a: Q, b: Q, op: str, expr: str) -> tuple[Q, Q]:
        if a.dim == b.dim:
            return a, b
        if a.dimless:
            self.notes.append(f"dimensionless operand of {op!r} taken as SI {dim_str(b.dim)} in {expr!r}")
            return Q(a.v, b.dim), b
        if b.dimless:
            self.notes.append(f"dimensionless operand of {op!r} taken as SI {dim_str(a.dim)} in {expr!r}")
            return a, Q(b.v, a.dim)
        raise ExprError(f"dimension clash {dim_str(a.dim)} {op} {dim_str(b.dim)}", expr, self._stack)

    def _ev(self, n, expr: str, local: Mapping[str, Q]) -> Q:
        if isinstance(n, Num):
            return n.q
        if isinstance(n, Name):
            if n.name in local:
                return local[n.name]
            if self.is_var(n.name):
                return self.var(n.name)
            if n.name in CONSTANTS:
                return Q(CONSTANTS[n.name])
            if n.name in UNITS and UNITS[n.name][2] == 0.0:
                f, dim, _ = UNITS[n.name]
                return Q(f, dim)
            raise ExprError(f"unknown name {n.name!r}", expr, self._stack)
        if isinstance(n, Unary):
            a = self._ev(n.a, expr, local)
            if n.op == "-":
                return Q(-a.v, a.dim)
            if n.op == "+":
                return a
            return Q(0.0 if a.v else 1.0)
        if isinstance(n, Bin):
            a = self._ev(n.a, expr, local)
            b = self._ev(n.b, expr, local)
            op = n.op
            if op in ("+", "-"):
                a, b = self._promote(a, b, op, expr)
                return Q(a.v + b.v if op == "+" else a.v - b.v, a.dim)
            if op == "*":
                return Q(a.v * b.v, tuple(x + y for x, y in zip(a.dim, b.dim)))
            if op == "/":
                if b.v == 0:
                    raise ExprError("division by zero", expr, self._stack)
                return Q(a.v / b.v, tuple(x - y for x, y in zip(a.dim, b.dim)))
            if op == "^":
                if not b.dimless:
                    raise ExprError("exponent must be dimensionless", expr, self._stack)
                p = Fraction(b.v).limit_denominator(1000)
                if not a.dimless and abs(float(p) - b.v) > 1e-12:
                    raise ExprError("non-rational power of a dimensioned value", expr, self._stack)
                if a.v < 0 and p.denominator != 1:
                    raise ExprError("fractional power of a negative value", expr, self._stack)
                return Q(a.v ** b.v, tuple(x * p for x in a.dim))
            if op in ("<", ">", "<=", ">=", "==", "!="):
                a, b = self._promote(a, b, op, expr)
                r = {"<": a.v < b.v, ">": a.v > b.v, "<=": a.v <= b.v, ">=": a.v >= b.v,
                     "==": a.v == b.v, "!=": a.v != b.v}[op]
                return Q(1.0 if r else 0.0)
            if op == "&&":
                return Q(1.0 if (a.v and b.v) else 0.0)
            if op == "||":
                return Q(1.0 if (a.v or b.v) else 0.0)
        if isinstance(n, Call):
            return self._call(n, expr, local)
        raise ExprError(f"cannot evaluate node {n!r}", expr, self._stack)

    def _call(self, n: Call, expr: str, local: Mapping[str, Q]) -> Q:
        fn = n.fn
        if fn == "if":
            if len(n.args) != 3:
                raise ExprError("if() takes 3 arguments", expr, self._stack)
            c = self._ev(n.args[0], expr, local)
            return self._ev(n.args[1] if c.v else n.args[2], expr, local)
        if fn == "pwl":
            if len(n.args) != 2 or not isinstance(n.args[0], Name):
                raise ExprError("pwl(dataset, x) expected", expr, self._stack)
            ds = self.datasets.get(n.args[0].name)
            if ds is None:
                raise ExprError(f"unknown dataset {n.args[0].name!r}", expr, self._stack)
            x = self._ev(n.args[1], expr, local)
            return Q(ds.pwl(x.v))
        args = [self._ev(a, expr, local) for a in n.args]
        try:
            return _apply(fn, args)
        except ExprError as e:
            raise ExprError(str(e), expr, self._stack) from None
        except (ValueError, OverflowError) as e:
            raise ExprError(f"math error in {fn}(): {e}", expr, self._stack) from None


def _need(fn: str, args: list[Q], n: int) -> None:
    if len(args) != n:
        raise ExprError(f"{fn}() takes {n} argument(s), got {len(args)}")


def _dimless(fn: str, a: Q) -> float:
    if not a.dimless:
        raise ExprError(f"{fn}() needs a dimensionless argument, got {dim_str(a.dim)}")
    return a.v


_TRIG: dict[str, Callable[[float], float]] = {
    "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "asin": math.asin, "acos": math.acos, "atan": math.atan,
    "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    "exp": math.exp, "ln": math.log, "log10": math.log10,
}


def _apply(fn: str, args: list[Q]) -> Q:
    if fn in _TRIG:
        _need(fn, args, 1)
        return Q(_TRIG[fn](_dimless(fn, args[0])))
    if fn == "log":
        _need(fn, args, 1)
        return Q(math.log(_dimless(fn, args[0])))
    if fn == "atan2":
        _need(fn, args, 2)
        y, x = args
        if y.dim != x.dim and not (y.dimless or x.dimless):
            raise ExprError("atan2() arguments of different dimension")
        return Q(math.atan2(y.v, x.v))
    if fn == "sqrt":
        _need(fn, args, 1)
        a = args[0]
        if a.v < 0:
            raise ExprError(f"sqrt() of negative value {a.v!r}")
        return Q(math.sqrt(a.v), tuple(x / 2 for x in a.dim))
    if fn == "abs":
        _need(fn, args, 1)
        return Q(abs(args[0].v), args[0].dim)
    if fn in ("sgn", "sign"):
        _need(fn, args, 1)
        return Q(float((args[0].v > 0) - (args[0].v < 0)))
    if fn in ("min", "max"):
        if len(args) < 1 or any(a.dim != args[0].dim for a in args):
            raise ExprError(f"{fn}() needs arguments of one dimension")
        pick = min if fn == "min" else max
        return Q(pick(a.v for a in args), args[0].dim)
    if fn in ("floor", "ceil", "int"):
        _need(fn, args, 1)
        f = {"floor": math.floor, "ceil": math.ceil, "int": math.trunc}[fn]
        return Q(float(f(args[0].v)), args[0].dim)
    if fn == "mod":
        _need(fn, args, 2)
        return Q(math.fmod(args[0].v, args[1].v), args[0].dim)
    if fn == "pow":
        _need(fn, args, 2)
        if not args[1].dimless or not args[0].dimless:
            raise ExprError("pow() supported for dimensionless arguments only")
        return Q(args[0].v ** args[1].v)
    raise ExprError(f"unknown function {fn!r}")


def names_in(expr: str) -> set[str]:
    """Variable-like names referenced by an expression (functions excluded)."""
    out: set[str] = set()

    def walk(n) -> None:
        if isinstance(n, Name):
            out.add(n.name)
        elif isinstance(n, Unary):
            walk(n.a)
        elif isinstance(n, Bin):
            walk(n.a)
            walk(n.b)
        elif isinstance(n, Call):
            for a in n.args:
                walk(a)

    walk(parse(expr))
    return out


def dependency_order(table: Mapping[str, str], known: set[str]) -> list[str]:
    """Topological order of ``table`` (name -> expression); raises on cycles."""
    order: list[str] = []
    state: dict[str, int] = {}

    def visit(name: str, chain: list[str]) -> None:
        st = state.get(name, 0)
        if st == 2:
            return
        if st == 1:
            raise ExprError(f"cyclic definition of {name!r}", "", chain + [name])
        state[name] = 1
        for dep in sorted(names_in(table[name])):
            if dep in table:
                visit(dep, chain + [name])
        state[name] = 2
        order.append(name)

    for nm in table:
        visit(nm, [])
    return order
