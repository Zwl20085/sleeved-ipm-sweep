"""Replay of the AEDT 2-D modeler history (``GeometryOperations``) with exact
boundaries (see geom2d).

Entity naming: AEDT refers to edges/vertices by integer IDs.  Primitive ops
number them sequentially from ``StartEdgeID``/``StartVertexID`` (edge k =
segment k, vertex k = start point of segment k); clones renumber in source ID
order; later ops publish positions of new entities in their
``OperationIdentity``.  In *nominal* mode those cached positions are attached
to the replayed geometry (and serve as an independent check), and the
resulting ID -> provenance-tag map is stored so that *variant* replays (other
variable values) resolve the same IDs by tag.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Any, Optional

from . import geom2d as g
from .aedt_expr import Evaluator, ExprError, Q
from .aedt_parse import Block
from .geom2d import Affine, Arc, Curve, Face, GeomError, Line, Sheet, Wire, dist

POS_TOL = 1e-5  # mm, cached TolVt tolerance is 5e-7 mm; we allow float noise of the replay


class Unresolved(RuntimeError):
    """Operation that cannot be reproduced from the text file."""


@dataclass
class PartRec:
    pid: int
    name: str
    attrs: dict
    ops: list  # Blocks
    toplevel: bool
    line: int


@dataclass
class Check:
    op_id: int
    part: str
    kind: str
    ent_id: Optional[int]
    cached: tuple
    mine: Optional[tuple]
    dev: float


def _op_params(op: Block) -> Block:
    for b in op.blocks():
        if b.name.endswith("Parameters"):
            return b
    raise KeyError(f"operation {op.get('ID')} has no parameter block")


def _identity_helper(op: Block) -> Optional[Block]:
    ident = op.opt_block("OperationIdentity")
    if ident is None:
        return None
    return ident.opt_block("GeomTopolBasedOperationIdentityHelper")


class GeometryModel:
    def __init__(self, design: Block, ev: Evaluator, naming: Optional[dict] = None,
                 donors: Optional[dict] = None):
        """``donors``: {part_id: (design_block, label)} - a design (e.g. in another
        project) holding the full construction history of a part whose history
        was purged here; it is replayed with THIS design's variables."""
        self.ev = ev
        self.donors = dict(donors or {})
        self._donor_models: dict = {}
        self.donor_used: dict = {}
        self.nominal = naming is None
        self.naming: dict = {} if naming is None else naming  # (op_id) -> {aedt_id: ('e'|'v', signature)}
        self.mirror_rule: dict = dict(self.naming.get("mirror_rules", {}))
        go = design.block("ModelSetup").block("GeometryCore").block("GeometryOperations")
        self.parts: dict[int, PartRec] = {}
        self.op_owner: dict[int, tuple[int, int]] = {}
        for grp, top in (("ToplevelParts", True), ("OperandParts", False)):
            for gp in go.block(grp).blocks("GeometryPart"):
                ops = [b for b in gp.block("Operations").blocks()]
                if not ops:
                    continue
                pid = ops[0]["ParentPartID"]
                attrs = gp.block("Attributes").to_dict()
                self.parts[pid] = PartRec(pid, attrs["Name"], attrs, ops, top, gp.line)
                for i, op in enumerate(ops):
                    self.op_owner[op["ID"]] = (pid, i)
        self.cs_ops: dict[int, Block] = {op["ID"]: op for op in go.block("CoordinateSystems").blocks()}
        self._cache: dict[tuple, Any] = {}
        self._cs_cache: dict[int, Affine] = {}
        self.checks: list[Check] = []
        self.offplane: dict = {}
        self.log: dict[str, list] = {}
        self.notes: list[str] = []

    # ------------------------------------------------------------------
    def part_by_name(self, name: str) -> PartRec:
        for p in self.parts.values():
            if p.name == name:
                return p
        raise KeyError(name)

    def final_body(self, pid: int):
        return self.body(pid, len(self.parts[pid].ops) - 1)

    def body(self, pid: int, upto: int):
        key = (pid, upto)
        if key in self._cache:
            val = self._cache[key]
            if isinstance(val, Exception):
                raise val
            return val
        part = self.parts[pid]
        prev = self.body(pid, upto - 1) if upto > 0 else None
        op = part.ops[upto]
        try:
            val = self._apply(part, op, prev)
            if self.nominal:
                val = self._attach_cached(part, op, val)
            else:
                val = self._attach_named(op, val)
        except (Unresolved, GeomError, ExprError, KeyError, ValueError, ZeroDivisionError) as e:
            err = e if isinstance(e, Unresolved) else Unresolved(
                f"{part.name}: op {op.get('OperationType')} id={op.get('ID')}: {type(e).__name__}: {e}")
            self._cache[key] = err
            raise err from None
        self._cache[key] = val
        return val

    def state_at_op(self, op_id: int):
        """Body state right after op ``op_id`` (clone placeholders: before)."""
        pid, i = self.op_owner[op_id]
        return self.body(pid, i)

    # ------------------------------------------------------------------
    # coordinate systems
    def cs(self, cs_id: int) -> Affine:
        if cs_id in (1, None):
            return Affine()
        if cs_id in self._cs_cache:
            return self._cs_cache[cs_id]
        op = self.cs_ops.get(cs_id)
        if op is None:
            raise Unresolved(f"coordinate system {cs_id} not found")
        typ = op["OperationType"]
        if typ == "CreateRelativeCoordinateSystem":
            prm = op.block("RelativeCSParameters")
            if prm["Mode"] != "Axis/Position":
                raise Unresolved(f"relative CS mode {prm['Mode']!r} not supported")
            L = self.ev.length_mm
            o = (L(prm["OriginX"]), L(prm["OriginY"]))
            if abs(L(prm["OriginZ"])) > 0 or abs(L(prm["XAxisZvec"])) > 0 or abs(L(prm["YAxisZvec"])) > 0:
                raise Unresolved(f"CS {cs_id} not in the XY plane")
            x = (L(prm["XAxisXvec"]), L(prm["XAxisYvec"]))
            y = (L(prm["YAxisXvec"]), L(prm["YAxisYvec"]))
            local = _frame(o, x, y, x_driven=True)
            parent = self.cs(op["ReferenceCoordSystemID"])
            T = local.then(parent)
        elif typ == "CreateObjectCoordinateSystem":
            T = self.object_cs(op)[0]
        else:
            raise Unresolved(f"CS type {typ}")
        self._cs_cache[cs_id] = T
        return T

    def object_cs(self, cs_op: Block) -> tuple[Affine, dict]:
        ph_id = cs_op["PlaceHolderOperationID"]
        if ph_id not in self.op_owner:
            raise Unresolved(f"object CS {cs_op.block('Attributes')['Name']}: placeholder {ph_id} not in history")
        pid, idx = self.op_owner[ph_id]
        ph = self.parts[pid].ops[idx]
        prm = ph.block("ObjectCSParameters")
        body = self.final_body(pid) if prm.get("MoveToEnd") else self.body(pid, idx)
        info: dict = {"part": self.parts[pid].name, "placeholder_op": ph_id,
                      "DrivenByXAxis": prm.get("DrivenByXAxis"), "MoveToEnd": prm.get("MoveToEnd")}
        o, o_cached = self._cs_point(body, prm.block("Origin"))
        info["origin_mm"], info["origin_cached_mm"] = o, o_cached
        xb = prm.opt_block("xAxisPos")
        if xb is None:
            raise Unresolved("object CS without xAxisPos")
        xp, xp_cached = self._cs_point(body, xb)
        info["xaxis_point_mm"], info["xaxis_point_cached_mm"] = xp, xp_cached
        yb = prm.opt_block("yAxis")
        if yb is None or yb["DirectionType"] != "AbsoluteDirection":
            raise Unresolved("object CS y axis not an absolute direction")
        yd = (float(yb["xDirection"]), float(yb["yDirection"]))
        info["yaxis_direction_stored"] = yd
        x_vec = (xp[0] - o[0], xp[1] - o[1])
        T_x = _frame(o, x_vec, yd, x_driven=True)
        T_y = _frame(o, x_vec, yd, x_driven=False)
        info["x_axis_deg_if_x_driven"] = math.degrees(math.atan2(T_x.c, T_x.a))
        info["x_axis_deg_if_y_driven"] = math.degrees(math.atan2(T_y.c, T_y.a))
        T = T_x if prm.get("DrivenByXAxis") else T_y
        if prm.get("ReverseXAxis") or prm.get("ReverseYAxis"):
            raise Unresolved("object CS axis reversal flags not supported")
        diff = abs((info["x_axis_deg_if_x_driven"] - info["x_axis_deg_if_y_driven"] + 180) % 360 - 180)
        info["x_axis_deg_reading_of_flag"] = math.degrees(math.atan2(T.c, T.a))
        info["axis_status"] = ("unambiguous (both readings agree)" if diff < 1e-6 else
                               f"AMBIGUOUS: the stored absolute y direction is not perpendicular to origin->xAxisPos; "
                               f"x-exact gives {info['x_axis_deg_if_x_driven']:.4f} deg, y-exact gives "
                               f"{info['x_axis_deg_if_y_driven']:.4f} deg (difference {diff:.4f} deg). "
                               f"DrivenByXAxis={prm.get('DrivenByXAxis')}; "
                               "its AEDT semantics are not stated in the file.")
        info["handedness_det"] = T.det
        return T, info

    def _cs_point(self, body, b: Block) -> tuple[tuple, tuple]:
        cached = (float(b["XPosition"]), float(b["YPosition"]))
        if not b.get("IsAttachedToEntity"):
            return cached, cached
        ent = b["EntityID"]
        pt = b["PositionType"]
        if pt == "OnVertex":
            p = self.vertex_pos(body, ent)
        elif pt == "EdgeCenter":
            li, i = self._one(g.find_edge(body, ent), "edge", ent)
            p = body.loops()[li][i].at(float(b.get("UParam", 0.5)))
        else:
            raise Unresolved(f"CS position type {pt}")
        return p, cached

    def vertex_pos(self, body, vid: int) -> tuple:
        hits = g.find_vertex(body, vid)
        li, i = self._one(hits, "vertex", vid)
        lp = body.loops()[li]
        return lp[i].start() if i < len(lp) else lp[-1].end()

    @staticmethod
    def _one(hits: list, kind: str, ent: int):
        if len(hits) != 1:
            raise Unresolved(f"{kind} id {ent} found {len(hits)} times")
        return hits[0]

    # ------------------------------------------------------------------
    def _apply(self, part: PartRec, op: Block, prev):
        typ = op["OperationType"]
        fn = getattr(self, "_op_" + typ, None)
        if fn is None:
            raise Unresolved(f"operation type {typ} not supported")
        return fn(part, op, prev)

    # primitives -------------------------------------------------------
    def _ids(self, op: Block) -> tuple[int, int]:
        ident = op.opt_block("OperationIdentity")
        if ident is None:
            return -1, -1
        return ident.get("StartEdgeID", -1), ident.get("StartVertexID", -1)

    def _number(self, segs: list, op: Block, closed: bool) -> list:
        e0, v0 = self._ids(op)
        out = []
        n = len(segs)
        for k, s in enumerate(segs):
            vs = frozenset({v0 + k}) if v0 > 0 else frozenset()
            ve = frozenset({v0 + (k + 1) % n if closed else v0 + k + 1}) if v0 > 0 else frozenset()
            out.append(s.with_meta(tag=f"e{op['ID']}.{k}", ids=frozenset({e0 + k}) if e0 > 0 else frozenset(),
                                   v0=vs, v1=ve))
        return out

    def _op_Polyline(self, part, op, prev):
        prm = op.block("PolylineParameters")
        T = self.cs(op["ReferenceCoordSystemID"])
        L = self.ev.length_mm
        pts = []
        zs = set()
        for p in prm.block("PolylinePoints").blocks("PLPoint"):
            zs.add(round(L(p["Z"]), 12))
            pts.append((L(p["X"]), L(p["Y"])))
        if len(zs) != 1:
            raise Unresolved(f"polyline not in a plane parallel to XY (z values {sorted(zs)})")
        z = zs.pop()
        if z != 0.0:
            # only acceptable if a later sweep along z brings it back through z=0
            self.offplane[part.pid] = z
        segs = []
        idx = 0
        for sb in prm.block("PolylineSegments").blocks("PLSegment"):
            st = sb.get("StartIndex", idx)
            npt = sb["NoOfPoints"]
            typ = sb["SegmentType"]
            if typ == "Line":
                if npt != 2:
                    raise Unresolved("line segment with NoOfPoints != 2")
                segs.append(Line(tag="", p0=pts[st], p1=pts[st + 1]))
            elif typ == "AngularArc":
                if sb["ArcPlane"] != "XY" or str(sb.get("NoOfSegments", "0")) != "0":
                    raise Unresolved(f"angular arc plane {sb['ArcPlane']} / faceted {sb.get('NoOfSegments')}")
                c = (L(sb["ArcCenterX"]), L(sb["ArcCenterY"]))
                ang = self.ev.angle_rad(sb["ArcAngle"])
                p0 = pts[st]
                r = dist(p0, c)
                a0 = math.atan2(p0[1] - c[1], p0[0] - c[0])
                arc = Arc(tag="", c=c, r=r, a0=a0, da=ang)
                stored_end = pts[st + npt - 1]
                self._record_stored(part, op, "angular-arc end (stored at creation)", stored_end, arc.end())
                pts[st + 1] = arc.at(0.5)
                pts[st + npt - 1] = arc.end()
                segs.append(arc)
            elif typ == "Arc":
                if npt != 3 or str(sb.get("NoOfSegments", "0")) != "0":
                    raise Unresolved("3-point arc with unexpected point count / faceting")
                segs.append(_arc3(pts[st], pts[st + 1], pts[st + 2]))
            else:
                raise Unresolved(f"polyline segment type {typ}")
            idx = st + npt - 1
        closed = bool(prm["IsPolylineClosed"])
        gap = dist(segs[-1].end(), segs[0].start())
        if closed and gap > 1e-9:
            raise Unresolved(f"closed polyline does not return to its first point (gap {gap:.3g} mm)")
        segs = [s.transform(T) for s in segs]
        # a polyline whose last point equals its first has as many vertices as segments
        # (AEDT NumNewVertices confirms this for the 'open' magnet outlines)
        topo_closed = closed or gap <= 1e-9
        ident = op.opt_block("OperationIdentity")
        if ident is not None and ident.has("NumNewVertices") and ident["NumNewVertices"] > 0:
            nv = ident["NumNewVertices"]
            if nv != (len(segs) if topo_closed else len(segs) + 1):
                raise Unresolved(f"polyline vertex count {nv} inconsistent with {len(segs)} segments")
        segs = self._number(segs, op, topo_closed)
        return Wire((tuple(segs),))

    def _record_stored(self, part, op, kind, stored, mine):
        self.checks.append(Check(op["ID"], part.name, kind, None, tuple(stored), tuple(mine), dist(stored, mine)))

    def _op_Rectangle(self, part, op, prev):
        prm = op.block("RectangleParameters")
        if prm["WhichAxis"] != "Z":
            raise Unresolved("rectangle not in XY")
        L = self.ev.length_mm
        x0, y0 = L(prm["XStart"]), L(prm["YStart"])
        w, h = L(prm["Width"]), L(prm["Height"])
        P = [(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)]
        segs = [Line(tag="", p0=P[k], p1=P[(k + 1) % 4]) for k in range(4)]
        T = self.cs(op["ReferenceCoordSystemID"])
        segs = self._number([s.transform(T) for s in segs], op, True)
        return Wire((tuple(segs),))

    def _op_Circle(self, part, op, prev):
        prm = op.block("CircleParameters")
        if prm["WhichAxis"] != "Z" or str(prm.get("NumSegments", "0")) != "0":
            raise Unresolved("circle not in XY or faceted")
        L = self.ev.length_mm
        arc = Arc(tag="", c=(L(prm["XCenter"]), L(prm["YCenter"])), r=L(prm["Radius"]), a0=0.0, da=g.TWO_PI)
        arc = arc.transform(self.cs(op["ReferenceCoordSystemID"]))
        e0, v0 = self._ids(op)
        arc = arc.with_meta(tag=f"e{op['ID']}.0", ids=frozenset({e0}) if e0 > 0 else frozenset())
        return Wire(((arc,),))

    def _op_EquationBasedCurve(self, part, op, prev):
        prm = op.block("EquationBasedCurveParameters")
        if prm["ZtFunction"].strip() not in ("0", "0mm"):
            raise Unresolved("equation curve with non-zero z")
        if int(prm.get("NumOfPointsOnCurve", 0) or 0) != 0:
            raise Unresolved("equation curve with NumOfPointsOnCurve != 0 (faceted)")
        xs, ys = prm["XtFunction"], prm["YtFunction"]
        t0, t1 = self.ev.si(prm["tStart"]), self.ev.si(prm["tEnd"])
        ev = self.ev

        # unit-checked evaluation at both ends, then a compiled numeric form
        for tt in (t0, t1, 0.5 * (t0 + t1)):
            ev.length_mm(xs, {"_t": Q(tt)}), ev.length_mm(ys, {"_t": Q(tt)})
        fx = ev.compile_si(xs, ["_t"], {"_t": Q(t0)})
        fy = ev.compile_si(ys, ["_t"], {"_t": Q(t0)})
        for tt in (t0, t1, 0.37 * t0 + 0.63 * t1):
            ref = (ev.length_mm(xs, {"_t": Q(tt)}), ev.length_mm(ys, {"_t": Q(tt)}))
            if abs(fx(tt) * 1e3 - ref[0]) > 1e-12 or abs(fy(tt) * 1e3 - ref[1]) > 1e-12:
                raise Unresolved("compiled equation curve disagrees with the evaluator")

        def f(t: float, _fx=fx, _fy=fy):
            return (_fx(t) * 1e3, _fy(t) * 1e3)
        T = self.cs(op["ReferenceCoordSystemID"])
        cur = Curve(tag="", f=f, t0=t0, t1=t1, T=T,
                    desc=(("x_of_t", xs), ("y_of_t", ys), ("parameter", "_t")))
        e0, v0 = self._ids(op)
        cur = cur.with_meta(tag=f"e{op['ID']}.0", ids=frozenset({e0}) if e0 > 0 else frozenset(),
                            v0=frozenset({v0}) if v0 > 0 else frozenset(),
                            v1=frozenset({v0 + 1}) if v0 > 0 else frozenset())
        return Wire(((cur,),))

    def _op_NativeBody(self, part, op, prev):
        prm = op.block("NativeBodyParameters")
        if part.pid in self.donors:
            dblock, label = self.donors[part.pid]
            if label not in self._donor_models:
                self._donor_models[label] = GeometryModel(dblock, self.ev, naming=None if self.nominal else {})
            dm = self._donor_models[label]
            if part.pid not in dm.parts:
                raise Unresolved(f"donor {label} has no part {part.pid}")
            body = dm.final_body(part.pid)
            self.donor_used[part.name] = {"donor": label, "donor_part": dm.parts[part.pid].name,
                                          "donor_checks": [c.__dict__ for c in dm.checks
                                                           if dm.op_owner.get(c.op_id, (None,))[0] == part.pid]}
            self.notes.append(f"{part.name}: history purged here (Parasolid only); geometry replayed from donor "
                              f"history {label} part {part.pid} with this design's variables")
            return body
        raise Unresolved(f"native (history-purged) body, source {prm.get('SourceFile')!r}, "
                         f"files {op.block('BodyBlock').get('BodyFileNamesVec')}: "
                         "geometry only in binary Parasolid data")

    def _op_CloneFrom(self, part, op, prev):
        prm = op.block("CloneFromParameters")
        src = prm["SourceID"]
        # source state: right before its CloneTo op
        sp = self.parts.get(src)
        if sp is None:
            raise Unresolved(f"clone source part {src} not in this design")
        idx = next(i for i, o in enumerate(sp.ops) if o["OperationType"] == "CloneTo")
        body = self.body(src, idx - 1)
        return body

    # sheet / wire ops ---------------------------------------------------
    def _op_CoverLines(self, part, op, prev):
        if not isinstance(prev, Wire):
            raise Unresolved("CoverLines on a non-wire body")
        return g.cover(prev)

    def _tool(self, pid: int, op_id: int):
        tp, ti = self.op_owner[op_id]
        if tp != pid:
            raise Unresolved(f"tool operation {op_id} does not belong to part {pid}")
        body = self.body(tp, ti)
        if tp in self.offplane:
            raise Unresolved(f"tool part {self.parts[tp].name} lies off the model plane (z={self.offplane[tp]} mm)")
        return body

    def _tools(self, prm: Block, op: Block) -> list:
        parts = prm["ToolParts"] if prm.has("ToolParts") else []
        ops_ = op["ToolOperations"] if op.has("ToolOperations") else []
        if len(parts) != len(ops_):
            raise Unresolved("tool parts/operations mismatch")
        return [self._tool(p, o) for p, o in zip(parts, ops_)]

    def _boolean(self, part, op, prev, kind):
        prm = op.blocks()[0]
        if prm.get("KeepOriginals"):
            # the tool parts are kept as separate objects; the blank result is unchanged
            self.notes.append(f"{part.name}: op {op['ID']} {kind} with KeepOriginals=true (tools kept)")
        tools = self._tools(prm, op)
        if isinstance(prev, Wire):
            if kind != "union" or not all(isinstance(t, Wire) for t in tools):
                raise Unresolved(f"{kind} on wire bodies")
            chains = list(prev.chains) + [c for t in tools for c in t.chains]
            return Wire(tuple(chains))
        if not all(isinstance(t, Sheet) for t in tools):
            raise Unresolved(f"{kind} of a sheet with a wire tool")
        log = self.log.setdefault(part.name, [])
        out = g.boolean(kind, prev, tools, log)
        return _unique_tags(out)

    def _op_Unite(self, part, op, prev):
        return self._boolean(part, op, prev, "union")

    def _op_Substract(self, part, op, prev):
        return self._boolean(part, op, prev, "difference")

    def _op_Intersect(self, part, op, prev):
        return self._boolean(part, op, prev, "intersection")

    def _op_Split(self, part, op, prev):
        prm = op.block("SplitToParameters")
        if prm["ToolType"] != "PlaneTool":
            raise Unresolved(f"split tool {prm['ToolType']}")
        T = self.cs(op["ReferenceCoordSystemID"])
        plane, side = prm["SplitPlane"], prm["WhichSide"]
        # keep the positive side: x>=0 for YZ, y>=0 for ZX (in the reference CS)
        if plane == "YZ":
            p, d = (0.0, 0.0), (0.0, -1.0)  # left of downward direction is +x
        elif plane == "ZX":
            p, d = (0.0, 0.0), (1.0, 0.0)   # left of +x direction is +y
        else:
            raise Unresolved(f"split plane {plane} in a 2-D XY design")
        if side == "NegativeOnly":
            d = (-d[0], -d[1])
        elif side != "PositiveOnly":
            raise Unresolved(f"split side {side}")
        hp = g.half_plane(T(p), T.vec(d), tag=f"split{op['ID']}")
        if isinstance(prev, Wire):
            raise Unresolved("split of a wire body")
        log = self.log.setdefault(part.name, [])
        out = g.boolean("intersection", prev, [hp], log)
        if any(s.tag == "farbox" for s in out.segments()):
            raise GeomError("split half-plane far box reached the body")
        return _unique_tags(out)

    def _op_SplitEdit(self, part, op, prev):
        return prev

    # transforms ------------------------------------------------------------
    def _rot(self, op, angle_expr: str, axis: str, times: float = 1.0) -> Affine:
        if axis != "Z":
            raise Unresolved(f"rotation about {axis} in a 2-D XY design")
        F = self.cs(op["ReferenceCoordSystemID"])
        R = Affine.rotation(times * self.ev.angle_rad(angle_expr))
        return F.inverse().then(R).then(F)

    def _op_Rotate(self, part, op, prev):
        prm = op.block("RotateParameters")
        return prev.transform(self._rot(op, prm["RotateAngle"], prm["RotateAxis"]))

    def _op_Move(self, part, op, prev):
        prm = op.block("TranslateParameters") if op.opt_block("TranslateParameters") else op.block("MoveParameters")
        L = self.ev.length_mm
        if abs(L(prm["TranslateVectorZ"])) > 0:
            raise Unresolved("move out of plane")
        F = self.cs(op["ReferenceCoordSystemID"])
        v = F.vec((L(prm["TranslateVectorX"]), L(prm["TranslateVectorY"])))
        return prev.transform(Affine.translation(*v))

    def _op_ObjectCSPlaceHolder(self, part, op, prev):
        return prev

    def _clone_source(self, op) -> tuple[Any, Block]:
        prm = op.blocks()[0]
        ph = prm["PlaceHolderOpnId"] if prm.has("PlaceHolderOpnId") else op["PlaceHolderOpnId"]
        sp, si = self.op_owner[ph]
        src_op = self.parts[sp].ops[si]
        if src_op.blocks()[0].get("CreateNewObjects") is False:
            raise Unresolved("duplicate with CreateNewObjects=false")
        return self.body(sp, si - 1), src_op

    def _renumber(self, body, op, mirrored: bool = False):
        """New IDs for a clone: source IDs in ascending order map to StartID+rank.

        For a mirror clone AEDT numbers the copy in reversed traversal order
        starting from the same first entity (rank r -> (-r) mod n); this was
        established from the cached positions of later operations on mirror
        clones (MoveEdges/Fillet checks) and is re-verified on every run.
        """
        e0, v0 = self._ids(op)
        loops = sorted(body.loops(), key=lambda lp: min((min(s.ids) for s in lp if s.ids), default=1 << 62))
        emap: dict = {}
        vmap: dict = {}
        eoff = voff = 0
        for lp in loops:
            closed = dist(lp[-1].end(), lp[0].start()) < 1e-9
            eids = sorted({min(s.ids) for s in lp if s.ids})
            vids = sorted({min(vs) for s in lp for vs in (s.v0, s.v1) if vs})

            def rank(k: int, n: int, closed: bool = closed) -> int:
                return (-k) % n if (mirrored and closed) else k
            if e0 > 0:
                emap.update({old: e0 + eoff + rank(k, len(eids)) for k, old in enumerate(eids)})
            if v0 > 0:
                vmap.update({old: v0 + voff + rank(k, len(vids)) for k, old in enumerate(vids)})
            eoff += len(eids)
            voff += len(vids)
        if mirrored and len(loops) > 1:
            self.notes.append(f"op {op['ID']}: mirror-clone renumbering of a multi-loop body is unverified")

        def fix(s):
            return s.with_meta(tag=f"c{op['ID']}:{s.tag}",
                               ids=frozenset({emap[min(s.ids)]}) if s.ids and emap else frozenset(),
                               v0=frozenset({vmap[min(s.v0)]}) if s.v0 and vmap else frozenset(),
                               v1=frozenset({vmap[min(s.v1)]}) if s.v1 and vmap else frozenset())
        return _map_segments(body, fix)

    def _op_DuplicateBodyAroundAxis(self, part, op, prev):
        body, src_op = self._clone_source(op)
        sprm = src_op.block("DuplicateAroundAxisParameters")
        k = op.blocks()[0]["WhichClone"] + 1
        T = self._rot(src_op, sprm["AngleStr"], sprm["WhichAxis"], times=k)
        return self._renumber(body.transform(T), op)

    def _op_DuplicateBodyMirror(self, part, op, prev):
        body, src_op = self._clone_source(op)
        sprm = src_op.block("DuplicateToMirrorParameters")
        L = self.ev.length_mm
        F = self.cs(src_op["ReferenceCoordSystemID"])
        if abs(L(sprm["DuplicateMirrorNormalZ"])) > 0:
            raise Unresolved("mirror normal out of plane")
        base = F((L(sprm["DuplicateMirrorBaseX"]), L(sprm["DuplicateMirrorBaseY"])))
        nrm = F.vec((L(sprm["DuplicateMirrorNormalX"]), L(sprm["DuplicateMirrorNormalY"])))
        rule = self.mirror_rule.get(op["ID"], "reversed")
        return self._renumber(body.transform(Affine.mirror(base, nrm)), op, mirrored=(rule == "reversed"))

    # ------------------------------------------------------------------
    def referenced_ids(self) -> set:
        """Entity IDs referenced by local ops / object CSs anywhere in the design."""
        out: set = set()
        for p in self.parts.values():
            for op in p.ops:
                for b in op.walk():
                    for key in ("EdgesToMove", "Vertices", "Edges"):
                        if b.has(key) and b.fields(key)[0].kind == "array":
                            out.update(v for v in b[key] if isinstance(v, int))
                    if b.name in ("Origin", "xAxisPos", "yAxisPos") and b.has("EntityID"):
                        out.add(b["EntityID"])
        return out

    def resolve_mirror_rules(self) -> list:
        """Choose, per mirror clone whose entities are referenced later, the
        numbering rule ('reversed' / 'identity') that agrees with AEDT's cached
        positions.  Step 1 uses evidence inside the clone part alone; ops left
        undecided are then searched jointly (all combinations) using the
        cached checks of every part that depends on them.  Nominal mode only."""
        refs = self.referenced_ids()
        users: dict = {}
        for p in self.parts.values():
            for dep in self._deps(p):
                users.setdefault(dep, set()).add(p.pid)
        relevant = []
        for p in sorted(self.parts.values(), key=lambda r: r.pid):
            for op in p.ops:
                if op["OperationType"] != "DuplicateBodyMirror":
                    continue
                ident = op.opt_block("OperationIdentity")
                e0, v0 = self._ids(op)
                ne = ident.get("NumNewEdges", 0) if ident else 0
                nv = ident.get("NumNewVertices", 0) if ident else 0
                mine = set(range(e0, e0 + ne)) | set(range(v0, v0 + nv))
                if mine & refs:
                    relevant.append((p, op))
        decisions: dict = {}
        shared: dict = {}

        def trial(rules: dict, parts: set) -> tuple:
            sc = GeometryModel.__new__(GeometryModel)
            sc.__dict__.update(self.__dict__)
            sc._cache = {k: v for k, v in shared.items() if k[0] not in parts}
            sc._cs_cache, sc.checks, sc.log, sc.notes, sc.naming, sc.offplane = {}, [], {}, [], {}, {}
            sc.mirror_rule = dict(self.mirror_rule)
            sc.mirror_rule.update(rules)
            unres = 0
            for a in sorted(parts):
                try:
                    sc.final_body(a)
                except Unresolved:
                    unres += 1
            rel = [c for c in sc.checks if "stored" not in c.kind and self.op_owner.get(c.op_id, (None,))[0] in parts]
            return unres, sum(c.dev > 1e-5 for c in rel), len(rel)

        undecided = []
        for p, op in relevant:  # step 1: local evidence
            sc = {}
            for rule in ("reversed", "identity"):
                sc[rule] = trial({op["ID"]: rule}, {p.pid})
            if sc["reversed"][2] > 0 and sc["reversed"][:2] != sc["identity"][:2]:
                best = min(sc, key=lambda r: sc[r][:2])
                self.mirror_rule[op["ID"]] = best
                decisions[op["ID"]] = {"op": op["ID"], "part": p.name, "rule": best, "evidence": "local",
                                       "scores_unresolved_bad_checks": sc}
            else:
                undecided.append((p, op))
        if undecided:  # step 2: joint search over the remaining ops
            parts = set()
            for p, _ in undecided:
                parts |= self._closure(p.pid, users)
            if len(undecided) > 6:
                raise Unresolved(f"too many coupled mirror clones ({len(undecided)}) for a joint search")
            results = {}
            for combo in itertools.product(("reversed", "identity"), repeat=len(undecided)):
                rules = {op["ID"]: r for (p, op), r in zip(undecided, combo)}
                results[combo] = trial(rules, parts)
            best_combo = min(results, key=lambda c: results[c][:2])
            best_score = results[best_combo][:2]
            winners = [c for c, v in results.items() if v[:2] == best_score]
            for k, ((p, op), r) in enumerate(zip(undecided, best_combo)):
                determined = len({c[k] for c in winners}) == 1
                self.mirror_rule[op["ID"]] = r
                decisions[op["ID"]] = {"op": op["ID"], "part": p.name, "rule": r,
                                       "evidence": "joint search over cached checks of dependent parts" if determined
                                       else "none: both rules agree equally with the cached data (choice unverified)",
                                       "joint_best_unresolved_bad_checks": results[best_combo],
                                       "worst_combination_bad_checks": max(v[1] for v in results.values())}
        self.naming["mirror_rules"] = dict(self.mirror_rule)
        return list(decisions.values())

    def _deps(self, p: PartRec) -> set:
        deps: set = set()
        for op in p.ops:
            for b in op.walk():
                for key in ("ToolParts",):
                    if b.has(key):
                        deps.update(b[key])
                if b.has("PlaceHolderOpnId"):
                    deps.add(self.op_owner.get(b["PlaceHolderOpnId"], (None,))[0])
            if op.has("PlaceHolderOpnId"):
                deps.add(self.op_owner.get(op["PlaceHolderOpnId"], (None,))[0])
            cs = op.get("ReferenceCoordSystemID")
            if cs in self.cs_ops and self.cs_ops[cs].has("PlaceHolderOperationID"):
                deps.add(self.op_owner.get(self.cs_ops[cs]["PlaceHolderOperationID"], (None,))[0])
        deps.discard(None)
        deps.discard(p.pid)
        return deps

    @staticmethod
    def _closure(pid: int, users: dict) -> set:
        out, stack = {pid}, [pid]
        while stack:
            for u in users.get(stack.pop(), ()):
                if u not in out:
                    out.add(u)
                    stack.append(u)
        return out

    def _op_DuplicateBodyAlongLine(self, part, op, prev):
        body, src_op = self._clone_source(op)
        sprm = src_op.block("DuplicateToAlongLineParameters")
        L = self.ev.length_mm
        F = self.cs(src_op["ReferenceCoordSystemID"])
        k = op.blocks()[0]["WhichClone"] + 1
        if abs(L(sprm["ZComponent"])) > 0:
            raise Unresolved("duplicate along z")
        v = F.vec((k * L(sprm["XComponent"]), k * L(sprm["YComponent"])))
        return self._renumber(body.transform(Affine.translation(*v)), op)

    def _op_DuplicateAroundAxis(self, part, op, prev):
        return prev

    _op_DuplicateMirror = _op_DuplicateAroundAxis
    _op_DuplicateAlongLine = _op_DuplicateAroundAxis
    _op_CloneTo = _op_DuplicateAroundAxis

    # local ops ---------------------------------------------------------------
    def _op_MoveEdges(self, part, op, prev):
        prm = op.block("MoveEdgesParameters")
        if not prm.get("MoveAlongNormalFlag"):
            raise Unresolved("MoveEdges along a vector")
        d = self.ev.length_mm(prm["OffsetDistance"])
        body = prev
        for eid in prm["EdgesToMove"]:
            li, i = self._one(g.find_edge(body, eid), "edge", eid)
            lp = body.loops()[li]
            body = g.replace_loop(body, li, g.offset_edge(lp, i, d))
        return body

    def _op_Fillet(self, part, op, prev):
        prm = op.block("FilletParameters")
        if prm["Edges"]:
            raise Unresolved("edge fillet in 2-D")
        if self.ev.length_mm(prm.get("Setback", "0mm")) != 0:
            raise Unresolved("fillet with setback")
        r = self.ev.length_mm(prm["Radius"])
        body = prev
        if not isinstance(body, Sheet):
            raise Unresolved("fillet of a wire")
        for vid in prm["Vertices"]:
            li, i = self._one(g.find_vertex(body, vid), "vertex", vid)
            lp = body.loops()[li]
            body = g.replace_loop(body, li, g.fillet_vertex(lp, i, r, tag=f"f{op['ID']}.v{vid}"))
        return body

    def _op_SweepAlongVector(self, part, op, prev):
        prm = op.block("VectorSweepParameters")
        L = self.ev.length_mm
        vx, vy, vz = L(prm["SweepVectorX"]), L(prm["SweepVectorY"]), L(prm["SweepVectorZ"])
        z0 = self.offplane.get(part.pid, 0.0)
        if vx or vy or self.ev.angle_rad(prm.get("DraftAngle", "0deg")) != 0 or not isinstance(prev, Sheet):
            raise Unresolved(f"sweep along ({vx}, {vy}, {vz}) mm with draft: not a prism along z")
        if not (min(z0, z0 + vz) <= 0.0 <= max(z0, z0 + vz)):
            raise Unresolved(f"swept prism z in [{z0}, {z0 + vz}] mm does not reach the model plane z=0")
        # the prism's cross-section with the model plane (here its end face) is the profile itself
        self.offplane.pop(part.pid, None)
        self.notes.append(f"{part.name}: profile at z={z0} mm swept by {vz} mm along z; the prism touches z=0 and "
                          "is used through its z=0 section (= the profile) as a 2-D tool")
        return prev

    # ------------------------------------------------------------------
    # entity naming
    def _cached_entities(self, op: Block) -> tuple[list, list, list]:
        helper = _identity_helper(op)
        edges, verts, faces = [], [], []
        if helper is None:
            return edges, verts, faces
        nb = helper.opt_block("NewEdges")
        for e in (nb.blocks("Edge") if nb else []):
            geo = e.opt_block("EdgeGeomTopol")
            src = geo.block("EdgeGeometry") if geo and geo.opt_block("EdgeGeometry") else e
            tv = [f.value()[:2] for b in src.walk() for f in b.fields("TolVt")]
            mids = [f.value()[:2] for b in src.walk() for f in b.fields("EdgeMidPoint")]
            edges.append((e["ID"], tv, mids[0] if mids else None))
        vb = helper.opt_block("NewVertices")
        for v in (vb.blocks("Vertex") if vb else []):
            pos = [f.value()[:2] for b in v.walk() for f in b.fields("VtPos")]
            if pos:
                verts.append((v["ID"], pos[0]))
        fb = helper.opt_block("NewFaces")
        for f in (fb.blocks("Face") if fb else []):
            area = [fl.value() for b in f.walk() for fl in b.fields("Area")]
            tv = [fl.value()[:2] for b in f.walk() for fl in b.fields("TolVt")]
            faces.append((f["ID"], area[0] if area else None, tv))
        return edges, verts, faces

    def _attach_cached(self, part, op, body):
        edges, verts, faces = self._cached_entities(op)
        if not (edges or verts or faces) or isinstance(body, Wire) and not (edges or verts):
            if faces and isinstance(body, Sheet):
                self._check_faces(part, op, body, faces)
            return body
        names: dict = {}
        segs_all = [(li, i, s) for li, lp in enumerate(body.loops()) for i, s in enumerate(lp)]
        new_ids: dict = {}
        for eid, tv, mid in edges:
            best = None
            for li, i, s in segs_all:
                ends = [s.start(), s.end()]
                if tv and len(tv) == 2:
                    dev = min(max(dist(ends[0], tv[0]), dist(ends[1], tv[1])),
                              max(dist(ends[0], tv[1]), dist(ends[1], tv[0])))
                else:
                    dev = 0.0
                if mid is not None:
                    dev = max(dev, _dist_to_seg(s, tuple(mid)))
                if best is None or dev < best[0]:
                    best = (dev, li, i, s)
            self.checks.append(Check(op["ID"], part.name, "edge", eid, tuple(mid or tv[0]),
                                     tuple(best[3].at(0.5)) if best else None, best[0] if best else math.inf))
            if best and best[0] < POS_TOL:
                new_ids[(best[1], best[2])] = eid
                names[eid] = ("e", best[3].tag)
        verts_found: dict = {}
        for vid, pos in verts:
            best = None
            for li, lp in enumerate(body.loops()):
                for i, s in enumerate(lp):
                    dev = dist(s.start(), pos)
                    if best is None or dev < best[0]:
                        best = (dev, li, i)
                    if i == len(lp) - 1:
                        dev = dist(s.end(), pos)
                        if dev < best[0]:
                            best = (dev, li, i + 1)
            lp = body.loops()[best[1]]
            mine = lp[best[2]].start() if best[2] < len(lp) else lp[-1].end()
            self.checks.append(Check(op["ID"], part.name, "vertex", vid, tuple(pos), tuple(mine), best[0]))
            if best[0] < POS_TOL:
                verts_found[(best[1], best[2])] = vid
                names[vid] = ("v", _vertex_sig(lp, best[2]))
        if faces and isinstance(body, Sheet):
            self._check_faces(part, op, body, faces)
        self.naming[op["ID"]] = names
        return _apply_ids(body, new_ids, verts_found)

    def _check_faces(self, part, op, body, faces):
        for fid, area, tv in faces:
            if area is None:
                continue
            best = min(body.faces, key=lambda f: abs(g.face_area(f) - area))
            mine = g.face_area(best)
            self.checks.append(Check(op["ID"], part.name, "face-area", fid, (area,), (mine,), abs(mine - area)))
            for p in tv:
                dev = min(dist(q, p) for q in g.body_vertices(Sheet((best,))))
                self.checks.append(Check(op["ID"], part.name, "face-vertex", fid, tuple(p), None, dev))

    def _attach_named(self, op, body):
        names = self.naming.get(op["ID"])
        if not names:
            return body
        new_ids, verts = {}, {}
        for aid, (kind, sig) in names.items():
            for li, lp in enumerate(body.loops()):
                for i, s in enumerate(lp):
                    if kind == "e" and s.tag == sig:
                        new_ids[(li, i)] = aid
                    if kind == "v" and _vertex_sig(lp, i) == sig:
                        verts[(li, i)] = aid
            if kind == "e" and aid not in new_ids.values():
                self.notes.append(f"op {op['ID']}: edge {aid} (tag {sig}) not found in variant")
            if kind == "v" and aid not in verts.values():
                self.notes.append(f"op {op['ID']}: vertex {aid} (between {sig}) not found in variant")
        return _apply_ids(body, new_ids, verts)


# ---------------------------------------------------------------------------
def _frame(o, x, y, x_driven: bool) -> Affine:
    """Local->global frame from origin, x direction and y direction (orthonormalised)."""
    def nrm(v):
        L = math.hypot(*v)
        if L == 0:
            raise Unresolved("zero-length CS axis")
        return (v[0] / L, v[1] / L)
    x, y = nrm(x), nrm(y)
    if x_driven:
        dot = x[0] * y[0] + x[1] * y[1]
        y = nrm((y[0] - dot * x[0], y[1] - dot * x[1]))
    else:
        dot = x[0] * y[0] + x[1] * y[1]
        x = nrm((x[0] - dot * y[0], x[1] - dot * y[1]))
    return Affine(x[0], y[0], x[1], y[1], o[0], o[1])


def _arc3(p0, p1, p2) -> Arc:
    ax, ay = p0
    bx, by = p1
    cx, cy = p2
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-14:
        raise Unresolved("collinear 3-point arc")
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    c = (ux, uy)
    r = dist(c, p0)
    a0 = math.atan2(ay - uy, ax - ux)
    a1 = math.atan2(by - uy, bx - ux)
    a2 = math.atan2(cy - uy, cx - ux)
    ccw = ((a1 - a0) % g.TWO_PI) < ((a2 - a0) % g.TWO_PI)
    da = (a2 - a0) % g.TWO_PI if ccw else -((a0 - a2) % g.TWO_PI)
    return Arc(tag="", c=c, r=r, a0=a0, da=da)


def _dist_to_seg(s, p) -> float:
    u = min(max(s.project(p), 0.0), 1.0)
    return dist(s.at(u), p)


def _vertex_sig(lp, i) -> tuple:
    if i >= len(lp):
        return (lp[-1].tag, None)
    return (lp[i - 1].tag, lp[i].tag)


def _map_segments(body, fn):
    if isinstance(body, Wire):
        return Wire(tuple(tuple(fn(s) for s in ch) for ch in body.chains))
    return Sheet(tuple(Face(tuple(fn(s) for s in f.outer), tuple(tuple(fn(s) for s in h) for h in f.holes))
                       for f in body.faces))


def _apply_ids(body, new_ids: dict, verts: dict):
    loops = [list(lp) for lp in body.loops()]
    for (li, i), eid in new_ids.items():
        loops[li][i] = loops[li][i].with_meta(ids=frozenset({eid}))
    for (li, i), vid in verts.items():
        lp = loops[li]
        if i < len(lp):
            lp[i] = lp[i].with_meta(v0=lp[i].v0 | {vid})
            prev = i - 1
            if prev >= 0 or isinstance(body, Sheet):
                lp[prev] = lp[prev].with_meta(v1=lp[prev].v1 | {vid})
        else:
            lp[-1] = lp[-1].with_meta(v1=lp[-1].v1 | {vid})
    out = body
    for li, lp in enumerate(loops):
        out = g.replace_loop(out, li, tuple(lp))
    return out


def _unique_tags(sh: Sheet) -> Sheet:
    """Pieces of one source segment get '~k' suffixes ordered by position (deterministic)."""
    groups: dict = {}
    for f in sh.faces:
        for lp in f.loops():
            for s in lp:
                groups.setdefault(s.tag, []).append(s)
    rename: dict = {}
    for tag, segs in groups.items():
        if len(segs) > 1:
            order = sorted(segs, key=lambda s: (round(s.at(0.5)[0], 6), round(s.at(0.5)[1], 6)))
            for k, s in enumerate(order):
                rename[id(s)] = f"{tag}~{k}"
    if not rename:
        return sh
    return _map_segments(sh, lambda s: s.with_meta(tag=rename.get(id(s), s.tag)))
