"""Command line interface.

    sipm em-sweep   CONFIG.toml [--stop-after N] [--retry-failed]
    sipm mech-sweep CONFIG.toml [--stop-after N] [--retry-failed]
    sipm sweep      CONFIG.toml            (physics taken from the file)
    sipm plot       CONFIG.toml | RESULTS.csv [--kind ...] [--out DIR]
    sipm params     [--markdown] [--topology T]
    sipm geometry   [--topology T] [--set k=v ...] --out FILE.png
    sipm deck       [--topology T] [--set k=v ...] --out FILE.inp
    sipm validate   ring|ring-liftoff|ring-multilayer|liftoff|mesh|open-circuit [...]
    sipm composite  --fibre F [--matrix M] [--vf V] [--angle DEG | --layers SPEC] [--model chamis|halpin_tsai|rom]
                    [--toml]
    sipm winding    --tension MPa[,MPa...] [--sublayers M] [--topology T] [--set k=v ...]

A sweep configuration (TOML):

    [sweep]
    name = "em_bridged_demo"          # folder under workdir
    physics = "em"                    # em | mech | acloss
    topology = "sleeved_bridged"
    workdir = "runs"                  # optional, else $SIPM_WORKDIR or <tempdir>/sipm_runs

    [base]                            # applied to every point
    em_span = 60.0

    [grid]                            # full factorial over these lists
    sleeve_thickness = [1.0, 2.0]
    gamma = [110, 120, 130, 140, 150]

    [[points]]                        # optional explicit points (added to the grid)
    sleeve_thickness = 1.5
    gamma = 125

    [ansys]                           # optional Config fields
    aedt_version = "2024.1"
    mapdl_cores = 2
    min_free_ram_gb = 4.0             # memory gate before every point
    ram_wait_timeout_s = 1800         # then stop cleanly (resumable)
    max_solver_ram_gb = 16.0          # optional: kill a solver above this RSS
    aedt_restart_every = 10           # release AEDT every N EM points

    [plot]                            # optional, used by `sipm plot CONFIG.toml`
    sigma_allow = 400.0
    interference_fraction = 0.2
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib


def load_config(path: str | Path) -> dict:
    with open(path, "rb") as fh:
        cfg = tomllib.load(fh)
    if "sweep" not in cfg or "name" not in cfg["sweep"]:
        raise ValueError(f"{path}: [sweep] section with 'name' is required")
    return cfg


def sweep_from_config(cfg: dict, physics: str | None = None):
    from .config import Config
    from .sweep import Sweep, grid

    s = cfg["sweep"]
    phys = physics or s.get("physics")
    if s.get("physics") and physics and s["physics"] != physics:
        raise ValueError(f"config physics is '{s['physics']}' but command asks for '{physics}'")
    table = grid(**cfg["grid"]) if cfg.get("grid") else []
    table += [dict(p) for p in cfg.get("points", [])]
    if not table:
        table = [{}]
    c = Config(**cfg.get("ansys", {}))
    if s.get("workdir"):
        c.workdir = Path(s["workdir"])
    return Sweep(s["name"], phys, table, base=cfg.get("base", {}), topology=s.get("topology"), cfg=c)


def _parse_sets(items):
    import dataclasses

    from .params import Params, coerce

    str_fields = {f.name for f in dataclasses.fields(Params) if isinstance(f.default, str)}
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        try:
            v = json.loads(v)
        except json.JSONDecodeError:
            pass
        k = k.strip()
        if k in str_fields:   # e.g. sleeve_layer_prestress=100 parses as a number: keep it a string
            v = coerce(k, v)
        out[k] = v
    return out


def _composite(a) -> int:
    """`sipm composite`: derived sleeve constants (through Params, so exactly what a sweep uses)."""
    import dataclasses

    from .composite import FIBRES, SLEEVE_KEYS
    from .params import Params, default_params

    spec = dict(sleeve_fibre=a.fibre, sleeve_matrix=a.matrix, sleeve_vf=a.vf, sleeve_micromechanics=a.model,
                sleeve_wind_angle=a.angle)
    if a.layers:
        spec["sleeve_layers"] = a.layers
    if a.fibre not in FIBRES:
        raise SystemExit(f"sipm composite: --fibre must be one of {tuple(FIBRES)}, got {a.fibre!r}")
    try:
        p = default_params("sleeved_bridged", **spec)
    except ValueError as e:
        raise SystemExit(f"sipm composite: {e}") from None
    units = {f.name: f.metadata.get("unit", "") for f in dataclasses.fields(Params)}
    if a.layers:
        return _composite_layers(a, p, spec, units)
    print(f"{a.fibre} / {a.matrix}, vf = {a.vf:g}, +/-{a.angle:g} deg from hoop, {a.model} micromechanics")
    for k in SLEEVE_KEYS:
        print(f"  {k:<24s} = {getattr(p, k):.6g} {units[k]}")
    if a.toml:
        print()
        print("# [base] lines: composite spec (sipm derives the sleeve constants from it at run time)")
        for k, v in spec.items():
            print(f'{k} = "{v}"' if isinstance(v, str) else f"{k} = {v!r}")
        print("# or the derived constants entered directly (leave sleeve_fibre unset); the failure indices")
        print("# then treat the sleeve as one hoop-wound ply, which is right only for a 0 deg winding")
        for k in SLEEVE_KEYS:
            print(f"# {k} = {getattr(p, k)!r}")
    return 0


def _composite_layers(a, p, spec: dict, units: dict) -> int:
    """`sipm composite --layers`: constants of every layer and of the stack (membrane homogenisation)."""
    from .composite import SLEEVE_KEYS, sleeve_params

    angles, fracs = p.sleeve_layer_angles(), p.sleeve_layer_fractions()
    consts = p.sleeve_layer_constants()
    stack = sleeve_params(p.sleeve_laminate())
    print(f"{a.fibre} / {a.matrix}, vf = {a.vf:g}, {a.model} micromechanics, {len(angles)} layers from the bore: "
          f"{a.layers}")
    head = "".join(f"{f'layer {i + 1}':>14s}" for i in range(len(angles))) + f"{'stack':>14s}"
    print(f"  {'':<24s} {'':<6s}{head}")
    print(f"  {'winding angle':<24s} {'deg':<6s}" + "".join(f"{f'+/-{x:g}':>14s}" for x in angles) + f"{'':>14s}")
    print(f"  {'thickness fraction':<24s} {'-':<6s}" + "".join(f"{x:14.4g}" for x in fracs) + f"{1.0:14.4g}")
    for k in SLEEVE_KEYS:
        print(f"  {k:<24s} {units[k]:<6s}" + "".join(f"{c[k]:14.6g}" for c in consts) + f"{stack[k]:14.6g}")
    print("  (stack = membrane homogenisation of the layers: the values the thirteen sleeve parameters take;")
    print("   the deck and the failure evaluation use each layer's own constants)")
    if a.toml:
        print()
        print("# [base] lines: layered composite spec (sipm derives every layer's constants at run time)")
        for k, v in spec.items():
            print(f'{k} = "{v}"' if isinstance(v, str) else f"{k} = {v!r}")
    return 0


def _winding(a) -> int:
    """`sipm winding`: equivalent sleeve_layer_prestress of a layer-by-layer winding (closed form)."""
    from .params import default_params
    from .validation import winding_prestress

    try:
        p = default_params(a.topology, **_parse_sets(a.set))
        tension = [float(x) for x in a.tension.split(",")]
        r = winding_prestress(p, tension[0] if len(tension) == 1 else tension, a.sublayers)
    except ValueError as e:
        raise SystemExit(f"sipm winding: {e}") from None
    radii, seq = r["radii"], r["sequential"]
    angles, consts = p.sleeve_layer_angles(), p.sleeve_layer_constants()
    print(f"{a.topology}: lamination as a solid disc of radius {radii[0]:g} mm, {len(angles)} sleeve layer(s) to "
          f"{radii[-1]:g} mm, wound one after the other"
          + (f" ({a.sublayers} sub-rings per layer)" if a.sublayers > 1 else ""))
    print(f"  {'layer':<6s}{'angle':>8s}{'r bore':>9s}{'r OD':>9s}{'E_hoop':>9s}{'winding':>10s}{'built-up':>10s}"
          f"{'equivalent':>12s}")
    print(f"  {'':<6s}{'deg':>8s}{'mm':>9s}{'mm':>9s}{'GPa':>9s}{'MPa':>10s}{'mean MPa':>10s}{'prestress':>12s}")
    tw = tension * len(angles) if len(tension) == 1 else tension
    for i, ang in enumerate(angles):
        print(f"  {i + 1:<6d}{ang:8g}{radii[i]:9.4f}{radii[i + 1]:9.4f}{consts[i]['sleeve_E_hoop']:9.4g}"
              f"{tw[i]:10.4g}{seq['layers'][i]['sh_mean']:10.4g}{r['prestress'][i]:12.5g}")
    print(f"  contact pressure on the lamination from the winding alone: {-seq['sr'][0]:.4g} MPa")
    print(f"  the equivalent prestresses reproduce every layer's hoop force and the interface pressures of the "
          f"build-up; hoop stresses within a layer differ by up to {r['max_hoop_dev']:.3g} MPa")
    print()
    print(f'sleeve_layer_prestress = "{r["spec"]}"')
    return 0


def _save_and_close(plot: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
    """Call a :mod:`sipm.post` function (which saves to ``out=``) and close its figure, so
    plotting after every sweep does not accumulate open figures."""
    import matplotlib.pyplot as plt

    fig = plot(*args, **kwargs)
    plt.close(fig)


def plot_results(df: pd.DataFrame, outdir: Path, physics: str, plot_cfg: dict | None = None,
                 prefix: str = "") -> list[Path]:
    """Standard figure set for a results table; returns the files written (figures are closed)."""
    from . import post

    plot_cfg = plot_cfg or {}
    outdir.mkdir(parents=True, exist_ok=True)
    files = []
    ok = df[df["status"] == "ok"] if "status" in df else df
    if ok.empty:
        return files
    note = post.conditions_note(ok, physics)  # shared topology / current / speed, in words
    if physics == "em":
        swept = [c for c in ("gamma", "sleeve_thickness") + post.GEOMETRY_SWEEP_COLUMNS
                 if c in ok and ok[c].nunique() > 1]
        if swept:
            f = outdir / f"{prefix}torque_scatter.png"
            _save_and_close(post.torque_scatter, df, swept, reference=plot_cfg.get("reference_torque"), out=f)
            files.append(f)
        if "gamma" in ok and ok["gamma"].nunique() > 1:
            grp = [c for c in ("sleeve_thickness",) + post.GEOMETRY_SWEEP_COLUMNS
                   if c in ok and ok[c].nunique() > 1] or ["topology"]
            f = outdir / f"{prefix}torque_vs_gamma.png"
            _save_and_close(post.mtpa_lines, df, grp, note=note, out=f)
            files.append(f)
    elif physics == "mech":
        # with load cases (mech_cases) also the worst case over them
        for q in ("svm_rotor_max", "sh_sleeve_max", "svm_lam_max__worst", "sh_sleeve_max__worst"):
            if q in ok:
                f = outdir / f"{prefix}{q}_vs_interference.png"
                _save_and_close(post.stress_lines, df, q, reference=plot_cfg.get(f"reference_{q}"), note=note, out=f)
                files.append(f)
        if "sigma_allow" in plot_cfg and ok["sleeve_thickness"].nunique() > 1 and ok["interference"].nunique() > 1:
            f = outdir / f"{prefix}stress_map.png"
            _save_and_close(post.stress_map, df, plot_cfg["sigma_allow"], plot_cfg.get("interference_fraction", 0.2),
                            note=note, out=f)
            files.append(f)
    elif physics == "acloss":
        f = outdir / f"{prefix}kac.png"
        _save_and_close(post.kac_lines, df, note=note, out=f)
        files.append(f)
    return files


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="sipm", description="Sleeved IPM rotor Maxwell/MAPDL sweeps")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("sweep", "em-sweep", "mech-sweep", "acloss-sweep"):
        s = sub.add_parser(name)
        s.add_argument("config")
        s.add_argument("--stop-after", type=int, default=None)
        s.add_argument("--retry-failed", action="store_true")
        s.add_argument("--no-plots", action="store_true")
        s.add_argument("--min-free-ram-gb", type=float, default=None,
                       help="wait before each point until this much RAM is free")
        s.add_argument("--ram-wait-timeout", type=float, default=None,
                       help="seconds to wait for free RAM before stopping the sweep (resumable)")
        s.add_argument("--max-solver-ram-gb", type=float, default=None,
                       help="terminate a solver process tree above this RSS; the point is recorded as failed")
        s.add_argument("--aedt-restart-every", type=int, default=None,
                       help="release and restart AEDT every N EM points (0 = never)")
    s = sub.add_parser("plot")
    s.add_argument("source", help="sweep config (.toml) or results .csv")
    s.add_argument("--physics", choices=["em", "mech", "acloss"])
    s.add_argument("--out", default=None)
    s.add_argument("--sigma-allow", type=float, default=None)
    s.add_argument("--interference-fraction", type=float, default=0.2)
    s = sub.add_parser("params")
    s.add_argument("--markdown", action="store_true")
    s.add_argument("--json", action="store_true")
    s.add_argument("--topology", default="sleeved_bridged")
    s.add_argument("--set", nargs="*", default=[])
    for name in ("geometry", "deck"):
        s = sub.add_parser(name)
        s.add_argument("--topology", default="sleeved_bridged")
        s.add_argument("--set", nargs="*", default=[])
        s.add_argument("--out", required=True)
    s = sub.add_parser("composite", help="sleeve constants from fibre, matrix, vf and winding angle (no Ansys)")
    s.add_argument("--fibre", required=True, help="fibre class of sipm.composite.FIBRES")
    s.add_argument("--matrix", default="epoxy", help="matrix class of sipm.composite.MATRICES")
    s.add_argument("--vf", type=float, default=0.55, help="fibre volume fraction")
    s.add_argument("--angle", type=float, default=0.0, help="winding angle from hoop (deg), balanced +/-")
    s.add_argument("--model", default="chamis", help="micromechanics: chamis, halpin_tsai or rom")
    s.add_argument("--layers", default="",
                   help="layered sleeve, 'angle:weight' items from the bore, e.g. '0:2,20:1,0:1' (--angle must be 0)")
    s.add_argument("--toml", action="store_true", help="also print the lines for a sweep's [base] table")
    s = sub.add_parser("winding", help="equivalent sleeve_layer_prestress of a layer-by-layer winding (no Ansys)")
    s.add_argument("--tension", required=True,
                   help="winding stress in MPa, one value or one per layer (mean hoop stress of a layer as it is laid)")
    s.add_argument("--sublayers", type=int, default=1, help="lay every layer as this many sub-rings in turn")
    s.add_argument("--topology", default="sleeved_bridged")
    s.add_argument("--set", nargs="*", default=[])
    s = sub.add_parser("validate")
    s.add_argument("what", choices=["ring", "ring-liftoff", "ring-multilayer", "liftoff", "mesh", "open-circuit"])
    s.add_argument("--case", default="orthotropic",
                   help="ring: isotropic|orthotropic|thermal; liftoff: cold|hot (default cold)")
    s.add_argument("--topology", default="sleeved_bridged")
    s.add_argument("--workdir", default=None)
    a = ap.parse_args(argv)

    if a.cmd in ("sweep", "em-sweep", "mech-sweep", "acloss-sweep"):
        phys = {"em-sweep": "em", "mech-sweep": "mech", "acloss-sweep": "acloss"}.get(a.cmd)
        cfg = load_config(a.config)
        sw = sweep_from_config(cfg, phys)
        for opt, fld in (("min_free_ram_gb", "min_free_ram_gb"), ("ram_wait_timeout", "ram_wait_timeout_s"),
                         ("max_solver_ram_gb", "max_solver_ram_gb"), ("aedt_restart_every", "aedt_restart_every")):
            v = getattr(a, opt)
            if v is not None:
                setattr(sw.cfg, fld, v)
        df = sw.run(stop_after=a.stop_after, retry_failed=a.retry_failed)
        print(f"results: {sw.dir / 'results.csv'}")
        print(df["status"].value_counts().to_string())
        if not a.no_plots:
            for f in plot_results(df, sw.dir / "figures", sw.physics, cfg.get("plot")):
                print(f"figure: {f}")
        if sw.stopped_reason:
            print(f"STOPPED: {sw.stopped_reason}")
            return 3
        return 0 if (df["status"] != "failed").all() else 2
    if a.cmd == "plot":
        import pandas as pd

        src = Path(a.source)
        plot_cfg = {}
        if src.suffix == ".toml":
            cfg = load_config(src)
            sw = sweep_from_config(cfg)
            df, phys, out = sw.results(), sw.physics, Path(a.out) if a.out else sw.dir / "figures"
            plot_cfg = cfg.get("plot", {})
        else:
            df = pd.read_csv(src)
            phys = a.physics or ("em" if "torque_avg" in df else "acloss" if "kac" in df else "mech")
            out = Path(a.out) if a.out else src.parent / "figures"
        if a.sigma_allow is not None:
            plot_cfg["sigma_allow"] = a.sigma_allow
            plot_cfg["interference_fraction"] = a.interference_fraction
        for f in plot_results(df, out, phys, plot_cfg):
            print(f"figure: {f}")
        return 0
    if a.cmd == "params":
        from .params import default_params, parameter_markdown

        if a.markdown:
            print(parameter_markdown())
        elif a.json:
            print(json.dumps(default_params(a.topology, **_parse_sets(a.set)).to_dict()))
        else:
            for k, v in default_params(a.topology, **_parse_sets(a.set)).to_dict().items():
                print(f"{k} = {v!r}")
        return 0
    if a.cmd == "composite":
        return _composite(a)
    if a.cmd == "winding":
        return _winding(a)
    if a.cmd in ("geometry", "deck"):
        from .params import default_params

        p = default_params(a.topology, **_parse_sets(a.set))
        if a.cmd == "geometry":
            from .post import geometry_plot

            _save_and_close(geometry_plot, p, out=a.out)
        else:
            from .mapdl import copy_deck_only

            copy_deck_only(p, a.out)
        print(a.out)
        return 0
    if a.cmd == "validate":
        from .config import Config

        c = Config(workdir=a.workdir) if a.workdir else Config()
        if a.what == "ring":
            from .validation import ring_check

            r = ring_check(a.case, cfg=c)
            for row in r["rows"]:
                print(f"{row['quantity']:<40s} FE {row['fe']:12.5f}  closed form {row['closed_form']:12.5f}  "
                      f"error {100 * row['rel_error']:+.3f} %")
        elif a.what == "ring-liftoff":
            from .validation import ring_liftoff_check

            for row in ring_liftoff_check(cfg=c)["rows"]:
                print(f"lift-off {row['label']:<4s} ({row['temperature']:g} degC)  FE {row['fe']:10.1f} rpm  "
                      f"closed form {row['closed_form']:10.1f} rpm  error {100 * row['rel_error']:+.4f} %  "
                      f"[{row['basis']}]")
        elif a.what == "ring-multilayer":
            from .validation import multilayer_ring_check

            r = multilayer_ring_check(cfg=c)
            for row in r["rows"]:
                print(f"{row['quantity']:<42s} FE {row['fe']:12.5f}  closed form {row['closed_form']:12.5f}  "
                      f"error {100 * row['rel_error']:+.3f} %")
            for row in r["liftoff"]:
                print(f"lift-off {row['label']:<4s} ({row['temperature']:g} degC)  FE {row['fe']:10.1f} rpm  "
                      f"closed form {row['closed_form']:10.1f} rpm  error {100 * row['rel_error']:+.4f} %")
        elif a.what == "liftoff":
            from .params import default_params
            from .validation import liftoff_check

            label = a.case if a.case in ("cold", "hot") else "cold"
            r = liftoff_check(default_params(a.topology), label, cfg=c)
            print(f"{a.topology} lift-off ({label}, {r['temperature']:g} degC): {r['liftoff_rpm']:.1f} rpm "
                  f"[{r['basis']}]")
            for tag, d in r["direct"].items():
                print(f"  direct run {tag:<5s} {d['rpm']:10.1f} rpm: contact p min {d['contact_p_min']:.4f} MPa, "
                      f"open fraction {d['sleeve_open_fraction']:.4f}")
        elif a.what == "mesh":
            from .params import default_params
            from .validation import mesh_study

            for row in mesh_study(default_params(a.topology), cfg=c):
                print(row)
        else:
            from .maxwell import MaxwellSession
            from .params import default_params
            from .validation import open_circuit_check

            p = default_params(a.topology)
            with MaxwellSession(c) as s:
                d = s.open_circuit(p)
            print(json.dumps(open_circuit_check(d, p), indent=1, default=float))
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
