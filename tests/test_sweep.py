import json

import pandas as pd
import pytest

from sipm.config import Config
from sipm.params import default_params
from sipm.sweep import Sweep, grid, mtpa, point_key, resolve


class FakeSolver:
    """Deterministic stand-in for a solver; fails on a chosen interference value."""

    def __init__(self, fail_at=None):
        self.calls = []
        self.fail_at = fail_at

    def __call__(self, p, wd):
        self.calls.append(p)
        if self.fail_at is not None and abs(p.interference - self.fail_at) < 1e-12:
            raise RuntimeError("synthetic solver failure")
        return dict(svm_rotor_max=500 - 100 * p.interference * p.sleeve_thickness, sh_sleeve_max=1000 * p.interference,
                    torque_avg=400 - (p.gamma - 125) ** 2 / 50, _wave={"Time": [0, 1], "Moving1.Torque": [1, 2]})


def test_grid():
    g = grid(a=[1, 2], b=[3, 4, 5])
    assert len(g) == 6 and g[0] == {"a": 1, "b": 3} and g[-1] == {"a": 2, "b": 5}


def test_resolve_and_key():
    p = resolve({"sleeve_thickness": 2, "topology": "sleeved_bridgeless"})
    assert p.topology == "sleeved_bridgeless" and p.sleeve_thickness == 2.0
    assert point_key(p, "mech") != point_key(p, "em")
    assert point_key(p, "mech") == point_key(default_params("sleeved_bridgeless", sleeve_thickness=2.0), "mech")
    with pytest.raises(KeyError):
        resolve({"not_a_param": 1})


def test_cache_resume_and_failures(tmp_path):
    cfg = Config(workdir=tmp_path, min_free_ram_gb=0)  # gate tested in test_memory.py
    table = grid(sleeve_thickness=[1.0, 2.0], interference=[0.0, 0.1, 0.2])
    table.append({"sleeve_thickness": 50.0, "interference": 0.1})  # invalid design
    sw = Sweep("t", "mech", table, topology="sleeved_bridged", cfg=cfg)
    s1 = FakeSolver(fail_at=0.2)
    df = sw.run(solver=s1, stop_after=3, progress=None)  # simulated interruption
    assert len(s1.calls) == 3
    assert (df["status"] == "pending").sum() == 3
    assert (df["status"] == "invalid").sum() == 1
    # resume: only the remaining points are solved
    s2 = FakeSolver(fail_at=0.2)
    df = sw.run(solver=s2, progress=None)
    assert len(s2.calls) == 3
    assert (df["status"] == "ok").sum() == 4 and (df["status"] == "failed").sum() == 2
    assert df.loc[df.status == "failed", "message"].str.contains("synthetic").all()
    # third run: nothing to do; failures kept unless retried
    s3 = FakeSolver()
    sw.run(solver=s3, progress=None)
    assert s3.calls == []
    df = sw.run(solver=s3, retry_failed=True, progress=None)
    assert len(s3.calls) == 2 and (df["status"] == "ok").sum() == 6
    # results file and cache files
    assert (tmp_path / "t" / "results.csv").exists()
    rec = json.loads(next((tmp_path / "t" / "cache").glob("*.json")).read_text())
    assert {"key", "params", "outputs", "status", "wall_s"} <= set(rec)
    back = pd.read_csv(tmp_path / "t" / "results.csv")
    assert len(back) == len(table)


def test_waveforms_and_mtpa(tmp_path):
    cfg = Config(workdir=tmp_path, min_free_ram_gb=0)  # gate tested in test_memory.py
    sw = Sweep("em", "em", grid(sleeve_thickness=[1.0, 2.0], gamma=[115, 125, 135]), topology="sleeved_bridged",
               cfg=cfg)
    df = sw.run(solver=FakeSolver(), progress=None)
    assert len(sw.waveforms()) == 6
    best = mtpa(df)
    assert len(best) == 2 and (best["gamma"] == 125).all()


def test_bad_physics():
    with pytest.raises(ValueError):
        Sweep("x", "thermal", [{}])
