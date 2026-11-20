"""Named stress sections (sipm.geometry.locations)."""

import math

import numpy as np
import pytest

from sipm.geometry import build_rotor
from sipm.geometry.checks import in_steel
from sipm.geometry.locations import QHOLE_BAND, section_peaks, stress_sections
from sipm.params import default_params

NAMES = {
    "sleeved_bridged": ["bridge1", "bridge2", "rib1", "web2", "qhole"],
    "sleeved_bridgeless": ["lip1_pole", "lip1_q", "lip2_pole", "lip2_q", "qhole"],
}


def _sample(sec, n=400):
    """Points spread over a section (for membership checks)."""
    rng = np.random.default_rng(3)
    u, v = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    if sec.kind == "sector":
        r = sec.r0 + 1e-3 + v * (sec.r1 - sec.r0 - 2e-3)
        a = sec.a0 + 1e-4 + u * (sec.a1 - sec.a0 - 2e-4)
        return np.c_[sec.c[0] + r * np.cos(a), sec.c[1] + r * np.sin(a)]
    o, d, nn = map(np.asarray, (sec.o, sec.d, sec.n))
    s = sec.s0 + u * (sec.s1 - sec.s0)
    t = 1e-3 + v * (sec.depth - 2e-3)
    return o + s[:, None] * d + t[:, None] * nn


@pytest.mark.parametrize("topo", list(NAMES))
def test_sections_named_and_inside_the_steel(topo):
    rp = build_rotor(default_params(topo))
    secs = stress_sections(rp)
    assert list(dict.fromkeys(s.name for s in secs)) == NAMES[topo]
    for s in secs:
        pts = _sample(s)
        assert s.contains(pts).all()
        assert in_steel(rp, pts).all(), s.name           # every section lies in the lamination
        assert (pts[:, 1] >= -1e-9).all()                # upper half pole (mechanical model)


def test_bridged_sections_follow_the_parameters():
    p = default_params("sleeved_bridged")
    rp = build_rotor(p)
    by = {s.name: s for s in stress_sections(rp)}
    for k in (1, 2):
        b = by[f"bridge{k}"]
        assert math.isclose(b.r0, rp.r_lam - getattr(p, f"pm{k}_bridge_thickness"), abs_tol=1e-9)
        assert math.isclose(b.r1, rp.r_lam, abs_tol=1e-9)
    assert math.isclose(by["rib1"].depth, p.pm1_rib_half_width + p.pm1_inner_offset, abs_tol=1e-9)
    assert math.isclose(by["web2"].depth, p.pm2_web_width)
    # the corner margin shortens each section by twice its value
    q = p.replace(stress_corner_margin=0.25)
    by2 = {s.name: s for s in stress_sections(build_rotor(q))}
    for n in ("bridge1", "bridge2", "rib1", "web2"):
        assert math.isclose(by2[n].span - by[n].span, 0.5, rel_tol=1e-6), n
    assert all(s.r1 - s.r0 == pytest.approx(QHOLE_BAND) for s in stress_sections(rp) if s.name == "qhole")


def test_section_peaks_pick_the_maximum_inside_each_section():
    rp = build_rotor(default_params("sleeved_bridged"))
    secs = stress_sections(rp)
    pts = np.vstack([_sample(s, 50) for s in secs])
    vals = np.arange(len(pts), dtype=float)
    out = section_peaks(secs, pts, vals)
    assert list(out) == NAMES["sleeved_bridged"]
    for name, pk in out.items():
        mask = np.zeros(len(pts), bool)
        for s in secs:
            if s.name == name:
                mask |= s.contains(pts)
        assert pk["value"] == vals[mask].max() and pk["nodes"] == mask.sum()
    with pytest.raises(Exception, match="no lamination node"):
        section_peaks(secs, np.array([[0.0, 0.0]]), np.array([1.0]))
