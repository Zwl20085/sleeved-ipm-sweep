from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sipm import post
from sipm.cli import main, plot_results
from sipm.maxwell import em_model_key, parse_report_csv
from sipm.params import default_params

ROOT = Path(__file__).resolve().parents[1]
FIX = Path(__file__).parent / "fixtures"
SAMPLE_EM = ROOT / "examples" / "data" / "sample_em_results.csv"
SAMPLE_MECH = ROOT / "examples" / "data" / "sample_mech_results.csv"


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


def _mech_df():
    rows = []
    for t in (0.5, 1.0, 1.5, 2.0):
        for d in (0.0, 0.1, 0.2, 0.3, 0.4):
            rows.append(dict(status="ok", sleeve_thickness=t, interference=d, svm_rotor_max=600 - 150 * d * t - 20 * t,
                             sh_sleeve_max=50 + 1500 * d))
    return pd.DataFrame(rows)


def test_stress_map_and_feasible(tmp_path):
    df = _mech_df()
    f = post.stress_map(df, sigma_allow=450, interference_fraction=0.2, out=tmp_path / "m.png")
    assert (tmp_path / "m.png").stat().st_size > 1000
    feas = post.feasible_points(df, 450, 0.2)
    assert ((feas["interference"] <= 0.2 * feas["sleeve_thickness"] + 1e-12) & (feas["svm_rotor_max"] <= 450)).all()
    assert len(feas) > 0
    f.clear()


def test_stress_map_needs_two_axes():
    with pytest.raises(ValueError):
        post.stress_map(_mech_df().query("sleeve_thickness == 1.0"), 400)


def test_line_and_scatter_plots(tmp_path):
    df = _mech_df()
    post.stress_lines(df, "svm_rotor_max", reference=420, out=tmp_path / "l.png")
    em = pd.DataFrame([dict(status="ok", gamma=g, sleeve_thickness=t, torque_avg=400 - (g - 125) ** 2 / 40 - 20 * t)
                       for g in (100, 115, 130, 145) for t in (1.0, 2.0)])
    post.torque_scatter(em, ["gamma", "sleeve_thickness"], reference=420, out=tmp_path / "s.png")
    post.mtpa_lines(em, "sleeve_thickness", out=tmp_path / "g.png")
    for n in ("l", "s", "g"):
        assert (tmp_path / f"{n}.png").exists()
    files = plot_results(em, tmp_path / "f", "em")
    assert len(files) == 2


def test_geometry_plot_and_cli(tmp_path):
    assert main(["geometry", "--topology", "sleeved_bridgeless", "--out", str(tmp_path / "g.png")]) == 0
    assert main(["deck", "--set", "interference=0.2", "--out", str(tmp_path / "d.inp")]) == 0
    assert "*SET,DELTA,0.2" in (tmp_path / "d.inp").read_text()
    assert main(["params", "--markdown"]) == 0


@pytest.mark.parametrize("sample", [SAMPLE_EM, SAMPLE_MECH])
def test_committed_samples_plot(tmp_path, sample):
    if not sample.exists():
        pytest.skip("sample not generated yet")
    assert main(["plot", str(sample), "--out", str(tmp_path), "--sigma-allow", "600"]) == 0
    assert any(tmp_path.glob("*.png"))


def _em_df():
    rows = [dict(status="ok", topology="sleeved_bridged", gamma=g, sleeve_thickness=t,
                 torque_avg=400 - (g - 125) ** 2 / 40 - 20 * t) for g in (100, 115, 130, 145) for t in (1.0, 2.0)]
    rows += [dict(status="ok", topology="conventional_baseline", gamma=g, sleeve_thickness=0.0,
                  torque_avg=330 - (g - 140) ** 2 / 40) for g in (120, 140, 160)]
    return pd.DataFrame(rows)


def _legend_texts(fig):
    return [t.get_text() for ax in fig.axes if ax.get_legend() for t in ax.get_legend().get_texts()]


def test_stress_map_legend_note_and_title():
    f = post.stress_map(_mech_df(), sigma_allow=450, interference_fraction=0.2, note="Test rotor, 20 000 rpm")
    texts = _legend_texts(f)
    assert texts == ["Stress = 450 MPa (allowable, user input)",
                     "Interference = 0.2 × sleeve thickness (limit)",
                     "Interference above limit (washed out)",
                     "Feasible region (both limits met)",
                     "Solved points"]
    leg = f.axes[0].get_legend()
    assert leg.get_title().get_text() == "Test rotor, 20 000 rpm"
    assert f.axes[0].get_title() == ""  # no title by default
    # the legend sits above the axes: it covers no data
    assert leg.get_window_extent().y0 >= f.axes[0].get_window_extent().y1
    g = post.stress_map(_mech_df(), 450, 0.2, title="custom title")
    assert g.axes[0].get_title() == "custom title"


def test_mtpa_lines_labels_and_note():
    f = post.mtpa_lines(_em_df(), ["topology", "sleeve_thickness"], note="420 A rms, 3000 rpm")
    texts = _legend_texts(f)
    assert texts == ["Conventional baseline (author's rotor, not in the package), no sleeve",
                     "Sleeved bridge-retained, 1.0 mm sleeve",
                     "Sleeved bridge-retained, 2.0 mm sleeve", "Maximum torque per curve (angle in deg)"]
    ax = f.axes[0]
    assert ax.get_title() == ""
    assert "420 A rms, 3000 rpm" in [t.get_text() for t in ax.texts]
    assert {"130°", "140°"} <= {t.get_text() for t in ax.texts}  # angle of each maximum
    # caller-supplied labels, one group column (scalar keys) and tuple keys
    g = post.mtpa_lines(_em_df().query("topology == 'sleeved_bridged'"), "sleeve_thickness",
                        labels={1.0: "thin sleeve", 2.0: "thick sleeve"})
    assert _legend_texts(g)[:2] == ["thin sleeve", "thick sleeve"]
    h = post.mtpa_lines(_em_df(), ["topology", "sleeve_thickness"], labels={("conventional_baseline", 0.0): "baseline"})
    assert _legend_texts(h)[0] == "baseline"
    assert _legend_texts(h)[1] == "Sleeved bridge-retained, 1.0 mm sleeve"  # default for keys not given
    # series differ by line style (not only colour) and use open markers
    series = [ln for ln in ax.get_lines() if not ln.get_label().startswith("_")]
    assert len(series) == 3
    assert len({str(ln._unscaled_dash_pattern) for ln in series}) == 3
    assert all(ln.get_markerfacecolor() == "none" for ln in series)


def test_series_label_defaults():
    assert post.series_label(["sleeve_thickness"], (0.0,)) == "No sleeve"
    assert post.series_label(["sleeve_thickness"], 1.5) == "1.5 mm sleeve"
    assert post.series_label(["topology", "rotor_slot_opening"], ("sleeved_bridgeless", 0.5)) == \
        "Sleeved bridgeless, rotor slot opening 0.5 mm"


def test_conditions_note():
    mech = pd.DataFrame(dict(topology=["sleeved_bridged"] * 2, mech_rpm=[20000.0] * 2))
    assert post.conditions_note(mech, "mech") == "Sleeved bridge-retained rotor, 20 000 rpm"
    em = pd.DataFrame(dict(topology=["conventional_baseline", "sleeved_bridged"], current_rms=[420.0, 420.0],
                           em_rpm=[np.nan, 3000.0]))
    assert post.conditions_note(em, "em") == "Phase current 420 A rms, 3000 rpm"
    assert post.conditions_note(pd.DataFrame(dict(x=[1])), "em") is None


def test_plots_leave_rcparams_unchanged(tmp_path):
    import matplotlib

    post.mtpa_lines(_em_df(), "sleeve_thickness")  # first call may switch the backend to Agg
    before = {k: v for k, v in matplotlib.rcParams.items() if k != "backend"}
    post.mtpa_lines(_em_df(), "sleeve_thickness", note="x")
    post.stress_map(_mech_df(), 450, 0.2)
    post.stress_lines(_mech_df(), reference=420, note="x")
    post.torque_scatter(_em_df(), ["gamma", "sleeve_thickness"], reference=420)
    post.kac_lines(pd.DataFrame(dict(topology=["sleeved_bridged"] * 2, em_rpm=[1000, 2000], kac=[1.1, 1.5])))
    post.torque_waveforms({"a": {"Time": [0, 1e-4], "Moving1.Torque": [1, 2]}}, 8, 3000)
    post.geometry_plot(default_params("sleeved_bridged"))
    post.topology_figure(out=tmp_path / "t.png")
    after = {k: v for k, v in matplotlib.rcParams.items() if k != "backend"}
    assert before == after
    assert (tmp_path / "t.png").stat().st_size > 1000
    import matplotlib.pyplot as plt

    plt.close("all")


def test_geometry_callouts_follow_parameters():
    from sipm._geomfig import callout_anchors
    from sipm.geometry import build_rotor, build_stator

    for topo, have, absent in [("sleeved_bridged", {"top_bridge", "centre_rib", "sleeve", "air_gap", "cavity",
                                                    "d_axis", "lightening_hole", "magnet_length"},
                                {"slot_opening", "lip"}),
                               ("sleeved_bridgeless", {"sleeve", "slot_opening", "lip", "air_gap"},
                                {"top_bridge", "centre_rib", "cavity"})]:
        p = default_params(topo)
        rp, sp = build_rotor(p), build_stator(p)
        a = callout_anchors(rp, sp)
        assert have <= set(a) and not (absent & set(a)), topo
        if "top_bridge" in a:  # inside the layer-2 outer bridge
            r = np.hypot(*a["top_bridge"])
            assert rp.r_lam - p.pm2_bridge_thickness < r < rp.r_lam
        if "sleeve" in a:
            assert rp.r_lam < np.hypot(*a["sleeve"]) < rp.r_sleeve_out
        if "slot_opening" in a:
            assert np.isclose(np.hypot(*a["slot_opening"]), rp.r_lam)
        assert p.r_sleeve_out < np.hypot(*a["air_gap"]) < p.r_stator_in
    # a thicker sleeve moves the sleeve callout end point with it
    p2 = default_params("sleeved_bridged", sleeve_thickness=3.0)
    a2 = callout_anchors(build_rotor(p2), build_stator(p2))
    assert build_rotor(p2).r_lam < np.hypot(*a2["sleeve"]) < p2.r_sleeve_out


def test_geometry_figure_legend():
    f = post.topology_figure()
    texts = [t.get_text() for leg in f.legends for t in leg.get_texts()]
    assert texts == ["Stator lamination", "Rotor lamination", "Permanent magnet (arrow: magnetisation direction)",
                     "Conductors (coil region)", "Carbon-fibre sleeve", "Air (magnet pocket / flux barrier)"]
    assert [ax.get_title().replace("\n", " ") for ax in f.axes] == [
        "(a) Sleeved, bridge-retained", "(b) Sleeved, bridgeless"]
    notes = [t.get_text() for t in f.texts]
    assert notes == ["48 slots / 8 poles, one pole (45°) shown; dimensions of the author's Maxwell models"]


def test_report_csv_units():
    # real AEDT 2024 R1 export: "Time [us]","Moving1.Torque [NewtonMeter]"
    d = parse_report_csv(FIX / "aedt_report.csv")
    assert np.isclose(d["Time"][1], 27.7777777778e-6)
    assert np.isclose(d["Moving1.Torque"][0], 416.658715682626)
    assert all(len(v) == len(d["Time"]) for v in d.values())


def test_report_csv_unknown_unit(tmp_path):
    f = tmp_path / "x.csv"
    f.write_text('"Time [fortnight]","Moving1.Torque [NewtonMeter]"\n0,1\n')
    with pytest.raises(ValueError):
        parse_report_csv(f)


def test_em_model_key_ignores_variables():
    p = default_params("sleeved_bridged")
    assert em_model_key(p) == em_model_key(p.replace(gamma=100.0, current_rms=200.0, interference=0.3))
    assert em_model_key(p) != em_model_key(p.replace(sleeve_thickness=1.0))
