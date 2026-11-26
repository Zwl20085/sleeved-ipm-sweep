"""Edge cases of the plot functions: degenerate tables, caller-supplied axes and labels,
the note fallback, backend selection and figure lifetime."""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from sipm import _style as st
from sipm import post
from sipm.cli import plot_results
from sipm.params import default_params


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


def _em(gammas=(100, 115, 130, 145), thicknesses=(1.0, 2.0)) -> pd.DataFrame:
    return pd.DataFrame([dict(status="ok", topology="sleeved_bridged", gamma=g, sleeve_thickness=t,
                              torque_avg=400 - (g - 125) ** 2 / 40 - 20 * t) for g in gammas for t in thicknesses])


def _mech(thicknesses=(0.0, 1.0, 2.0), interferences=(0.0, 0.1, 0.2)) -> pd.DataFrame:
    return pd.DataFrame([dict(status="ok", sleeve_thickness=t, interference=d, svm_rotor_max=600 - 150 * d * t - 20 * t,
                              sh_sleeve_max=50 + 1500 * d) for t in thicknesses for d in interferences])


def test_pyplot_does_not_probe_gui_backends():
    """Resolving the default backend imports Qt/Tk to probe them; on a machine with a broken
    Qt install that prints 'Windows fatal exception: code 0xc0000139' under faulthandler."""
    code = ("import faulthandler, sys; faulthandler.enable()\n"
            "from sipm import _style\n"
            "plt = _style.pyplot()\n"
            "import matplotlib\n"
            "probed = sorted(m for m in sys.modules if m.startswith(('PyQt', 'PySide', 'tkinter', "
            "'matplotlib.backends.backend_qt', 'matplotlib.backends.backend_tk')))\n"
            "print(matplotlib.get_backend().lower(), probed)\n")
    env = {k: v for k, v in os.environ.items() if k not in ("MPLBACKEND", "SIPM_INTERACTIVE")}
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, timeout=120)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "agg []"
    assert "fatal exception" not in r.stderr.lower()


def test_mtpa_lines_duplicate_index():
    """Tables joined with pd.concat (no ignore_index) repeat index labels."""
    df = pd.concat([_em(thicknesses=(1.0,)), _em(thicknesses=(1.0,))])
    assert df.index.has_duplicates
    f = post.mtpa_lines(df, "sleeve_thickness")
    assert "130°" in [t.get_text() for t in f.axes[0].texts]


@pytest.mark.parametrize("rows", ["empty", "single", "nan"])
def test_plots_accept_degenerate_tables(rows):
    em, mech = _em(), _mech()
    kac = pd.DataFrame(dict(topology=["sleeved_bridged"] * 2, em_rpm=[1000.0, 2000.0], kac=[1.1, 1.5]))
    pick = {"empty": lambda d: d.iloc[:0], "single": lambda d: d.iloc[:1],
            "nan": lambda d: d.assign(**{c: np.nan for c in ("torque_avg", "svm_rotor_max", "kac") if c in d})}[rows]
    post.mtpa_lines(pick(em), "sleeve_thickness", note="x")
    post.torque_scatter(pick(em), ["gamma", "sleeve_thickness"], reference=400)
    post.stress_lines(pick(mech), note="x", reference=300)
    post.stress_lines(pick(mech), note="x")  # no reference line: the legend may have no entry at all
    post.kac_lines(pick(kac), note="x")


def test_stress_lines_zero_sleeve_reads_no_sleeve():
    f = post.stress_lines(_mech())
    assert [t.get_text() for t in f.axes[0].get_legend().get_texts()] == \
        ["No sleeve", "1.0 mm sleeve", "2.0 mm sleeve"]


def test_place_note_without_free_spot_goes_to_legend_title():
    plt = st.pyplot()
    with st.style():
        fig, ax = plt.subplots(figsize=(st.COLUMN_W_IN, st.DEFAULT_H_IN))
        ax.plot([0, 1], [0, 1], label="A long legend entry that spans most of the axes width")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        leg = st.legend_above(ax, ncol=1)
        xx, yy = np.meshgrid(np.linspace(0, 1, 60), np.linspace(0, 1, 60))
        t = st.place_note(ax, "NOTE TEXT", np.c_[xx.ravel(), yy.ravel()])
        fig.draw_without_rendering()
        assert leg.get_title().get_text() == "NOTE TEXT"
        assert t is leg.get_title()
        assert "NOTE TEXT" not in [x.get_text() for x in ax.texts]


def test_plot_results_closes_figures(tmp_path):
    import matplotlib.pyplot as plt

    plt.close("all")
    files = plot_results(_em(), tmp_path, "em")
    assert len(files) == 2
    assert plt.get_fignums() == []


def test_stress_map_collinear_points_clear_error():
    df = pd.DataFrame(dict(sleeve_thickness=[1.0, 2.0, 3.0], interference=[0.1, 0.2, 0.3],
                           svm_rotor_max=[300.0, 350.0, 400.0]))
    with pytest.raises(ValueError, match="collinear"):
        post.stress_map(df, 350)


def test_unmatched_labels_warn():
    with pytest.warns(UserWarning, match="5.0"):
        f = post.mtpa_lines(_em(), "sleeve_thickness", labels={1.0: "thin", 5.0: "typo"})
    assert [t.get_text() for t in f.axes[0].get_legend().get_texts()][:2] == ["thin", "2.0 mm sleeve"]


def test_geometry_plot_rejects_unknown_callout():
    with pytest.raises(ValueError, match="bogus"):
        post.geometry_plot(default_params("sleeved_bridged"), callouts=["bogus"])


def test_geometry_plot_into_caller_axes():
    import matplotlib.pyplot as plt

    fig, (a, b) = plt.subplots(1, 2)
    n = len(plt.get_fignums())
    got = post.geometry_plot(default_params("sleeved_bridgeless"), ax=b, heading="mine")
    assert got is fig and len(plt.get_fignums()) == n
    assert b.get_title() == "mine" and a.get_title() == ""
    assert fig.legends == []  # no legend by default when the caller owns the figure


def test_inputs_not_mutated():
    em, labels = _em(), {1.0: "thin"}
    em_before, labels_before = em.copy(), dict(labels)
    post.mtpa_lines(em, "sleeve_thickness", labels=labels, note="x")
    post.torque_scatter(em, ["gamma"])
    pd.testing.assert_frame_equal(em, em_before)
    assert labels == labels_before
    mech = _mech()
    mech_before = mech.copy()
    post.stress_map(mech, 500, note="x")
    post.feasible_points(mech, 500)
    pd.testing.assert_frame_equal(mech, mech_before)


def _callout_texts(fig) -> list[str]:
    from sipm._geomfig import CALLOUT_TEXT

    known = set(CALLOUT_TEXT.values())
    return [t.get_text() for ax in fig.axes for t in ax.texts if t.get_text() in known]


def test_topology_figure_custom_headings():
    f = post.topology_figure(headings=["h = 1", "h = 2"], note=False)
    assert [ax.get_title() for ax in f.axes] == ["(a) h = 1", "(b) h = 2"]


@pytest.mark.parametrize("headings", [["only one"], ["a", "b", "c"], []])
def test_topology_figure_headings_wrong_length(headings):
    with pytest.raises(ValueError, match="headings"):
        post.topology_figure(headings=headings)


def test_topology_figure_callouts_off():
    assert _callout_texts(post.topology_figure(note=False))  # default draws some callouts
    assert _callout_texts(post.topology_figure(callouts=False, note=False)) == []
