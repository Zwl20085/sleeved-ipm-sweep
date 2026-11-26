import math

import numpy as np
import pytest

from sipm import composite as c
from sipm.composite import (
    FIBRES,
    MATRICES,
    MODELS,
    SLEEVE_KEYS,
    STRENGTHS,
    Fibre,
    Lamina,
    Matrix,
    Strength,
    angle_ply,
    failure_indices,
    lamina,
    sleeve_params,
)

CF, EP = FIBRES["carbon_sm"], MATRICES["epoxy"]
ANGLES = [0.0, 7.5, 15.0, 30.0, 45.0, 54.7, 60.0, 75.0, 90.0]


def _iso_fibre(m: Matrix) -> Fibre:
    return Fibre(E1=m.E, E2=m.E, G12=m.G, nu12=m.nu, nu23=m.nu, density=m.density, cte1=m.cte, cte2=m.cte)


def _is_spd(S):
    return np.abs(S - S.T).max() <= 1e-12 * np.abs(S).max() and np.linalg.eigvalsh(S).min() > 0


# ---------------------------------------------------------------------------------------------- limits
@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("vf", [0.0, 0.3, 0.6, 1.0])
def test_fibre_equal_to_matrix_returns_matrix(model, vf):
    lam = lamina(_iso_fibre(EP), EP, vf, model)
    for n in ("E1", "E2", "E3"):
        assert math.isclose(getattr(lam, n), EP.E, rel_tol=1e-12)
    for n in ("G12", "G13", "G23"):
        assert math.isclose(getattr(lam, n), EP.G, rel_tol=1e-12)
    for n in ("nu12", "nu13", "nu23"):
        assert math.isclose(getattr(lam, n), EP.nu, rel_tol=1e-12)
    for n in ("cte1", "cte2", "cte3"):
        assert math.isclose(getattr(lam, n), EP.cte, rel_tol=1e-12)
    assert math.isclose(lam.density, EP.density)


@pytest.mark.parametrize("model", MODELS)
def test_vf_zero_returns_matrix(model):
    lam = lamina(CF, EP, 0.0, model)
    assert np.allclose([lam.E1, lam.E2, lam.E3], EP.E, rtol=1e-12)
    assert np.allclose([lam.G12, lam.G13, lam.G23], EP.G, rtol=1e-12)
    assert np.allclose([lam.nu12, lam.nu13, lam.nu23], EP.nu, rtol=1e-12)
    assert np.allclose([lam.cte1, lam.cte2, lam.cte3], EP.cte, rtol=1e-12)
    assert lam.density == EP.density
    # and continuously so
    lam = lamina(CF, EP, 1e-9, model)
    assert math.isclose(lam.E2, EP.E, rel_tol=1e-3)


@pytest.mark.parametrize("model", MODELS)
def test_vf_one_returns_fibre(model):
    lam = lamina(CF, EP, 1.0, model)
    assert np.allclose([lam.E1, lam.E2, lam.G12, lam.G23, lam.nu12, lam.nu23],
                       [CF.E1, CF.E2, CF.G12, CF.G23, CF.nu12, CF.nu23], rtol=1e-12)
    assert np.allclose([lam.cte1, lam.cte2], [CF.cte1, CF.cte2], rtol=1e-12)


@pytest.mark.parametrize("model", MODELS)
def test_rule_of_mixtures_E1_nu12_density_exact(model):
    for vf in np.linspace(0, 1, 11):
        lam = lamina(CF, EP, vf, model)
        assert lam.E1 == pytest.approx(vf * CF.E1 + (1 - vf) * EP.E, rel=1e-14)
        assert lam.nu12 == pytest.approx(vf * CF.nu12 + (1 - vf) * EP.nu, rel=1e-14)
        assert lam.density == pytest.approx(vf * CF.density + (1 - vf) * EP.density, rel=1e-14)


def test_rom_transverse_is_reuss_and_bounds():
    vf = 0.55
    rom, ch, ht = (lamina(CF, EP, vf, m) for m in ("rom", "chamis", "halpin_tsai"))
    assert rom.E2 == pytest.approx(1 / (vf / CF.E2 + (1 - vf) / EP.E), rel=1e-14)
    # Reuss is the lower bound of the transverse estimates; Voigt the upper
    voigt = vf * CF.E2 + (1 - vf) * EP.E
    for lam in (ch, ht):
        assert rom.E2 < lam.E2 < voigt


def test_schapery_cte_independent_derivation():
    """Phase-wise 1D check: common axial strain, zero transverse stress, sum of transverse strains."""
    vf, f, m = 0.6, CF, EP
    lam = lamina(f, m, vf)
    # axial: force balance with common strain a1
    sf = f.E1 * (lam.cte1 - f.cte1)
    sm = m.E * (lam.cte1 - m.cte)
    assert vf * sf + (1 - vf) * sm == pytest.approx(0.0, abs=1e-12)
    eps_f2 = f.cte2 - f.nu12 * sf / f.E1
    eps_m2 = m.cte - m.nu * sm / m.E
    assert lam.cte2 == pytest.approx(vf * eps_f2 + (1 - vf) * eps_m2, rel=1e-12)


def test_halpin_tsai_published_value():
    """Halpin-Tsai E2 of glass/epoxy, Ef = 85 GPa, nuf = 0.2, Em = 3.4 GPa, num = 0.3, vf = 0.7 (xi = 2):
    20.20 GPa, quoted as a textbook value (attributed to Kaw, Mechanics of Composite Materials, 2nd ed.,
    2006; not checked against the book). The assertion is the Halpin-Tsai formula evaluated by hand."""
    g = Fibre(E1=85.0, E2=85.0, G12=85.0 / 2.4, nu12=0.2, nu23=0.2, density=2500.0, cte1=5e-6, cte2=5e-6)
    m = Matrix(E=3.4, nu=0.3, density=1200.0, cte=60e-6)
    lam = lamina(g, m, 0.7, "halpin_tsai")
    assert round(lam.E2, 2) == 20.20
    assert round(lam.E1, 2) == 60.52
    assert round(lamina(g, m, 0.7, "rom").E2, 2) == 10.37


# ---------------------------------------------------------------------------------------------- SPD
@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("fname", list(FIBRES))
@pytest.mark.parametrize("mname", list(MATRICES))
def test_all_compliances_spd(model, fname, mname):
    for vf in np.linspace(0, 1, 21):
        lam = lamina(FIBRES[fname], MATRICES[mname], vf, model)
        assert _is_spd(lam.compliance())
        for th in ANGLES:
            L = angle_ply(lam, th)
            assert _is_spd(L.S)
            assert _is_spd(L.C)


# ---------------------------------------------------------------------------------------------- rotation
def test_bond_matrix_properties():
    rng = np.random.default_rng(1)
    for _ in range(5):
        Q, _r = np.linalg.qr(rng.normal(size=(3, 3)))
        M = c._bond(Q)
        sig = rng.normal(size=(3, 3))
        sig = sig + sig.T
        v = np.array([sig[i, j] for i, j in c._VOIGT])
        sp = Q @ sig @ Q.T
        assert np.allclose(M @ v, [sp[i, j] for i, j in c._VOIGT])
        assert np.allclose(c._bond(Q.T), np.linalg.inv(M))


def _qbar_closed_form(E1, E2, nu12, G12, th_deg):
    """Transformed reduced stiffness Qbar of an off-axis lamina (the standard formulas, e.g. R.M. Jones,
    Mechanics of Composite Materials, 2nd ed., 1999, ch. 2)."""
    nu21 = nu12 * E2 / E1
    d = 1 - nu12 * nu21
    Q11, Q22, Q12, Q66 = E1 / d, E2 / d, nu12 * E2 / d, G12
    t = math.radians(th_deg)
    m, n = math.cos(t), math.sin(t)
    q11 = Q11 * m ** 4 + 2 * (Q12 + 2 * Q66) * n * n * m * m + Q22 * n ** 4
    q22 = Q11 * n ** 4 + 2 * (Q12 + 2 * Q66) * n * n * m * m + Q22 * m ** 4
    q12 = (Q11 + Q22 - 4 * Q66) * n * n * m * m + Q12 * (n ** 4 + m ** 4)
    q66 = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * n * n * m * m + Q66 * (n ** 4 + m ** 4)
    q16 = (Q11 - Q12 - 2 * Q66) * n * m ** 3 + (Q12 - Q22 + 2 * Q66) * n ** 3 * m
    q26 = (Q11 - Q12 - 2 * Q66) * n ** 3 * m + (Q12 - Q22 + 2 * Q66) * n * m ** 3
    return np.array([[q11, q12, q16], [q12, q22, q26], [q16, q26, q66]])


def _rotated_qbar(lam, th):
    """Plane-stress reduced stiffness (hoop, axial, hoop-axial) of the 3D-rotated ply."""
    M = c._bond(c._ply_to_cyl(th))
    Minv = np.linalg.inv(M)
    S = Minv.T @ lam.compliance() @ Minv
    return np.linalg.inv(S[np.ix_([0, 1, 5], [0, 1, 5])])


def test_rotation_matches_closed_form_qbar():
    lam = lamina(CF, EP, 0.6)
    for th in np.linspace(-90, 90, 37):
        assert np.allclose(_rotated_qbar(lam, th), _qbar_closed_form(lam.E1, lam.E2, lam.nu12, lam.G12, th),
                           rtol=1e-12, atol=1e-10)


def test_published_off_axis_example():
    """Graphite/epoxy E1 = 181, E2 = 10.3, G12 = 7.17 GPa, nu12 = 0.28 at 60 deg:
    Qbar = [[23.65, 32.46, 20.05], [32.46, 109.4, 54.19], [20.05, 54.19, 36.74]] GPa, the worked example
    of an off-axis lamina attributed to Kaw, Mechanics of Composite Materials, 2nd ed. (2006) (recalled
    from memory, not checked against the book). The values are reproduced independently by the
    closed-form reduced stiffness of :func:`_qbar_closed_form` (Jones' transformation formulas)."""
    lam = Lamina(E1=181.0, E2=10.3, E3=10.3, nu12=0.28, nu13=0.28, nu23=0.4, G12=7.17, G13=7.17,
                 G23=10.3 / 2.8, density=1600.0, cte1=0.0, cte2=0.0, cte3=0.0)
    ref = np.array([[23.65, 32.46, 20.05], [32.46, 109.4, 54.19], [20.05, 54.19, 36.74]])
    assert np.allclose(_rotated_qbar(lam, 60.0), ref, rtol=6e-4, atol=0)
    assert np.allclose(_qbar_closed_form(181.0, 10.3, 0.28, 7.17, 60.0), ref, rtol=6e-4, atol=0)


# ---------------------------------------------------------------------------------------------- angle ply
def test_angle_ply_zero_reproduces_lamina():
    lam = lamina(CF, EP, 0.6)
    L = angle_ply(lam, 0.0)
    assert np.allclose(L.C, lam.stiffness(), rtol=1e-12, atol=1e-12)
    assert np.allclose(L.alpha, lam.alpha, rtol=1e-12, atol=1e-20)
    assert L.E_hoop == pytest.approx(lam.E1) and L.E_axial == pytest.approx(lam.E2)
    assert L.E_radial == pytest.approx(lam.E3)
    assert L.nu_hoop_radial == pytest.approx(lam.nu13) and L.nu_hoop_axial == pytest.approx(lam.nu12)
    assert L.G_hoop_axial == pytest.approx(lam.G12) and L.G_radial_axial == pytest.approx(lam.G23)


def test_angle_ply_ninety_swaps_hoop_and_axial():
    lam = lamina(CF, EP, 0.6, "halpin_tsai")
    C0, C90 = angle_ply(lam, 0.0).C, angle_ply(lam, 90.0).C
    p = [1, 0, 2, 4, 3, 5]   # swap hoop<->axial and the two transverse shears
    assert np.allclose(C90, C0[np.ix_(p, p)], rtol=1e-12, atol=1e-9)
    L0, L90 = angle_ply(lam, 0.0), angle_ply(lam, 90.0)
    assert L90.E_hoop == pytest.approx(L0.E_axial) and L90.E_axial == pytest.approx(L0.E_hoop)
    assert L90.cte_hoop == pytest.approx(L0.cte_axial) and L90.cte_axial == pytest.approx(L0.cte_hoop)
    assert L90.G_hoop_radial == pytest.approx(L0.G_radial_axial)


def test_angle_ply_even_in_theta_and_orthotropic():
    lam = lamina(FIBRES["carbon_hm"], EP, 0.55)
    for th in ANGLES:
        Lp, Lm, Ls = angle_ply(lam, th), angle_ply(lam, -th), angle_ply(lam, 180.0 - th)
        for L in (Lm, Ls):
            assert np.allclose(L.C, Lp.C, rtol=1e-12, atol=1e-9)
            assert np.allclose(L.alpha, Lp.alpha, rtol=1e-10, atol=1e-18)
        # no monoclinic coupling: normal <-> hoop-axial shear and the two transverse shears
        scale = np.abs(Lp.C).max()
        for i, j in [(0, 5), (1, 5), (2, 5), (3, 4), (0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5),
                     (4, 5)]:
            assert abs(Lp.C[i, j]) < 1e-12 * scale
        assert np.allclose(Lp.alpha[3:], 0.0, atol=1e-20)
    L45 = angle_ply(lam, 45.0)
    assert L45.E_hoop == pytest.approx(L45.E_axial, rel=1e-12)


def _direct_homogenisation(Cp, Cm, eps):
    """Independent route: given the macroscopic strain, find ply out-of-plane strains such that the
    out-of-plane stresses are continuous and their average is the macroscopic one; return mean stress."""
    I, O = [0, 1, 5], [2, 3, 4]
    # unknowns z = [eo_plus (3), eo_minus (3)]
    A = np.zeros((6, 6))
    b = np.zeros(6)
    A[:3, :3], A[:3, 3:] = Cp[np.ix_(O, O)], -Cm[np.ix_(O, O)]
    b[:3] = -(Cp[np.ix_(O, I)] - Cm[np.ix_(O, I)]) @ eps[I]
    A[3:, :3], A[3:, 3:] = 0.5 * np.eye(3), 0.5 * np.eye(3)
    b[3:] = eps[O]
    z = np.linalg.solve(A, b)
    sig = []
    for C, eo in ((Cp, z[:3]), (Cm, z[3:])):
        e = np.array(eps, dtype=float)
        e[O] = eo
        sig.append(C @ e)
    return 0.5 * (sig[0] + sig[1])


def test_homogenisation_matches_direct_layer_solution():
    lam = lamina(CF, EP, 0.6)
    rng = np.random.default_rng(7)
    for th in (15.0, 30.0, 52.0, 80.0):
        L = angle_ply(lam, th)
        Cs = []
        for s in (+1, -1):
            M = c._bond(c._ply_to_cyl(s * th))
            Cs.append(M @ lam.stiffness() @ M.T)
        for _ in range(4):
            eps = rng.normal(size=6) * 1e-3
            assert np.allclose(L.C @ eps, _direct_homogenisation(*Cs, eps), rtol=1e-10, atol=1e-14)
        # and it is not the plain Voigt average
        voigt = 0.5 * (Cs[0] + Cs[1])
        assert np.abs(voigt - L.C).max() > 1e-6 * np.abs(L.C).max()


# ---------------------------------------------------------------------------------------------- ply stresses
@pytest.mark.parametrize("dT", [0.0, -120.0])
def test_ply_stresses_average_back_to_input(dT):
    rng = np.random.default_rng(3)
    lam = lamina(FIBRES["carbon_im"], EP, 0.58, "chamis")
    for th in rng.uniform(-90, 90, 8):
        L = angle_ply(lam, th)
        sig = rng.normal(scale=200.0, size=(50, 6))
        sp, sm = L.ply_stresses(sig, dT)
        assert sp.shape == sm.shape == (50, 6)
        back_p = sp @ c._bond(c._ply_to_cyl(th)).T
        back_m = sm @ c._bond(c._ply_to_cyl(-th)).T
        assert np.allclose(0.5 * (back_p + back_m), sig, rtol=1e-10, atol=1e-8)
        # out-of-plane stresses are common to both plies
        assert np.allclose(back_p[:, [2, 3, 4]], sig[:, [2, 3, 4]], atol=1e-8)
        assert np.allclose(back_m[:, [2, 3, 4]], sig[:, [2, 3, 4]], atol=1e-8)
        # in-plane strains common: compute each ply's strain from its own law
        Cm = lam.stiffness() * 1e3
        for s_ply, sgn in ((sp, 1), (sm, -1)):
            eps_mat = (s_ply @ np.linalg.inv(Cm).T) + lam.alpha * dT
            N = np.linalg.inv(c._bond(c._ply_to_cyl(sgn * th))).T
            eps_cyl = eps_mat @ N.T
            eps_eff = sig @ (L.S / 1e3).T + L.alpha * dT
            assert np.allclose(eps_cyl[:, [0, 1, 5]], eps_eff[:, [0, 1, 5]], rtol=1e-9, atol=1e-13)


def test_ply_stresses_single_state_and_simple_cases():
    lam = lamina(CF, EP, 0.6)
    L0 = angle_ply(lam, 0.0)
    sig = np.array([500.0, 20.0, -10.0, 1.0, 2.0, 3.0])
    sp, sm = L0.ply_stresses(sig)
    assert sp.shape == (6,)
    assert np.allclose(sp, sig) and np.allclose(sm, sig)
    # pure hoop winding, free thermal expansion: no residual stress
    sp, sm = L0.ply_stresses(np.zeros(6), delta_T=-150.0)
    assert np.allclose(sp, 0, atol=1e-9) and np.allclose(sm, 0, atol=1e-9)
    # +/-30 winding cooled from cure: self-equilibrated residual stresses, fibre tension/compression mirror
    L = angle_ply(lam, 30.0)
    sp, sm = L.ply_stresses(np.zeros(6), delta_T=-150.0)
    assert abs(sp[0]) > 1.0
    assert np.allclose(sp[[0, 1, 2]], sm[[0, 1, 2]]) and np.allclose(sp[5], -sm[5])
    # vectorised delta_T
    sp, _ = L.ply_stresses(np.zeros((3, 6)), delta_T=np.array([0.0, -50.0, -100.0]))
    assert np.allclose(sp[2], 2 * sp[1]) and np.allclose(sp[0], 0)


# ---------------------------------------------------------------------------------------------- failure
@pytest.mark.parametrize("sname", list(STRENGTHS))
def test_failure_indices_unity_at_uniaxial_strengths(sname):
    X = STRENGTHS[sname]
    cases = {  # (s1, s2, t12): criteria that must give exactly 1
        "Xt": ((X.Xt, 0, 0), ("max_stress", "tsai_wu", "hashin_fibre")),
        "Xc": ((-X.Xc, 0, 0), ("max_stress", "tsai_wu", "hashin_fibre")),
        "Yt": ((0, X.Yt, 0), ("max_stress", "tsai_wu", "hashin_matrix")),
        "Yc": ((0, -X.Yc, 0), ("max_stress", "tsai_wu", "hashin_matrix")),
        "+S": ((0, 0, X.S), ("max_stress", "tsai_wu", "hashin_matrix", "hashin_fibre")),
        "-S": ((0, 0, -X.S), ("max_stress", "tsai_wu", "hashin_matrix", "hashin_fibre")),
    }
    for (s1, s2, t), names in cases.values():
        r = failure_indices(s1, s2, t, X)
        for n in names:
            assert r[n] == pytest.approx(1.0, rel=1e-12), n
            assert r[n + "_sr"] == pytest.approx(1.0, rel=1e-12), n
    # and Hashin with a separate transverse shear strength
    Xs = Strength(X.Xt, X.Xc, X.Yt, X.Yc, X.S, St=0.4 * X.Yc)
    assert failure_indices(0, -X.Yc, 0, Xs)["hashin_matrix"] == pytest.approx(1.0, rel=1e-12)


def test_strength_ratio_scales_and_lands_on_envelope():
    rng = np.random.default_rng(11)
    X = STRENGTHS["carbon_sm_epoxy"]
    s = rng.normal(size=(3, 400)) * np.array([[800.0], [60.0], [40.0]])
    r = failure_indices(*s, X)
    for n in ("max_stress", "tsai_wu", "hashin_fibre", "hashin_matrix"):
        R = r[n + "_sr"]
        on = failure_indices(*(s * R), X)[n]
        assert np.allclose(on, 1.0, rtol=1e-10), n
        assert np.all((R > 1) == (r[n] < 1)), n   # safety factor > 1 exactly when the index is < 1
    k = 1.7
    rk = failure_indices(*(k * s), X)
    for n in ("max_stress", "tsai_wu", "hashin_fibre", "hashin_matrix"):
        assert np.allclose(rk[n + "_sr"], r[n + "_sr"] / k, rtol=1e-12)
    assert np.allclose(r["max_stress_sr"], 1 / r["max_stress"])


def test_failure_zero_stress_and_tsai_wu_f12():
    X = STRENGTHS["carbon_sm_epoxy"]
    r = failure_indices(0.0, 0.0, 0.0, X)
    for n in ("max_stress", "tsai_wu", "hashin_fibre", "hashin_matrix"):
        assert r[n] == 0 and np.isinf(r[n + "_sr"])
    # equal biaxial s1 = s2 = s: F12 = -0.5 sqrt(F11 F22) enters as 2 F12 s^2
    s = 30.0
    F1, F2 = 1 / X.Xt - 1 / X.Xc, 1 / X.Yt - 1 / X.Yc
    F11, F22 = 1 / (X.Xt * X.Xc), 1 / (X.Yt * X.Yc)
    exp = (F1 + F2) * s + (F11 + F22 - math.sqrt(F11 * F22)) * s * s
    assert failure_indices(s, s, 0.0, X)["tsai_wu"] == pytest.approx(exp, rel=1e-12)


# ---------------------------------------------------------------------------------------------- params link
def test_sleeve_params_keys_and_convention():
    lam = lamina(CF, EP, 0.6)
    d = sleeve_params(lam)
    assert tuple(d) == SLEEVE_KEYS
    assert d == sleeve_params(angle_ply(lam, 0.0))
    assert all(isinstance(v, float) for v in d.values())
    # rebuild the compliance the way Params/MAPDL read the constants and compare with the laminate
    L = angle_ply(lam, 35.0)
    p = sleeve_params(L)
    Eh, Ea, Er = p["sleeve_E_hoop"], p["sleeve_E_axial"], p["sleeve_E_radial"]
    S = np.zeros((6, 6))
    S[0, 0], S[1, 1], S[2, 2] = 1 / Eh, 1 / Ea, 1 / Er
    S[0, 1] = S[1, 0] = -p["sleeve_nu_hoop_axial"] / Eh
    S[0, 2] = S[2, 0] = -p["sleeve_nu_hoop_radial"] / Eh
    S[1, 2] = S[2, 1] = -p["sleeve_nu_radial_axial"] / Er
    S[3, 3], S[4, 4], S[5, 5] = 1 / p["sleeve_G_radial_axial"], 1 / p["sleeve_G_hoop_radial"], \
        1 / p["sleeve_G_hoop_axial"]
    assert np.allclose(S, L.S, rtol=1e-10, atol=1e-14)
    assert [p["sleeve_cte_hoop"], p["sleeve_cte_axial"], p["sleeve_cte_radial"]] == pytest.approx(L.alpha[:3])


def test_standard_modulus_preset_against_table_ii():
    """Documented link to the Params defaults (paper Table II: E_hoop 125, E_radial 8.8 GPa, nu 0.28,
    rho 1800). Reported, not tuned: the preset reaches E1 = 125 GPa at vf ~ 0.536."""
    vf = (125.0 - EP.E) / (CF.E1 - EP.E)
    assert 0.5 < vf < 0.6
    lam = lamina(CF, EP, vf, "chamis")
    assert lam.E1 == pytest.approx(125.0, rel=1e-12)
    assert lam.E2 == pytest.approx(7.98, abs=0.01)
    assert lam.nu12 == pytest.approx(0.270, abs=0.001)
    assert lam.density == pytest.approx(1522.0, abs=1.0)


# ---------------------------------------------------------------------------------------------- input checks
def test_input_validation():
    with pytest.raises(ValueError):
        lamina(CF, EP, 1.2)
    with pytest.raises(ValueError):
        lamina(CF, EP, 0.5, "mori_tanaka")
    with pytest.raises(ValueError):
        Fibre(E1=230, E2=15, G12=15, nu12=0.2, density=1800, cte1=0, cte2=0)
    with pytest.raises(ValueError):
        Fibre(E1=230, E2=15, G12=15, nu12=0.2, density=1800, cte1=0, cte2=0, G23=5, nu23=0.3)
    assert Fibre(E1=230, E2=15, G12=15, nu12=0.2, density=1800, cte1=0, cte2=0, G23=6.0).nu23 == \
        pytest.approx(0.25)
    with pytest.raises(ValueError):
        Strength(Xt=1, Xc=-1, Yt=1, Yc=1, S=1)
    with pytest.raises(TypeError):
        sleeve_params(CF)


# ---------------------------------------------------------------------------------------------- ply/sleeve failure
def test_ply_failure_both_transverse_planes():
    st = STRENGTHS["carbon_sm_epoxy"]
    z = np.zeros(4)
    s = np.zeros((4, 6))
    s[0, 1] = st.Yt          # transverse tension in the 12 plane
    s[1, 2] = st.Yt          # the same through the thickness (33), seen only by the 13 plane
    s[2, 4] = st.S           # t13 shear, seen only by the 13 plane
    s[3, 3] = 1e3            # t23 is ignored by construction
    fi = c.ply_failure([s], st)
    assert np.allclose(fi["max_stress"], [1, 1, 1, 0])
    assert np.allclose(fi["hashin_matrix"], [1, 1, 1, 0])
    assert np.allclose(fi["max_stress_sr"][:3], 1) and np.isinf(fi["max_stress_sr"][3])
    one = c.failure_indices(z[0], st.Yt, z[0], st)     # uniaxial transverse strength: on the envelope
    assert np.allclose(fi["tsai_wu"][[0, 1]], one["tsai_wu"]) and np.isclose(one["tsai_wu"], 1.0)
    assert set(fi) == {k for n in c.CRITERIA for k in (n, f"{n}_sr")}


def test_ply_failure_worst_over_plies_and_shape():
    st = STRENGTHS["carbon_sm_epoxy"]
    rng = np.random.default_rng(3)
    a, b = rng.normal(scale=40.0, size=(5, 6)), rng.normal(scale=40.0, size=(5, 6))
    fab, fa, fb = c.ply_failure([a, b], st), c.ply_failure([a], st), c.ply_failure([b], st)
    for n in c.CRITERIA:
        assert np.array_equal(fab[n], np.maximum(fa[n], fb[n]))
        assert np.array_equal(fab[f"{n}_sr"], np.minimum(fa[f"{n}_sr"], fb[f"{n}_sr"]))
        assert fab[n].shape == (5,)
    with pytest.raises(ValueError):
        c.ply_failure([np.zeros((3, 5))], st)
    with pytest.raises(ValueError):
        c.ply_failure([], st)


def test_sleeve_failure_hoop_ply_and_laminate_zero_agree():
    st = STRENGTHS["carbon_sm_epoxy"]
    sig = np.array([[300.0, 0.0, -5.0, 0.0, 2.0, 0.0], [-100.0, 0.0, 1.0, 0.0, -3.0, 0.0]])
    lam0 = angle_ply(lamina(CF, EP, 0.6), 0.0)
    a, b = c.sleeve_failure(sig, st), c.sleeve_failure(sig, st, lam0, delta_T=80.0)
    for k in a:
        assert np.allclose(a[k], b[k], rtol=1e-9, atol=1e-12), k
    # hoop ply: fibre index from the hoop stress, compression branch for the second node
    assert a["max_stress"][0] == pytest.approx(max(300 / st.Xt, 5 / st.Yc, 2 / st.S))
    assert a["hashin_fibre"][1] == pytest.approx((100 / st.Xc) ** 2)


def test_sleeve_failure_angle_ply_uses_ply_stresses():
    st = STRENGTHS["carbon_sm_epoxy"]
    L = angle_ply(lamina(CF, EP, 0.6), 30.0)
    sig = np.array([[200.0, 0.0, -4.0, 0.0, 1.0, 0.0]])
    plus, minus = L.ply_stresses(sig, 50.0)
    ref = c.ply_failure([plus, minus], st)
    out = c.sleeve_failure(sig, st, L, delta_T=50.0)
    for k in out:
        assert np.array_equal(out[k], ref[k])
    # the hoop-axial shear in the 30 deg plies is not zero, unlike the homogenised stress
    assert abs(plus[0, 5]) > 1.0


@pytest.mark.parametrize("f12", [-1.0, 1.0, -1.5, 2.0])
def test_tsai_wu_interaction_outside_closed_envelope_rejected(f12):
    with pytest.raises(ValueError, match="f12_star"):
        failure_indices(100.0, -150.0, 10.0, STRENGTHS["carbon_sm_epoxy"], f12_star=f12)
    # just inside, the strength ratio is finite for a biaxial state the open quadric would miss
    fi = failure_indices(100.0, -150.0, 10.0, STRENGTHS["carbon_sm_epoxy"], f12_star=0.999 * np.sign(f12))
    assert np.isfinite(fi["tsai_wu_sr"]) and fi["tsai_wu_sr"] > 0
