import math

import numpy as np
from scipy.integrate import solve_bvp

from sipm.analytic import disc_ring, isotropic, orthotropic, summary

A, B = 71.1, 72.6
W = 20000 * 2 * math.pi / 60


def test_isotropic_matches_textbook():
    E, nu, rho, d = 200e3, 0.3, 7850e-12, 0.1
    D = R = isotropic(E, nu, rho)
    s = disc_ring(A, B, D, R, W, d)
    a, b = A, B
    st_ring = (3 + nu) / 8 * rho * W ** 2 * (a * a + 2 * b * b - (1 + 3 * nu) / (3 + nu) * a * a)
    st_disc = (3 + nu) / 8 * rho * W ** 2 * (a * a - (1 + 3 * nu) / (3 + nu) * a * a)
    gap = d - a / E * (st_ring - st_disc)
    p = gap / (a / E * ((b * b + a * a) / (b * b - a * a) + nu) + a * (1 - nu) / E)
    assert math.isclose(s["p_interface"], p, rel_tol=1e-10)
    # solid disc centre: sigma_r = sigma_t
    _, sr0, st0 = s["disc"](np.array([1e-9]))
    assert math.isclose(sr0[0], st0[0], rel_tol=1e-6)


def test_orthotropic_satisfies_ode_and_bcs():
    D = isotropic(200e3, 0.3, 7850e-12)
    R = orthotropic(8.8e3, 125e3, 0.28, 1800e-12)
    s = disc_ring(A, B, D, R, W, 0.1)
    p_int = s["p_interface"]
    Q11, Q22, Q12, rho = R.Q11, R.Q22, R.Q12, R.rho

    def f(r, y):
        u, sr = y
        up = (sr - Q12 * u / r) / Q11
        st = Q12 * up + Q22 * u / r
        return np.vstack([up, -(sr - st) / r - rho * W * W * r])

    sol = solve_bvp(f, lambda ya, yb: np.array([ya[1] + p_int, yb[1]]), np.linspace(A, B, 50), np.zeros((2, 50)),
                    tol=1e-10)
    ud = s["disc"](np.array([A]))[0][0]
    assert math.isclose(sol.sol(A)[0] - ud, 0.1, rel_tol=1e-7)
    u, sr, st = s["ring"](np.array([B]))
    assert abs(sr[0]) < 1e-9


def test_separation_detected():
    D = isotropic(200e3, 0.3, 7850e-12)
    R = isotropic(70e3, 0.33, 7850e-12)
    s = disc_ring(A, B, D, R, 4 * W, 0.0)  # no interference, fast: ring lifts off
    assert not s["in_contact"] and s["p_interface"] == 0.0


def test_summary_keys():
    D = isotropic(200e3, 0.3, 7850e-12)
    R = orthotropic(8.8e3, 125e3, 0.28, 1800e-12)
    out = summary(A, B, D, R, W, 0.1)
    assert out["in_contact"] and out["sh_ring_max"] > 0 and out["svm_disc_max"] > 0
