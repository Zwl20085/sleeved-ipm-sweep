"""Closed-form reference: rotating solid disc with a shrink-fitted ring (plane stress).

Used to validate the MAPDL path. The disc is isotropic; the ring is cylindrically
orthotropic (isotropic as a special case). With radial displacement u(r) and
plane-stress reduced stiffnesses Q11 (radial), Q22 (hoop), Q12,

    sigma_r = Q11 u' + Q12 u/r,      sigma_t = Q12 u' + Q22 u/r,

equilibrium d(sigma_r)/dr + (sigma_r - sigma_t)/r + rho w^2 r = 0 becomes

    Q11 (u'' + u'/r) - Q22 u/r^2 = -rho w^2 r,

whose solution is u = A r^k + B r^-k + C r^3 with k = sqrt(Q22/Q11) and
C = -rho w^2 / (9 Q11 - Q22). The solid disc has B = 0. The constants follow from
a traction-free ring OD, continuity of sigma_r at the interface and the
interference condition u_ring(a) - u_disc(a) = delta. If the computed contact
pressure is negative the parts separate and are solved independently.

Units: mm, MPa, t/mm^3 (= kg/m^3 * 1e-12), rad/s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class Body:
    Q11: float
    Q22: float
    Q12: float
    rho: float  # t/mm^3

    @property
    def k(self):
        return math.sqrt(self.Q22 / self.Q11)

    def C(self, w):
        return -self.rho * w * w / (9 * self.Q11 - self.Q22)

    def basis(self, r, w):
        """Rows: u, sigma_r, sigma_t for (A-term, B-term, particular)."""
        k = self.k
        r = np.asarray(r, float)
        uA, uB = r ** k, r ** (-k)
        dA, dB = k * r ** (k - 1), -k * r ** (-k - 1)
        C = self.C(w)
        uP, dP = C * r ** 3, 3 * C * r ** 2
        sr = lambda u, d: self.Q11 * d + self.Q12 * u / r  # noqa: E731
        st = lambda u, d: self.Q12 * d + self.Q22 * u / r  # noqa: E731
        return (np.array([uA, uB, uP]), np.array([sr(uA, dA), sr(uB, dB), sr(uP, dP)]),
                np.array([st(uA, dA), st(uB, dB), st(uP, dP)]))


def isotropic(E, nu, rho) -> Body:
    d = 1 - nu * nu
    return Body(E / d, E / d, nu * E / d, rho)


def orthotropic(E_r, E_t, nu_tr, rho) -> Body:
    """nu_tr = -eps_r/eps_t under hoop stress (major ratio for a hoop-wound ring)."""
    nu_rt = nu_tr * E_r / E_t
    d = 1 - nu_rt * nu_tr
    return Body(E_r / d, E_t / d, nu_rt * E_t / d, rho)


def disc_ring(a, b, disc: Body, ring: Body, w, delta):
    """Solve the disc (0..a) + ring (a..b) problem.

    Returns a dict with callables ``disc(r)`` and ``ring(r)`` giving (u, sr, st), the
    interface pressure and a flag telling whether the parts stay in contact.
    """
    def solve(couple: bool):
        # unknowns x = [A_d, A_r, B_r]
        _, srd, _ = disc.basis(a, w)
        ud, _, _ = disc.basis(a, w)
        ura, sra, _ = ring.basis(a, w)
        _, srb, _ = ring.basis(b, w)
        M = np.zeros((3, 3))
        rhs = np.zeros(3)
        # 1) ring OD traction free
        M[0, 1], M[0, 2] = srb[0], srb[1]
        rhs[0] = -srb[2]
        if couple:
            # 2) sigma_r continuity at a
            M[1, 0] = -srd[0]
            M[1, 1], M[1, 2] = sra[0], sra[1]
            rhs[1] = srd[2] - sra[2]
            # 3) u_ring(a) - u_disc(a) = delta
            M[2, 0] = -ud[0]
            M[2, 1], M[2, 2] = ura[0], ura[1]
            rhs[2] = delta - ura[2] + ud[2]
        else:
            # free disc OD and free ring bore
            M[1, 0] = srd[0]
            rhs[1] = -srd[2]
            M[2, 1], M[2, 2] = sra[0], sra[1]
            rhs[2] = -sra[2]
        return np.linalg.solve(M, rhs)

    x = solve(True)
    _, sra, _ = ring.basis(a, w)
    p_int = -(sra[0] * x[1] + sra[1] * x[2] + sra[2])
    contact = p_int >= 0
    if not contact:
        x = solve(False)
        p_int = 0.0

    def fd(r):
        u, sr, st = disc.basis(r, w)
        return (u[0] * x[0] + u[2], sr[0] * x[0] + sr[2], st[0] * x[0] + st[2])

    def fr(r):
        u, sr, st = ring.basis(r, w)
        return (u[0] * x[1] + u[1] * x[2] + u[2], sr[0] * x[1] + sr[1] * x[2] + sr[2],
                st[0] * x[1] + st[1] * x[2] + st[2])

    return dict(disc=fd, ring=fr, p_interface=float(p_int), in_contact=bool(contact), coeffs=x)


def von_mises_plane_stress(sr, st):
    return np.sqrt(sr * sr - sr * st + st * st)


def summary(a, b, disc: Body, ring: Body, w, delta, n=2001) -> dict:
    """Peak quantities comparable with the FE results file."""
    sol = disc_ring(a, b, disc, ring, w, delta)
    rd = np.linspace(1e-6 * a, a, n)
    rr = np.linspace(a, b, n)
    _, srd, std = sol["disc"](rd)
    _, srr, str_ = sol["ring"](rr)
    vm = von_mises_plane_stress(srd, std)
    return dict(svm_disc_max=float(vm.max()), sh_ring_max=float(str_.max()), sh_ring_at_bore=float(str_[0]),
                sr_ring_min=float(srr.min()), p_interface=sol["p_interface"], in_contact=sol["in_contact"])
