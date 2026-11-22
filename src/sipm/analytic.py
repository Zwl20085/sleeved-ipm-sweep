"""Closed-form reference: rotating solid disc with a shrink-fitted ring (plane stress).

Used to validate the MAPDL path. The disc is isotropic; the ring is cylindrically
orthotropic (isotropic as a special case). Both may carry a uniform temperature change
dT from the reference (fitted) state, with thermal expansion coefficients alpha_r
(radial) and alpha_t (hoop) per body. With radial displacement u(r) and plane-stress
reduced stiffnesses Q11 (radial), Q22 (hoop), Q12,

    sigma_r = Q11 u' + Q12 u/r - beta_r dT,      sigma_t = Q12 u' + Q22 u/r - beta_t dT,
    beta_r = Q11 alpha_r + Q12 alpha_t,          beta_t = Q12 alpha_r + Q22 alpha_t

(the stress-free thermal strains alpha_r dT, alpha_t dT subtracted from the total strains
u', u/r; plane stress: sigma_z = 0, so alpha_z does not enter). Equilibrium
d(sigma_r)/dr + (sigma_r - sigma_t)/r + rho w^2 r = 0 becomes, dT uniform,

    Q11 (u'' + u'/r) - Q22 u/r^2 = -rho w^2 r + (beta_r - beta_t) dT / r,

whose solution is u = A r^k + B r^-k + C r^3 + u_T with k = sqrt(Q22/Q11) and
C = -rho w^2 / (9 Q11 - Q22) (singular for Q22 = 9 Q11, k = 3; not handled).

Thermal particular term. The forcing F/r, F = (beta_r - beta_t) dT, is matched by D r with
Q11 (1 - k^2) D = F, i.e. D = F / (Q11 - Q22). This fails at k = 1, where r is itself a
homogeneous solution, and is ill-conditioned near it. Adding the homogeneous solution
-D r^k (absorbed into A) gives a form valid for every k > 0:

    u_T = c (r^k - r) / (k - 1) = c r L(r),      c = F / (Q11 (1 + k)),
    L(r) = (r^(k-1) - 1) / (k - 1) = ln(r) * expm1(x) / x,   x = (k - 1) ln(r),
    u_T' = c (k L(r) + 1),

(r in mm), which tends to the resonant solution u_T = F/(2 Q11) r ln r at k = 1 (check:
u'' + u'/r - u/r^2 = 2/r for u = r ln r). An isotropic body, or any body with
beta_r = beta_t, has F = 0: its thermal strain is a free uniform expansion u = alpha dT r,
carried by the A r term with zero stress. A solid disc with k = 1 and F != 0 has a
logarithmic stress singularity at the centre (physical, not a defect of the formula); the
validation disc is isotropic.

The solid disc has B = 0. The constants follow from a traction-free ring OD, continuity
of sigma_r at the interface and the interference condition u_ring(a) - u_disc(a) = delta,
where u is measured from the stress-free fitted dimensions at the reference temperature
(so u includes the free thermal expansion and delta is the interference defined at the
reference temperature). If the computed contact pressure is negative the parts separate
and are solved independently.

Units: mm, MPa, t/mm^3 (= kg/m^3 * 1e-12), rad/s, 1/K, K.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def _log_ratio_term(r, k):
    """L(r) = (r^(k-1) - 1)/(k - 1), evaluated stably (= ln r at k = 1)."""
    lr = np.log(r)
    x = (k - 1.0) * lr
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = np.where(np.abs(x) < 1e-8, 1.0 + x / 2.0, np.expm1(x) / np.where(x == 0, 1.0, x))
    return lr * phi


@dataclass
class Body:
    Q11: float
    Q22: float
    Q12: float
    rho: float  # t/mm^3
    alpha_r: float = 0.0  # 1/K, radial thermal expansion
    alpha_t: float = 0.0  # 1/K, hoop thermal expansion

    @property
    def k(self):
        return math.sqrt(self.Q22 / self.Q11)

    @property
    def beta_r(self):
        return self.Q11 * self.alpha_r + self.Q12 * self.alpha_t

    @property
    def beta_t(self):
        return self.Q12 * self.alpha_r + self.Q22 * self.alpha_t

    def C(self, w):
        return -self.rho * w * w / (9 * self.Q11 - self.Q22)

    def c_thermal(self, dT):
        """Coefficient c of the thermal particular term u_T = c r L(r)."""
        return (self.beta_r - self.beta_t) * dT / (self.Q11 * (1 + self.k))

    def basis(self, r, w, dT=0.0):
        """Rows: u, sigma_r, sigma_t for (A-term, B-term, particular).

        The particular column holds the centrifugal (C r^3) and thermal (u_T and the
        -beta dT stress terms) parts for speed ``w`` and uniform temperature change ``dT``.
        """
        k = self.k
        r = np.asarray(r, float)
        uA, uB = r ** k, r ** (-k)
        dA, dB = k * r ** (k - 1), -k * r ** (-k - 1)
        C = self.C(w)
        uP, dP = C * r ** 3, 3 * C * r ** 2
        if dT:
            c = self.c_thermal(dT)
            if c:
                L = _log_ratio_term(r, k)
                uP = uP + c * r * L
                dP = dP + c * (k * L + 1.0)
        sr = lambda u, d: self.Q11 * d + self.Q12 * u / r  # noqa: E731
        st = lambda u, d: self.Q12 * d + self.Q22 * u / r  # noqa: E731
        srP = sr(uP, dP) - self.beta_r * dT
        stP = st(uP, dP) - self.beta_t * dT
        return (np.array([uA, uB, uP]), np.array([sr(uA, dA), sr(uB, dB), srP]),
                np.array([st(uA, dA), st(uB, dB), stP]))


def isotropic(E, nu, rho, alpha=0.0) -> Body:
    d = 1 - nu * nu
    return Body(E / d, E / d, nu * E / d, rho, alpha, alpha)


def orthotropic(E_r, E_t, nu_tr, rho, alpha_r=0.0, alpha_t=0.0) -> Body:
    """nu_tr = -eps_r/eps_t under hoop stress (major ratio for a hoop-wound ring)."""
    nu_rt = nu_tr * E_r / E_t
    d = 1 - nu_rt * nu_tr
    return Body(E_r / d, E_t / d, nu_rt * E_t / d, rho, alpha_r, alpha_t)


def disc_ring(a, b, disc: Body, ring: Body, w, delta, dT=0.0):
    """Solve the disc (0..a) + ring (a..b) problem.

    ``delta`` is the radial interference at the reference temperature, ``dT`` the uniform
    temperature change of both bodies from it. Returns a dict with callables ``disc(r)``
    and ``ring(r)`` giving (u, sr, st), the interface pressure and a flag telling whether
    the parts stay in contact.
    """
    def solve(couple: bool):
        # unknowns x = [A_d, A_r, B_r]
        ud, srd, _ = disc.basis(a, w, dT)
        ura, sra, _ = ring.basis(a, w, dT)
        _, srb, _ = ring.basis(b, w, dT)
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
    _, sra, _ = ring.basis(a, w, dT)
    p_int = -(sra[0] * x[1] + sra[1] * x[2] + sra[2])
    contact = p_int >= 0
    if not contact:
        x = solve(False)
        p_int = 0.0

    def fd(r):
        u, sr, st = disc.basis(r, w, dT)
        return (u[0] * x[0] + u[2], sr[0] * x[0] + sr[2], st[0] * x[0] + st[2])

    def fr(r):
        u, sr, st = ring.basis(r, w, dT)
        return (u[0] * x[1] + u[1] * x[2] + u[2], sr[0] * x[1] + sr[1] * x[2] + sr[2],
                st[0] * x[1] + st[1] * x[2] + st[2])

    return dict(disc=fd, ring=fr, p_interface=float(p_int), in_contact=bool(contact), coeffs=x)


def coupled_pressure(a, b, disc: Body, ring: Body, w, delta, dT=0.0) -> float:
    """Signed interface pressure of the coupled (always in contact) solution, MPa.

    Equal to ``disc_ring(...)["p_interface"]`` while the parts stay in contact; negative where
    :func:`disc_ring` finds them separated (and reports 0). The coupled problem is linear in its
    loads, so this pressure is exactly p(w) = p0(delta, dT) - c w^2."""
    sol = disc_ring(a, b, disc, ring, w, delta, dT)
    if sol["in_contact"]:
        return sol["p_interface"]
    # disc_ring re-solved the separated parts; repeat the coupled solve for the signed value
    x = _coupled_coeffs(a, b, disc, ring, w, delta, dT)
    _, sra, _ = ring.basis(a, w, dT)
    return float(-(sra[0] * x[1] + sra[1] * x[2] + sra[2]))


def _coupled_coeffs(a, b, disc: Body, ring: Body, w, delta, dT):
    ud, srd, _ = disc.basis(a, w, dT)
    ura, sra, _ = ring.basis(a, w, dT)
    _, srb, _ = ring.basis(b, w, dT)
    M = np.array([[0.0, srb[0], srb[1]], [-srd[0], sra[0], sra[1]], [-ud[0], ura[0], ura[1]]])
    rhs = np.array([-srb[2], srd[2] - sra[2], delta - ura[2] + ud[2]])
    return np.linalg.solve(M, rhs)


def liftoff_speed(a, b, disc: Body, ring: Body, delta, dT=0.0) -> float:
    """Lift-off speed (rad/s) of the disc + ring: the speed at which the interface pressure of the
    coupled solution falls to zero, at uniform temperature change ``dT`` from the fitted state.

    By superposition p(w) = p0 - c w^2 with p0 = p(0; delta, dT) and -c = p(1 rad/s; 0, 0) (pure
    centrifugal load), so w_lift = sqrt(p0 / c). Returns 0 if the parts are already apart at
    standstill (p0 <= 0) and ``inf`` if the pressure does not fall with speed (c <= 0)."""
    p0 = coupled_pressure(a, b, disc, ring, 0.0, delta, dT)
    c = -coupled_pressure(a, b, disc, ring, 1.0, 0.0, 0.0)
    if p0 <= 0:
        return 0.0
    if c <= 0:
        return math.inf
    return math.sqrt(p0 / c)


def von_mises_plane_stress(sr, st):
    return np.sqrt(sr * sr - sr * st + st * st)


def summary(a, b, disc: Body, ring: Body, w, delta, n=2001, dT=0.0) -> dict:
    """Peak quantities comparable with the FE results file."""
    sol = disc_ring(a, b, disc, ring, w, delta, dT)
    rd = np.linspace(1e-6 * a, a, n)
    rr = np.linspace(a, b, n)
    _, srd, std = sol["disc"](rd)
    _, srr, str_ = sol["ring"](rr)
    vm = von_mises_plane_stress(srd, std)
    return dict(svm_disc_max=float(vm.max()), sh_ring_max=float(str_.max()), sh_ring_at_bore=float(str_[0]),
                sr_ring_min=float(srr.min()), sr_ring_max=float(srr.max()), p_interface=sol["p_interface"],
                in_contact=sol["in_contact"])
