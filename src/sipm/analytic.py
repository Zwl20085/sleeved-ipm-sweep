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

Layered sleeve (:func:`disc_rings`): the disc with N perfectly bonded rings, each with its own
constants, thermal strain and a uniform hoop eigenstrain (stress-free hoop strain -eps0, the
winding prestress of the deck). Each ring has the solution above plus the eigenstrain's
particular term (:func:`free_strain_terms`, the thermal term with the free strains replaced);
the 2N + 1 constants follow from the traction-free OD, the contact conditions at the disc and
continuity of u and sigma_r at every ring interface. N = 1 is :func:`disc_ring`, same linear
system. :func:`sequential_winding` converts a layer-by-layer winding into equivalent
simultaneous eigenstrains.

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


# ============================================================================================
# disc + N bonded rings with hoop eigenstrain
# ============================================================================================
def free_strain_terms(body: Body, r, f_r: float, f_t: float):
    """(u, sigma_r, sigma_t) of the particular solution for uniform stress-free strains ``f_r``
    (radial) and ``f_t`` (hoop) in ``body`` (no rotation): the thermal term of the module docstring
    with alpha_r dT, alpha_t dT replaced by f_r, f_t. Any homogeneous part is left to A, B."""
    r = np.asarray(r, float)
    k = body.k
    b_r = body.Q11 * f_r + body.Q12 * f_t
    b_t = body.Q12 * f_r + body.Q22 * f_t
    c = (b_r - b_t) / (body.Q11 * (1 + k))
    if c:
        L = _log_ratio_term(r, k)
        u, d = c * r * L, c * (k * L + 1.0)
    else:
        u, d = np.zeros_like(r), np.zeros_like(r)
    sr = body.Q11 * d + body.Q12 * u / r - b_r
    st = body.Q12 * d + body.Q22 * u / r - b_t
    return u, sr, st


def hoop_modulus(body: Body) -> float:
    """Plane-stress hoop modulus of a body, sigma_t / eps_t at sigma_r = 0 (= E_t of
    :func:`orthotropic`)."""
    return body.Q22 - body.Q12 * body.Q12 / body.Q11


def _ring_basis(body: Body, r, w, dT, eps0):
    """:meth:`Body.basis` plus, in the particular column, a uniform hoop eigenstrain: stress-free
    hoop strain -``eps0`` (added only when eps0 != 0, so eps0 = 0 is :meth:`Body.basis` exactly)."""
    u, sr, st = body.basis(r, w, dT)
    if eps0:
        ue, sre, ste = free_strain_terms(body, r, 0.0, -eps0)
        u, sr, st = u.copy(), sr.copy(), st.copy()
        u[2] = u[2] + ue
        sr[2] = sr[2] + sre
        st[2] = st[2] + ste
    return u, sr, st


class RingStack:
    """Solution of a solid disc (0..radii[0]) with N perfectly bonded rings (radii[i-1]..radii[i]).

    Built by :func:`disc_rings`. ``disc(r)`` and ``ring(i, r)`` (i = 0 .. N-1, from the bore outward)
    give (u, sigma_r, sigma_t); ``p_interface`` is the disc/ring contact pressure (0 if the parts
    separate, ``in_contact`` False), ``sr_interfaces`` the radial stress at the N - 1 ring/ring
    interfaces (perfectly bonded: continuous u and sigma_r), ``coeffs`` = [A_disc, A_1, B_1, ...]."""

    def __init__(self, radii, disc, rings, w, dT, eps0, coeffs, p_interface, in_contact):
        self.radii, self.disc_body, self.rings = tuple(radii), disc, tuple(rings)
        self.w, self.dT, self.eps0 = w, dT, tuple(eps0)
        self.coeffs = coeffs
        self.p_interface, self.in_contact = float(p_interface), bool(in_contact)

    @property
    def n(self) -> int:
        return len(self.rings)

    def disc(self, r):
        u, sr, st = self.disc_body.basis(r, self.w, self.dT)
        x = self.coeffs
        return (u[0] * x[0] + u[2], sr[0] * x[0] + sr[2], st[0] * x[0] + st[2])

    def ring(self, i: int, r):
        u, sr, st = _ring_basis(self.rings[i], r, self.w, self.dT, self.eps0[i])
        A, B = self.coeffs[1 + 2 * i], self.coeffs[2 + 2 * i]
        return (u[0] * A + u[1] * B + u[2], sr[0] * A + sr[1] * B + sr[2], st[0] * A + st[1] * B + st[2])

    @property
    def sr_interfaces(self) -> list[float]:
        return [float(self.ring(i, self.radii[i + 1])[1]) for i in range(self.n - 1)]

    def mean_hoop(self, i: int, n: int = 0) -> float:
        """Thickness-average hoop stress of ring i. Without rotation it follows exactly from
        equilibrium, int sigma_t dr = [r sigma_r] over the ring; with rotation (or ``n`` > 0) by
        Gauss-Legendre quadrature with ``n`` (default 32) points."""
        r0, r1 = self.radii[i], self.radii[i + 1]
        if self.w == 0 and n == 0:
            s0, s1 = self.ring(i, r0)[1], self.ring(i, r1)[1]
            return float((r1 * s1 - r0 * s0) / (r1 - r0))
        xg, wg = np.polynomial.legendre.leggauss(n or 32)
        rr = 0.5 * (r1 - r0) * xg + 0.5 * (r1 + r0)
        return float(0.5 * np.dot(wg, self.ring(i, rr)[2]))


def _solve_rings(radii, disc: Body, rings, w, delta, dT, eps0, couple: bool):
    """Coefficients [A_disc, A_1, B_1, ..., A_N, B_N]. Rows: ring-N OD traction free; at the disc
    either sigma_r continuity and u_ring1 - u_disc = delta (``couple``) or a free disc OD and a free
    ring bore; then sigma_r and u continuity at every ring/ring interface. For N = 1 the system is
    the one of :func:`disc_ring`, entry by entry."""
    N = len(rings)
    n = 2 * N + 1
    a = radii[0]
    ud, srd, _ = disc.basis(a, w, dT)
    ura, sra, _ = _ring_basis(rings[0], a, w, dT, eps0[0])
    _, srb, _ = _ring_basis(rings[-1], radii[-1], w, dT, eps0[-1])
    M = np.zeros((n, n))
    rhs = np.zeros(n)
    cN = 2 * N - 1
    M[0, cN], M[0, cN + 1] = srb[0], srb[1]
    rhs[0] = -srb[2]
    if couple:
        M[1, 0] = -srd[0]
        M[1, 1], M[1, 2] = sra[0], sra[1]
        rhs[1] = srd[2] - sra[2]
        M[2, 0] = -ud[0]
        M[2, 1], M[2, 2] = ura[0], ura[1]
        rhs[2] = delta - ura[2] + ud[2]
    else:
        M[1, 0] = srd[0]
        rhs[1] = -srd[2]
        M[2, 1], M[2, 2] = sra[0], sra[1]
        rhs[2] = -sra[2]
    row = 3
    for i in range(N - 1):
        r = radii[i + 1]
        ui, si, _ = _ring_basis(rings[i], r, w, dT, eps0[i])
        uo, so, _ = _ring_basis(rings[i + 1], r, w, dT, eps0[i + 1])
        ci, co = 1 + 2 * i, 3 + 2 * i
        for vin, vout in ((si, so), (ui, uo)):
            M[row, ci], M[row, ci + 1] = vin[0], vin[1]
            M[row, co], M[row, co + 1] = -vout[0], -vout[1]
            rhs[row] = vout[2] - vin[2]
            row += 1
    return np.linalg.solve(M, rhs), sra


def _check_stack(radii, rings, prestrain):
    radii = [float(x) for x in radii]
    N = len(rings)
    if N < 1 or len(radii) != N + 1:
        raise ValueError("disc_rings: need N >= 1 rings and N + 1 radii (disc OD, then each ring OD)")
    if not (radii[0] > 0 and all(r1 > r0 for r0, r1 in zip(radii, radii[1:]))):
        raise ValueError("disc_rings: radii must be positive and increasing")
    eps0 = [0.0] * N if prestrain is None else [float(e) for e in prestrain]
    if len(eps0) != N:
        raise ValueError("disc_rings: need one prestrain per ring")
    return radii, eps0


def disc_rings(radii, disc: Body, rings, w, delta, dT=0.0, prestrain=None) -> RingStack:
    """Solve the disc (0..radii[0]) + N bonded rings problem (plane stress), N = len(rings) >= 1.

    ``radii`` = [disc OD = bore of ring 1, OD of ring 1, ..., OD of ring N] (mm). ``delta`` is the
    radial interference between the disc and ring 1 at the reference temperature, ``dT`` the uniform
    temperature change of all bodies, ``w`` the speed (rad/s). ``prestrain`` gives each ring a uniform
    hoop eigenstrain eps0_i >= 0: its stress-free hoop strain is -eps0_i (it is "too short" by
    eps0_i), so held at zero hoop strain with no radial stress it carries the hoop stress
    E_t eps0_i (E_t = :func:`hoop_modulus`). The rings are perfectly bonded to each other
    (continuous u and sigma_r); the disc/ring-1 interface is a contact that separates if its
    pressure would be negative, as in :func:`disc_ring`. N = 1 without prestrain is
    :func:`disc_ring` (same linear system, the same numbers to the last bit)."""
    radii, eps0 = _check_stack(radii, rings, prestrain)
    x, sra = _solve_rings(radii, disc, rings, w, delta, dT, eps0, True)
    p_int = -(sra[0] * x[1] + sra[1] * x[2] + sra[2])
    contact = p_int >= 0
    if not contact:
        x, _ = _solve_rings(radii, disc, rings, w, delta, dT, eps0, False)
        p_int = 0.0
    return RingStack(radii, disc, rings, w, dT, eps0, x, p_int, contact)


def rings_coupled_pressure(radii, disc: Body, rings, w, delta, dT=0.0, prestrain=None) -> float:
    """Signed disc/ring-1 pressure of the coupled (always in contact) N-ring solution, MPa
    (:func:`coupled_pressure` for N rings)."""
    radii, eps0 = _check_stack(radii, rings, prestrain)
    x, sra = _solve_rings(radii, disc, rings, w, delta, dT, eps0, True)
    return float(-(sra[0] * x[1] + sra[1] * x[2] + sra[2]))


def rings_liftoff_speed(radii, disc: Body, rings, delta, dT=0.0, prestrain=None) -> float:
    """Lift-off speed (rad/s) of the disc + N rings (:func:`liftoff_speed` for N rings): p(w) =
    p0 - c w^2 with p0 from interference, temperature and prestrain at standstill and c from the
    centrifugal load alone."""
    p0 = rings_coupled_pressure(radii, disc, rings, 0.0, delta, dT, prestrain)
    c = -rings_coupled_pressure(radii, disc, rings, 1.0, 0.0, 0.0, None)
    if p0 <= 0:
        return 0.0
    if c <= 0:
        return math.inf
    return math.sqrt(p0 / c)


def rings_summary(radii, disc: Body, rings, w, delta, dT=0.0, prestrain=None, n=2001) -> dict:
    """Peak quantities of the N-ring solution: ``svm_disc_max``, ``p_interface``, ``in_contact``,
    ``sh_ring_max`` (over all rings), ``sr_interfaces`` and per ring (``layers``) the max / min hoop
    stress, the hoop stress at its bore and OD, the min / max radial stress and the mean hoop
    stress."""
    sol = disc_rings(radii, disc, rings, w, delta, dT, prestrain)
    a = sol.radii[0]
    rd = np.linspace(1e-6 * a, a, n)
    _, srd, std = sol.disc(rd)
    layers = []
    for i in range(sol.n):
        rr = np.linspace(sol.radii[i], sol.radii[i + 1], n)
        _, sr, st = sol.ring(i, rr)
        layers.append(dict(sh_max=float(st.max()), sh_min=float(st.min()), sh_bore=float(st[0]),
                           sh_od=float(st[-1]), sr_min=float(sr.min()), sr_max=float(sr.max()),
                           sh_mean=sol.mean_hoop(i)))
    return dict(svm_disc_max=float(von_mises_plane_stress(srd, std).max()), p_interface=sol.p_interface,
                in_contact=sol.in_contact, sh_ring_max=max(d["sh_max"] for d in layers),
                sr_interfaces=sol.sr_interfaces, layers=layers)


# ============================================================================================
# sequential winding -> equivalent simultaneous eigenstrains
# ============================================================================================
def sequential_winding(radii, disc: Body, rings, sigma_w, sublayers: int = 1) -> dict:
    """Equivalent simultaneous hoop eigenstrains of a sleeve wound layer by layer under tension.

    Build-up model (no rotation, reference temperature, no interference): ring j is laid onto the
    disc and the rings below it, which are already stressed by the earlier rings, and has the
    thickness-mean hoop stress ``sigma_w[j]`` (MPa) the moment it is laid, i.e. after the structure
    beneath has given way under it. Each lay-up step is an incremental problem on the structure
    built so far (disc + rings 1..j, displacements counted from the moment ring j is laid), loaded
    by a uniform hoop eigenstrain in ring j alone whose size makes ring j's mean hoop stress
    ``sigma_w[j]``; the later steps relax it. The final state is the sum of the steps. With
    ``sublayers`` = m > 1 every layer is laid as m equal sub-rings in turn, each with the layer's
    ``sigma_w`` (closer to continuous winding).

    The simultaneous model of the deck (:func:`disc_rings` with ``prestrain``) applies one uniform
    eigenstrain per layer to the finished stack at once. The build-up is not exactly one of those
    states (a ring laid on an already displaced surface is equivalent to a misfit at its bore, not
    to a uniform hoop eigenstrain), so the equivalent eigenstrains are chosen to reproduce the
    resultant hoop force of every layer. By equilibrium (int sigma_t dr = [r sigma_r] over a ring
    without rotation) that is the same as reproducing the radial stress at the disc interface and
    at every layer interface, so the contact pressure on the disc and the interlaminar pressures are
    exact; what is approximate is only the hoop-stress gradient through each layer (``max_hoop_dev``).

    ``sigma_w``: one value or one per layer. Returns a dict: ``prestrain`` (eps0 per layer),
    ``prestress`` (eps0 x E_t per layer, MPa: the values for ``sleeve_layer_prestress``),
    ``sequential`` (the build-up: ``sr`` at [radii[0], ..., radii[N-1]], per-layer ``sh_bore``,
    ``sh_od`` and ``sh_mean``, and ``steps``, the sub-ring eigenstrains), ``simultaneous`` (the
    :class:`RingStack` of the equivalent eigenstrains) and ``max_hoop_dev`` (largest |difference| of
    the hoop stress between the two over the sampled radii, MPa).
    """
    radii, _ = _check_stack(radii, rings, None)
    N = len(rings)
    sw = np.broadcast_to(np.asarray(sigma_w, float), (N,)) if np.ndim(sigma_w) == 0 else np.asarray(sigma_w, float)
    if sw.shape != (N,):
        raise ValueError("sequential_winding: give one winding stress or one per layer")
    m = int(sublayers)
    if m < 1:
        raise ValueError("sequential_winding: sublayers must be >= 1")
    # sub-rings of the build-up
    sub_r, sub_b, sub_layer = [radii[0]], [], []
    for i in range(N):
        r0, r1 = radii[i], radii[i + 1]
        for q in range(1, m + 1):
            sub_r.append(r0 + (r1 - r0) * q / m if q < m else r1)
            sub_b.append(rings[i])
            sub_layer.append(i)
    M = len(sub_b)
    steps = []
    for j in range(M):
        rr, bb = sub_r[:j + 2], sub_b[:j + 1]
        unit = disc_rings(rr, disc, bb, 0.0, 0.0, 0.0, [0.0] * j + [1.0])
        e = sw[sub_layer[j]] / unit.mean_hoop(j)
        steps.append(disc_rings(rr, disc, bb, 0.0, 0.0, 0.0, [0.0] * j + [e]))

    def seq_field(r, layer):
        """Final build-up (sigma_r, sigma_t) at radii r inside ``layer`` (sum of the steps)."""
        r = np.atleast_1d(np.asarray(r, float))
        sr, st = np.zeros_like(r), np.zeros_like(r)
        # each radius belongs to one sub-ring (a sub-ring boundary to the lower one: the hoop stress
        # of the build-up jumps there)
        bounds = np.asarray(sub_r[layer * m: layer * m + m + 1])
        loc = np.clip(np.searchsorted(bounds, r, side="left") - 1, 0, m - 1)
        for ql in range(m):
            q = layer * m + ql
            msk = loc == ql
            if not msk.any():
                continue
            for j in range(q, M):
                _, a, b = steps[j].ring(q, r[msk])
                sr[msk] += a
                st[msk] += b
        return sr, st

    seq_sr = []
    for i in range(N):
        seq_sr.append(float(seq_field(radii[i], i)[0][0]))
    # influence of a unit eigenstrain of each layer on the radial stress at [radii[0] .. radii[N-1]]
    G = np.zeros((N, N))
    for k in range(N):
        unit = disc_rings(radii, disc, rings, 0.0, 0.0, 0.0, [1.0 if i == k else 0.0 for i in range(N)])
        G[:, k] = [float(unit.ring(i, radii[i])[1]) for i in range(N)]
    eps = np.linalg.solve(G, np.asarray(seq_sr))
    sim = disc_rings(radii, disc, rings, 0.0, 0.0, 0.0, eps)
    layers, dev = [], 0.0
    for i in range(N):
        rr = np.linspace(radii[i], radii[i + 1], 8 * m + 1)
        _, st_seq = seq_field(rr, i)
        _, _, st_sim = sim.ring(i, rr)
        dev = max(dev, float(np.max(np.abs(st_seq - st_sim))))
        t = radii[i + 1] - radii[i]
        r_sr = (radii[i + 1] * (seq_sr[i + 1] if i + 1 < N else 0.0) - radii[i] * seq_sr[i]) / t
        layers.append(dict(sh_bore=float(st_seq[0]), sh_od=float(st_seq[-1]), sh_mean=float(r_sr)))
    return dict(prestrain=[float(e) for e in eps],
                prestress=[float(e * hoop_modulus(b)) for e, b in zip(eps, rings)],
                sequential=dict(sr=seq_sr, layers=layers, steps=[s.eps0[-1] for s in steps]),
                simultaneous=sim, max_hoop_dev=dev)


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
