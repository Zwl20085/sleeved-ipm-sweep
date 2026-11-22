"""Composite sleeve properties: micromechanics, +/-theta winding homogenisation and ply failure.

Pure numerics (numpy only). The chain is

    Fibre + Matrix + vf  --lamina()-->  Lamina  --angle_ply()-->  Laminate  --sleeve_params()-->  dict

and ``sleeve_params`` emits the ``Params`` sleeve fields (``sleeve_E_hoop`` ... ``sleeve_density``) plus
``sleeve_cte_hoop``, ``sleeve_cte_radial`` and ``sleeve_cte_axial``. A layered sleeve is a list of
``Laminate`` objects; ``layer_stack()`` gives their membrane homogenisation (``LayerStack``), the source
of the ``Params`` constants of a layered sleeve.

Units: moduli in GPa, strengths and stresses in MPa, density in kg/m3, CTE in 1/K, angles in degrees.

No vendor data are shipped (same policy as ``sipm.materials``). The presets in ``FIBRES``, ``MATRICES``
and ``STRENGTHS`` are *illustrative, class-typical round numbers* for a fibre class (standard-, intermediate-
and high-modulus carbon, E-glass) or a matrix class (epoxy, PEEK-class thermoplastic). They are not the
datasheet values of any commercial grade. Supply your own ``Fibre``/``Matrix``/``Strength`` for real work.

Conventions (all of them; nothing else is assumed anywhere in this module)
-------------------------------------------------------------------------
Poisson ratios
    ``nu_ab = -eps_b / eps_a`` under a uniaxial stress along ``a``, so the compliance term is
    ``S[a, b] = -nu_ab / E_a``. For a lamina ``nu12`` is the major ratio (load along the fibre). This is
    the convention of ``Params`` (``sleeve_nu_hoop_radial = -eps_r/eps_hoop`` under hoop stress) and of the
    MAPDL deck (PRYZ = ``sleeve_nu_hoop_axial``, PRXZ = ``sleeve_nu_radial_axial`` with X radial, Y hoop,
    Z axial), so ``sleeve_nu_radial_axial = -eps_axial/eps_radial`` under radial stress.
Ply (material) axes
    1 = fibre, 2 = transverse in the plane of the winding, 3 = through the thickness. Voigt order
    ``(11, 22, 33, 23, 13, 12)``; shear strains are engineering strains (gamma = 2 eps).
Cylinder axes
    Voigt order ``(hoop, axial, radial, axial-radial, hoop-radial, hoop-axial)``, i.e. the usual laminate
    order with x = hoop, y = axial, z = radial. (hoop, axial, radial) is a right-handed triad.
Winding angle
    ``theta`` is measured from the hoop direction in the hoop-axial plane; ``theta = 0`` is a pure hoop
    winding (fibre along hoop, ply 2 along axial, ply 3 along radial), ``theta = 90`` puts the fibre along
    the axis. A *positive* angle turns the fibre from +hoop towards +axial (right-handed about +radial).
    For a balanced +/-theta winding the effective constants do not depend on the sign; the sign only labels
    the two plies in ``Laminate.ply_stresses``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "Fibre", "Matrix", "Strength", "Lamina", "Laminate", "LayerStack",
    "lamina", "angle_ply", "layer_stack", "sleeve_params", "failure_indices", "ply_failure", "sleeve_failure",
    "FIBRES", "MATRICES", "STRENGTHS", "MODELS", "SLEEVE_KEYS", "CRITERIA",
]

MODELS = ("chamis", "halpin_tsai", "rom")

#: Keys produced by :func:`sleeve_params` (ten existing ``Params`` fields plus three CTE fields).
SLEEVE_KEYS = (
    "sleeve_E_hoop", "sleeve_E_radial", "sleeve_E_axial",
    "sleeve_nu_hoop_radial", "sleeve_nu_hoop_axial", "sleeve_nu_radial_axial",
    "sleeve_G_hoop_radial", "sleeve_G_hoop_axial", "sleeve_G_radial_axial",
    "sleeve_density", "sleeve_cte_hoop", "sleeve_cte_radial", "sleeve_cte_axial",
)


# ---------------------------------------------------------------------------------------------------------
# constituents
# ---------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Fibre:
    """Transversely isotropic fibre, axis 1 along the fibre.

    E1, E2, G12 in GPa; nu12 = -eps2/eps1 under axial stress (major); density kg/m3; cte1 (axial) and
    cte2 (transverse) in 1/K. Give exactly one of ``G23`` (GPa) or ``nu23``; the other follows from
    transverse isotropy, G23 = E2 / (2 (1 + nu23)). An isotropic fibre has E1 = E2, G12 = E/(2(1+nu)),
    nu23 = nu12 and cte1 = cte2.
    """

    E1: float
    E2: float
    G12: float
    nu12: float
    density: float
    cte1: float
    cte2: float
    G23: float | None = None
    nu23: float | None = None

    def __post_init__(self):
        if (self.G23 is None) == (self.nu23 is None):
            raise ValueError("Fibre: give exactly one of G23 or nu23")
        if self.G23 is None:
            object.__setattr__(self, "G23", self.E2 / (2.0 * (1.0 + self.nu23)))
        else:
            object.__setattr__(self, "nu23", self.E2 / (2.0 * self.G23) - 1.0)
        for n in ("E1", "E2", "G12", "G23", "density"):
            if not getattr(self, n) > 0:
                raise ValueError(f"Fibre.{n} must be > 0")
        if not -1.0 < self.nu23 < 1.0:
            raise ValueError(f"Fibre: nu23 = {self.nu23:g} is outside (-1, 1)")


@dataclass(frozen=True)
class Matrix:
    """Isotropic matrix: E in GPa, nu, density kg/m3, cte 1/K."""

    E: float
    nu: float
    density: float
    cte: float

    def __post_init__(self):
        if not (self.E > 0 and self.density > 0 and -1.0 < self.nu < 0.5):
            raise ValueError("Matrix: need E > 0, density > 0 and -1 < nu < 0.5")

    @property
    def G(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))


@dataclass(frozen=True)
class Strength:
    """Unidirectional ply strengths in MPa, all given as positive magnitudes.

    Xt, Xc: fibre-direction tension / compression; Yt, Yc: transverse tension / compression; S: in-plane
    (12) shear. ``St`` is the transverse (23) shear strength used only by the Hashin matrix-compression
    branch; when it is None the in-plane value S is used (an assumption, see :func:`failure_indices`).
    """

    Xt: float
    Xc: float
    Yt: float
    Yc: float
    S: float
    St: float | None = None

    def __post_init__(self):
        for n in ("Xt", "Xc", "Yt", "Yc", "S"):
            if not getattr(self, n) > 0:
                raise ValueError(f"Strength.{n} must be a positive magnitude")
        if self.St is not None and not self.St > 0:
            raise ValueError("Strength.St must be > 0")


@dataclass(frozen=True)
class Lamina:
    """Unidirectional ply in material axes (1 fibre, 2 in-plane transverse, 3 through-thickness).

    Nine orthotropic constants (GPa; ``nu_ij = -eps_j/eps_i`` under stress along i), density kg/m3,
    CTEs in 1/K and the fibre volume fraction it was built at (``nan`` if entered directly).
    """

    E1: float
    E2: float
    E3: float
    nu12: float
    nu13: float
    nu23: float
    G12: float
    G13: float
    G23: float
    density: float
    cte1: float
    cte2: float
    cte3: float
    vf: float = math.nan

    def compliance(self) -> np.ndarray:
        """6x6 compliance [1/GPa], Voigt order (11, 22, 33, 23, 13, 12), engineering shear strain."""
        S = np.zeros((6, 6))
        S[0, 0], S[1, 1], S[2, 2] = 1 / self.E1, 1 / self.E2, 1 / self.E3
        S[0, 1] = S[1, 0] = -self.nu12 / self.E1
        S[0, 2] = S[2, 0] = -self.nu13 / self.E1
        S[1, 2] = S[2, 1] = -self.nu23 / self.E2
        S[3, 3], S[4, 4], S[5, 5] = 1 / self.G23, 1 / self.G13, 1 / self.G12
        return S

    def stiffness(self) -> np.ndarray:
        """6x6 stiffness [GPa], inverse of :meth:`compliance`."""
        return np.linalg.inv(self.compliance())

    @property
    def alpha(self) -> np.ndarray:
        """CTE as a Voigt strain vector [1/K]."""
        return np.array([self.cte1, self.cte2, self.cte3, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------------------------------------
# micromechanics
# ---------------------------------------------------------------------------------------------------------
def _halpin_tsai(Mf: float, Mm: float, vf: float, xi: float) -> float:
    eta = (Mf / Mm - 1.0) / (Mf / Mm + xi)
    return Mm * (1.0 + xi * eta * vf) / (1.0 - eta * vf)


def _chamis(Mf: float, Mm: float, vf: float) -> float:
    return Mm / (1.0 - math.sqrt(vf) * (1.0 - Mm / Mf))


def _reuss(Mf: float, Mm: float, vf: float) -> float:
    return 1.0 / (vf / Mf + (1.0 - vf) / Mm)


def lamina(fibre: Fibre, matrix: Matrix, vf: float, model: str = "chamis") -> Lamina:
    """Unidirectional ply properties from fibre, matrix and fibre volume fraction ``vf`` (0..1).

    All models share (vm = 1 - vf, Voigt rule of mixtures, exact in the 1D strain-compatible limit)::

        E1   = vf Ef1  + vm Em
        nu12 = vf nuf12 + vm num
        rho  = vf rhof + vm rhom

    and transverse isotropy of the ply: E3 = E2, nu13 = nu12, G13 = G12, nu23 = E2 / (2 G23) - 1.

    ``model`` selects E2, G12, G23:

    ``"chamis"``  (C.C. Chamis, "Simplified composite micromechanics equations for hygral, thermal and
        mechanical properties", SAMPE Quarterly 15(3):14-23, 1984; NASA TM-83320, 1983)::

            M = Mm / (1 - sqrt(vf) (1 - Mm/Mf))      with (Mf, Mm) = (Ef2, Em), (Gf12, Gm), (Gf23, Gm)

    ``"halpin_tsai"``  (J.C. Halpin, J.L. Kardos, "The Halpin-Tsai equations: a review", Polym. Eng. Sci.
        16(5):344-352, 1976)::

            M / Mm = (1 + xi eta vf) / (1 - eta vf),   eta = (Mf/Mm - 1) / (Mf/Mm + xi)

        with xi = 2 for E2 and xi = 1 for G12 (circular fibres, square array). Halpin-Tsai has no
        canonical xi for G23; this module uses xi = 1 there too (a choice, stated here).

    ``"rom"``  Voigt for E1 and nu12 (above), Reuss (series) for E2, G12 (Gf12) and G23 (Gf23)::

            1/M = vf/Mf + vm/Mm

    Thermal expansion (all models; R.A. Schapery, "Thermal expansion coefficients of composite materials
    based on energy principles", J. Compos. Mater. 2(3):380-404, 1968, written for a transversely isotropic
    fibre; the fibre-anisotropy terms follow from the same 1D argument: common axial strain, zero transverse
    stress in each phase)::

        cte1 = (vf Ef1 af1 + vm Em am) / (vf Ef1 + vm Em)
        cte2 = vf (af2 + nuf12 af1) + vm (1 + num) am - cte1 nu12,       cte3 = cte2

    where nu12 is the rule-of-mixtures value above (Schapery's), whatever ``model`` is.
    """
    if model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}, got {model!r}")
    vf = float(vf)
    if not 0.0 <= vf <= 1.0:
        raise ValueError(f"vf must be in [0, 1], got {vf}")
    vm = 1.0 - vf
    f, m = fibre, matrix
    Gm = m.G

    E1 = vf * f.E1 + vm * m.E
    nu12 = vf * f.nu12 + vm * m.nu
    rho = vf * f.density + vm * m.density
    if model == "chamis":
        E2, G12, G23 = _chamis(f.E2, m.E, vf), _chamis(f.G12, Gm, vf), _chamis(f.G23, Gm, vf)
    elif model == "halpin_tsai":
        E2, G12, G23 = _halpin_tsai(f.E2, m.E, vf, 2.0), _halpin_tsai(f.G12, Gm, vf, 1.0), \
            _halpin_tsai(f.G23, Gm, vf, 1.0)
    else:
        E2, G12, G23 = _reuss(f.E2, m.E, vf), _reuss(f.G12, Gm, vf), _reuss(f.G23, Gm, vf)
    nu23 = E2 / (2.0 * G23) - 1.0

    a1 = (vf * f.E1 * f.cte1 + vm * m.E * m.cte) / (vf * f.E1 + vm * m.E)
    a2 = vf * (f.cte2 + f.nu12 * f.cte1) + vm * (1.0 + m.nu) * m.cte - a1 * nu12

    lam = Lamina(E1=E1, E2=E2, E3=E2, nu12=nu12, nu13=nu12, nu23=nu23, G12=G12, G13=G12, G23=G23,
                 density=rho, cte1=a1, cte2=a2, cte3=a2, vf=vf)
    _require_spd(lam.compliance(), f"lamina({model!r}, vf={vf:g})")
    return lam


def _require_spd(S: np.ndarray, what: str) -> None:
    if np.abs(S - S.T).max() > 1e-10 * np.abs(S).max() or np.linalg.eigvalsh(0.5 * (S + S.T)).min() <= 0:
        raise ValueError(f"{what}: compliance is not symmetric positive definite")


# ---------------------------------------------------------------------------------------------------------
# rotation and +/-theta homogenisation
# ---------------------------------------------------------------------------------------------------------
_VOIGT = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
_I = [0, 1, 5]   # in-plane components (hoop, axial, hoop-axial): strains common to the plies
_O = [2, 3, 4]   # out-of-plane components (radial, axial-radial, hoop-radial): stresses common to the plies


def _bond(Q: np.ndarray) -> np.ndarray:
    """6x6 stress transformation M with sigma' = M sigma for the tensor rotation sigma' = Q sigma Q^T.

    Built column by column from the tensor rule, so it carries no transcription risk. With engineering
    shear strains the strain transformation is inv(M).T and stiffness transforms as C' = M C M^T.
    """
    M = np.empty((6, 6))
    for b, (k, l) in enumerate(_VOIGT):
        T = np.zeros((3, 3))
        T[k, l] = T[l, k] = 1.0
        Tp = Q @ T @ Q.T
        M[:, b] = [Tp[i, j] for i, j in _VOIGT]
    return M


def _ply_to_cyl(theta_deg: float) -> np.ndarray:
    """Rotation whose columns are the ply axes (1, 2, 3) in cylinder axes (hoop, axial, radial)."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _mixed(C: np.ndarray, beta: np.ndarray):
    """Partial inversion of sigma = C eps - beta dT into [sig_I; eps_O] = H [eps_I; sig_O] + g dT."""
    Coo_inv = np.linalg.inv(C[np.ix_(_O, _O)])
    Cio, Coi = C[np.ix_(_I, _O)], C[np.ix_(_O, _I)]
    H = np.block([[C[np.ix_(_I, _I)] - Cio @ Coo_inv @ Coi, Cio @ Coo_inv],
                  [-Coo_inv @ Coi, Coo_inv]])
    g = np.concatenate([Cio @ Coo_inv @ beta[_O] - beta[_I], Coo_inv @ beta[_O]])
    return H, g


def _invert_mixed(Hb: np.ndarray, gb: np.ndarray):
    """Invert an averaged mixed law [sig_I; eps_O] = Hb [eps_I; sig_O] + gb dT back to
    sigma = C (eps - alpha dT). Returns (C, alpha) with C symmetrised, in the unit of Hb."""
    hii, hio, hoi, hoo = Hb[:3, :3], Hb[:3, 3:], Hb[3:, :3], Hb[3:, 3:]
    hoo_inv = np.linalg.inv(hoo)
    C = np.zeros((6, 6))
    C[np.ix_(_I, _I)] = hii - hio @ hoo_inv @ hoi
    C[np.ix_(_I, _O)] = hio @ hoo_inv
    C[np.ix_(_O, _I)] = -hoo_inv @ hoi
    C[np.ix_(_O, _O)] = hoo_inv
    C = 0.5 * (C + C.T)
    # free thermal strain: zero averaged in-plane stress and zero (common) out-of-plane stress
    eps_I = -np.linalg.solve(hii, gb[:3])
    eps_O = hoi @ eps_I + gb[3:]
    alpha = np.zeros(6)
    alpha[_I], alpha[_O] = eps_I, eps_O
    return C, alpha


class _CylOrthotropic:
    """Engineering constants (cylinder axes) of an object with a compliance ``S`` [1/GPa] in cylinder
    Voigt order (hoop, axial, radial, axial-radial, hoop-radial, hoop-axial) and a CTE vector ``alpha``,
    with the module's Poisson convention ``nu_ab = -eps_b/eps_a`` under stress along a."""

    S: np.ndarray
    alpha: np.ndarray

    def _nu(self, a: int, b: int) -> float:
        return -self.S[a, b] / self.S[a, a]

    E_hoop = property(lambda self: 1.0 / self.S[0, 0])
    E_axial = property(lambda self: 1.0 / self.S[1, 1])
    E_radial = property(lambda self: 1.0 / self.S[2, 2])
    nu_hoop_axial = property(lambda self: self._nu(0, 1))
    nu_hoop_radial = property(lambda self: self._nu(0, 2))
    nu_radial_axial = property(lambda self: self._nu(2, 1))
    G_radial_axial = property(lambda self: 1.0 / self.S[3, 3])
    G_hoop_radial = property(lambda self: 1.0 / self.S[4, 4])
    G_hoop_axial = property(lambda self: 1.0 / self.S[5, 5])
    cte_hoop = property(lambda self: float(self.alpha[0]))
    cte_axial = property(lambda self: float(self.alpha[1]))
    cte_radial = property(lambda self: float(self.alpha[2]))


class Laminate(_CylOrthotropic):
    """Homogenised balanced +/-theta winding (equal ply thicknesses) in cylinder axes.

    Built by :func:`angle_ply`. Attributes: ``ply`` (the :class:`Lamina`), ``angle_deg``, ``C`` and ``S``
    (6x6 effective stiffness [GPa] / compliance [1/GPa] in cylinder Voigt order
    (hoop, axial, radial, axial-radial, hoop-radial, hoop-axial)), ``alpha`` (effective CTE vector, 1/K)
    and ``density``. The engineering constants are exposed as properties (``E_hoop`` ... ``cte_axial``)
    with the module's Poisson convention ``nu_ab = -eps_b/eps_a`` under stress along a.
    """

    def __init__(self, ply: Lamina, angle_deg: float):
        self.ply = ply
        self.angle_deg = float(angle_deg)
        Cm, am = ply.stiffness() * 1e3, ply.alpha   # internal work in MPa
        self._M, self._H, self._g = {}, {}, {}
        Hs, gs = [], []
        for sgn in (+1, -1):
            M = _bond(_ply_to_cyl(sgn * self.angle_deg))
            C = M @ Cm @ M.T
            a = np.linalg.inv(M).T @ am
            H, g = _mixed(C, C @ a)
            self._M[sgn], self._H[sgn], self._g[sgn] = M, H, g
            Hs.append(H)
            gs.append(g)
        Hb, gb = 0.5 * (Hs[0] + Hs[1]), 0.5 * (gs[0] + gs[1])
        # invert the averaged mixed law back to sigma = C (eps - alpha dT)
        C, alpha = _invert_mixed(Hb, gb)
        S = np.linalg.inv(C)
        self._S_mpa = 0.5 * (S + S.T)
        self.C = C / 1e3
        self.S = self._S_mpa * 1e3
        self.alpha = alpha
        self.density = ply.density
        _require_spd(self.S, f"angle_ply({self.angle_deg:g} deg)")

    # engineering constants (cylinder axes): inherited from _CylOrthotropic

    def ply_stresses(self, sig, delta_T=0.0):
        """Stresses in the +theta and -theta plies, in each ply's material axes.

        ``sig``: homogenised stress, shape (..., 6), cylinder Voigt order
        (hoop, axial, radial, axial-radial, hoop-radial, hoop-axial), any consistent unit (MPa for
        :func:`failure_indices`). ``delta_T``: temperature change [K], scalar or broadcastable to
        ``sig[..., 0]``; when it is not zero, ``sig`` and the result are in MPa.

        The effective strain is eps = S sig + alpha dT. Within each ply the in-plane strains
        (hoop, axial, hoop-axial) equal the effective ones and the out-of-plane stresses (radial and the two
        transverse shears) equal the homogenised ones; the ply's in-plane stresses then follow from its own
        stiffness. The thickness average of the two ply stresses rotated back to cylinder axes is ``sig``.

        Returns ``(sig_plus, sig_minus)``, each shape (..., 6), ply Voigt order (11, 22, 33, 23, 13, 12),
        with 1 the local fibre direction. ``sig_plus`` is the ply at +theta (fibre turned from hoop towards
        +axial). Use components 0, 1, 5 for :func:`failure_indices`.
        """
        sig = np.asarray(sig, dtype=float)
        dT = np.asarray(delta_T, dtype=float)[..., None]
        eps = sig @ self._S_mpa.T + self.alpha * dT
        x = np.concatenate([eps[..., _I], sig[..., _O]], axis=-1)
        out = []
        for sgn in (+1, -1):
            H, g = self._H[sgn], self._g[sgn]
            s_in = x @ H[:3].T + g[:3] * dT
            s_cyl = np.empty_like(sig)
            s_cyl[..., _I] = s_in
            s_cyl[..., _O] = sig[..., _O]
            out.append(s_cyl @ np.linalg.inv(self._M[sgn]).T)   # cylinder -> ply axes
        return out[0], out[1]

    def __repr__(self) -> str:
        return (f"Laminate(+/-{self.angle_deg:g} deg, E_hoop={self.E_hoop:.4g}, E_axial={self.E_axial:.4g}, "
                f"E_radial={self.E_radial:.4g} GPa)")


def angle_ply(lam: Lamina, angle_deg: float) -> Laminate:
    """Effective orthotropic constants of a balanced +/-theta winding in cylinder axes.

    Each ply's 3D stiffness is rotated about the radial axis by +theta and -theta (``C' = M C M^T`` with
    the Bond matrix M of the rotation; CTE as a strain vector). The two plies (equal thickness, stacked
    through the radius) are homogenised with the exact layered-medium average: in-plane strains
    (hoop, axial, hoop-axial) and out-of-plane stresses (radial, axial-radial, hoop-radial) are common to
    the plies, and in-plane stresses and out-of-plane strains are thickness-averaged (C.T. Sun, S. Li,
    "Three-dimensional effective elastic constants for thick laminates", J. Compos. Mater. 22:629-639,
    1988). This is done by partially inverting each ply law to [sig_I; eps_O] = H [eps_I; sig_O] + g dT,
    averaging H and g and inverting back. It is not a plain Voigt average of C: the out-of-plane couplings
    C36 and C45 change sign between the plies and the plain average is wrong by terms of second order in
    them. For +/-theta the monoclinic terms cancel exactly and the result is orthotropic in cylinder axes.
    The effective CTE is the free strain at zero homogenised stress.

    ``angle_deg = 0`` reproduces the lamina (E_hoop = E1, E_axial = E2, E_radial = E3).
    """
    return Laminate(lam, angle_deg)


class LayerStack(_CylOrthotropic):
    """Membrane (thin-wall) homogenisation of a stack of perfectly bonded sleeve layers.

    Built by :func:`layer_stack`. ``layers`` are :class:`Laminate` objects (each a balanced +/- angle
    winding, orthotropic in cylinder axes) stacked through the radius with thickness ``fractions``
    (normalised to sum 1). The layers share the in-plane strains (hoop, axial, hoop-axial) and carry
    the same out-of-plane stresses (radial and the two transverse shears); the other components are
    thickness averages. This is the same layered-medium average :func:`angle_ply` applies to the two
    plies of a winding (Sun & Li 1988), applied once more to the layers. It is exact for a flat
    laminate and, for a cylinder, neglects the variation of the hoop strain u/r through the
    thickness, i.e. it is the description of a sleeve whose thickness is small against its radius.
    The CTE is the free strain at zero homogenised stress and the density the thickness average.

    Attributes as :class:`Laminate` (``C``, ``S``, ``alpha``, ``density`` and the engineering constants),
    but there are no plies: the failure evaluation needs the individual ``layers``.
    """

    def __init__(self, layers, fractions):
        layers, w = list(layers), np.asarray(fractions, dtype=float)
        if not layers or len(layers) != len(w):
            raise ValueError("layer_stack: need one fraction per layer")
        if not all(isinstance(lay, Laminate) for lay in layers):
            raise TypeError("layer_stack: every layer must be a Laminate")
        if not (np.all(np.isfinite(w)) and np.all(w > 0)):
            raise ValueError("layer_stack: fractions must be finite and > 0")
        w = w / w.sum()
        self.layers, self.fractions = layers, tuple(float(x) for x in w)
        Hb, gb = np.zeros((6, 6)), np.zeros(6)
        for lay, wi in zip(layers, w):
            C = lay.C * 1e3   # MPa
            H, g = _mixed(C, C @ lay.alpha)
            Hb += wi * H
            gb += wi * g
        C, alpha = _invert_mixed(Hb, gb)
        S = np.linalg.inv(C)
        S = 0.5 * (S + S.T)
        self.C = C / 1e3
        self.S = S * 1e3
        self.alpha = alpha
        self.density = float(sum(wi * lay.density for lay, wi in zip(layers, w)))
        _require_spd(self.S, "layer_stack")

    def __repr__(self) -> str:
        spec = ", ".join(f"+/-{lay.angle_deg:g}:{f:.3g}" for lay, f in zip(self.layers, self.fractions))
        return f"LayerStack({spec}; E_hoop={self.E_hoop:.4g}, E_axial={self.E_axial:.4g} GPa)"


def layer_stack(layers, fractions) -> LayerStack:
    """Membrane homogenisation of bonded sleeve layers (:class:`LayerStack`); ``layers`` from the bore
    outward, ``fractions`` their thickness shares (normalised)."""
    return LayerStack(layers, fractions)


def sleeve_params(x: Lamina | Laminate | LayerStack) -> dict[str, float]:
    """``Params`` sleeve fields from a winding; a bare :class:`Lamina` is taken as a pure hoop winding,
    a :class:`LayerStack` gives the membrane homogenisation of a layered sleeve.

    Keys: :data:`SLEEVE_KEYS` (GPa, -, kg/m3, 1/K). ``sleeve_nu_hoop_radial = -eps_r/eps_hoop`` under hoop
    stress, ``sleeve_nu_hoop_axial = -eps_axial/eps_hoop`` under hoop stress and
    ``sleeve_nu_radial_axial = -eps_axial/eps_r`` under radial stress, as ``Params`` and the MAPDL deck read
    them. Note that a +/-theta winding can give Poisson ratios outside [0, 0.5) (e.g. nu_hoop_axial > 1
    near +/-30 deg), which is physical for a laminate.
    """
    if isinstance(x, Lamina):
        x = angle_ply(x, 0.0)
    if not isinstance(x, (Laminate, LayerStack)):
        raise TypeError("sleeve_params needs a Lamina, a Laminate or a LayerStack")
    return {
        "sleeve_E_hoop": float(x.E_hoop),
        "sleeve_E_radial": float(x.E_radial),
        "sleeve_E_axial": float(x.E_axial),
        "sleeve_nu_hoop_radial": float(x.nu_hoop_radial),
        "sleeve_nu_hoop_axial": float(x.nu_hoop_axial),
        "sleeve_nu_radial_axial": float(x.nu_radial_axial),
        "sleeve_G_hoop_radial": float(x.G_hoop_radial),
        "sleeve_G_hoop_axial": float(x.G_hoop_axial),
        "sleeve_G_radial_axial": float(x.G_radial_axial),
        "sleeve_density": float(x.density),
        "sleeve_cte_hoop": x.cte_hoop,
        "sleeve_cte_radial": x.cte_radial,
        "sleeve_cte_axial": x.cte_axial,
    }


# ---------------------------------------------------------------------------------------------------------
# failure criteria
# ---------------------------------------------------------------------------------------------------------
def _ratio(a, b):
    """Positive root R of a R^2 + b R - 1 = 0 (a >= 0), written as 2 / (b + sqrt(b^2 + 4a)) for stability.

    R = inf where the load is zero or never reaches the envelope.
    """
    den = b + np.sqrt(b * b + 4.0 * a)
    with np.errstate(divide="ignore"):
        return np.where(den > 0, 2.0 / np.where(den > 0, den, 1.0), np.inf)


def failure_indices(s1, s2, t12, strength: Strength, f12_star: float = -0.5) -> dict[str, np.ndarray]:
    """Plane-stress ply failure criteria (vectorised; MPa; ply axes, 1 = fibre).

    For every criterion ``<name>`` the failure index (``>= 1`` means failure) and ``<name>_sr``, the
    strength ratio (the factor by which the whole stress state may be scaled to reach the envelope, i.e.
    the safety factor) are returned:

    ``max_stress``: max(s1/Xt or -s1/Xc, s2/Yt or -s2/Yc, |t12|/S). Linear, so sr = 1/index.

    ``tsai_wu`` (S.W. Tsai, E.M. Wu, J. Compos. Mater. 5:58-80, 1971)::

        F = F1 s1 + F2 s2 + F11 s1^2 + F22 s2^2 + F66 t12^2 + 2 F12 s1 s2
        F1 = 1/Xt - 1/Xc, F2 = 1/Yt - 1/Yc, F11 = 1/(Xt Xc), F22 = 1/(Yt Yc), F66 = 1/S^2
        F12 = f12_star sqrt(F11 F22),  f12_star = -1/2 by default (Tsai & Hahn, Introduction to Composite
        Materials, 1980; the value that makes the quadric a generalised von Mises ellipse).

    ``tsai_wu`` is the quadratic F itself. ``tsai_wu_sr`` is the positive root R of
    (quadratic part) R^2 + (linear part) R = 1, solved in a cancellation-free form; both are 1 exactly on
    the envelope, but only 1/sr scales linearly with the load.

    ``hashin_fibre`` / ``hashin_matrix`` (Z. Hashin, "Failure criteria for unidirectional fiber composites",
    J. Appl. Mech. 47:329-334, 1980, plane-stress form)::

        fibre,  s1 >= 0: (s1/Xt)^2 + (t12/S)^2        fibre,  s1 < 0: (s1/Xc)^2
        matrix, s2 >= 0: (s2/Yt)^2 + (t12/S)^2
        matrix, s2 < 0:  (s2/(2 St))^2 + ((Yc/(2 St))^2 - 1) s2/Yc + (t12/S)^2

    St is ``strength.St`` or, if None, S (assumption). Their ``_sr`` solve the same quadratic in R as
    Tsai-Wu (= 1/sqrt(index) on the homogeneous branches).
    """
    s1, s2, t12 = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (s1, s2, t12)))
    X = strength
    St = X.S if X.St is None else X.St

    ms = np.maximum.reduce([np.where(s1 >= 0, s1 / X.Xt, -s1 / X.Xc),
                            np.where(s2 >= 0, s2 / X.Yt, -s2 / X.Yc),
                            np.abs(t12) / X.S])
    with np.errstate(divide="ignore"):
        ms_sr = np.where(ms > 0, 1.0 / np.where(ms > 0, ms, 1.0), np.inf)

    F1, F2 = 1 / X.Xt - 1 / X.Xc, 1 / X.Yt - 1 / X.Yc
    F11, F22, F66 = 1 / (X.Xt * X.Xc), 1 / (X.Yt * X.Yc), 1 / X.S ** 2
    F12 = f12_star * math.sqrt(F11 * F22)
    quad = F11 * s1 ** 2 + F22 * s2 ** 2 + F66 * t12 ** 2 + 2 * F12 * s1 * s2
    lin = F1 * s1 + F2 * s2

    shear = (t12 / X.S) ** 2
    hf_q = np.where(s1 >= 0, (s1 / X.Xt) ** 2 + shear, (s1 / X.Xc) ** 2)
    hm_q = np.where(s2 >= 0, (s2 / X.Yt) ** 2 + shear, (s2 / (2 * St)) ** 2 + shear)
    hm_l = np.where(s2 >= 0, 0.0, ((X.Yc / (2 * St)) ** 2 - 1.0) * s2 / X.Yc)

    return {
        "max_stress": ms, "max_stress_sr": ms_sr,
        "tsai_wu": quad + lin, "tsai_wu_sr": _ratio(quad, lin),
        "hashin_fibre": hf_q, "hashin_fibre_sr": _ratio(hf_q, np.zeros_like(hf_q)),
        "hashin_matrix": hm_q + hm_l, "hashin_matrix_sr": _ratio(hm_q, hm_l),
    }


CRITERIA = ("max_stress", "tsai_wu", "hashin_fibre", "hashin_matrix")
# the two transverse planes of a ply, as (s1, s_transverse, t_fibre-transverse) in ply Voigt order
_PLANES = {"12": (0, 1, 5), "13": (0, 2, 4)}


def ply_failure(plies, strength: Strength, f12_star: float = -0.5) -> dict[str, np.ndarray]:
    """Plane-stress failure indices of one or more plies, worst over plies and transverse planes.

    ``plies``: sequence of ply stress arrays, each shape (..., 6), ply Voigt order (11, 22, 33, 23, 13, 12)
    in MPa (e.g. the pair from :meth:`Laminate.ply_stresses`). For each ply the criteria of
    :func:`failure_indices` are evaluated on both transverse planes, (s11, s22, t12) and (s11, s33, t13),
    with the same strengths: this relies on the ply's transverse isotropy (Yt, Yc, S are the same across the
    two transverse directions). It ignores t23 and any interaction of s22 with s33, i.e. it is not a 3D
    criterion. Returns, for every criterion of :data:`CRITERIA`, the index (max over plies and planes) and
    ``<name>_sr`` (min strength ratio over plies and planes), each of shape (...,).
    """
    plies = [np.asarray(s, dtype=float) for s in plies]
    if not plies or any(s.shape[-1] != 6 for s in plies):
        raise ValueError("ply_failure needs one or more arrays of shape (..., 6)")
    out: dict[str, np.ndarray] = {}
    for s in plies:
        for i1, i2, i12 in _PLANES.values():
            fi = failure_indices(s[..., i1], s[..., i2], s[..., i12], strength, f12_star)
            for name in CRITERIA:
                for key, red in ((name, np.maximum), (f"{name}_sr", np.minimum)):
                    out[key] = fi[key] if key not in out else red(out[key], fi[key])
    return out


def sleeve_failure(sig_cyl, strength: Strength, laminate: Laminate | None = None, delta_T=0.0,
                   f12_star: float = -0.5) -> dict[str, np.ndarray]:
    """Failure indices of the sleeve from homogenised stresses (MPa) in cylinder Voigt order
    (hoop, axial, radial, axial-radial, hoop-radial, hoop-axial), shape (..., 6).

    With a ``laminate`` the stresses are split into the +theta and -theta ply stresses
    (:meth:`Laminate.ply_stresses`, including the ply-to-ply thermal mismatch for ``delta_T``); without one
    the sleeve is taken as a single hoop-wound ply (fibre = hoop, ply axes 2 = axial, 3 = radial), whose ply
    stresses are the homogenised ones. Then :func:`ply_failure`. Cure residual stresses are not included.
    """
    sig = np.asarray(sig_cyl, dtype=float)
    if laminate is not None and not isinstance(laminate, Laminate):
        raise TypeError("sleeve_failure needs the Laminate of one winding (a LayerStack has no plies: "
                        "evaluate each layer with its own Laminate)")
    plies = [sig] if laminate is None else list(laminate.ply_stresses(sig, delta_T))
    return ply_failure(plies, strength, f12_star)


# ---------------------------------------------------------------------------------------------------------
# illustrative presets (class-typical round numbers, NOT vendor grades)
# ---------------------------------------------------------------------------------------------------------
#: Illustrative fibre classes. Round numbers typical of the class; not any commercial grade.
FIBRES: dict[str, Fibre] = {
    "carbon_sm": Fibre(E1=230.0, E2=15.0, G12=15.0, nu12=0.20, nu23=0.25, density=1800.0,
                       cte1=-0.5e-6, cte2=7.0e-6),      # standard-modulus PAN carbon
    "carbon_im": Fibre(E1=290.0, E2=15.0, G12=15.0, nu12=0.20, nu23=0.25, density=1800.0,
                       cte1=-0.6e-6, cte2=7.0e-6),      # intermediate-modulus PAN carbon
    "carbon_hm": Fibre(E1=390.0, E2=7.0, G12=15.0, nu12=0.20, nu23=0.30, density=1900.0,
                       cte1=-1.0e-6, cte2=8.0e-6),      # high-modulus carbon
    "e_glass": Fibre(E1=72.0, E2=72.0, G12=72.0 / 2.44, nu12=0.22, nu23=0.22, density=2550.0,
                     cte1=5.0e-6, cte2=5.0e-6),         # isotropic E-glass
}

#: Illustrative matrix classes. Round numbers typical of the class; not any commercial grade.
MATRICES: dict[str, Matrix] = {
    "epoxy": Matrix(E=3.5, nu=0.35, density=1200.0, cte=60e-6),
    "peek": Matrix(E=3.8, nu=0.40, density=1300.0, cte=50e-6),   # PEEK-class semicrystalline thermoplastic
}

#: Illustrative unidirectional ply strengths [MPa] for the class at vf ~ 0.6 with an epoxy matrix.
#: They are NOT computed from vf; scale them yourself if you change vf substantially.
STRENGTHS: dict[str, Strength] = {
    "carbon_sm_epoxy": Strength(Xt=2000.0, Xc=1200.0, Yt=50.0, Yc=200.0, S=70.0),
    "carbon_im_epoxy": Strength(Xt=2500.0, Xc=1400.0, Yt=60.0, Yc=220.0, S=90.0),
    "carbon_hm_epoxy": Strength(Xt=1500.0, Xc=800.0, Yt=40.0, Yc=170.0, S=60.0),
    "e_glass_epoxy": Strength(Xt=1000.0, Xc=600.0, Yt=40.0, Yc=140.0, S=70.0),
}
