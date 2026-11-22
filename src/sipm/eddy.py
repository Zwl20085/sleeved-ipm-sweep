"""Sleeve eddy-current loss: averaging, resolution checks and a thin-shell estimate (numpy only).

Used by the Maxwell driver for the 2D sleeve loss (docs/methods.md, 2.1) and by its validation.
Everything here is licence free and unit tested.

Conventions: 2D, axial (z) current only. In the rotor frame the sleeve is at rest, so its current
density is J_z = sigma E_z with E_z = -dA_z/dt + C(t), A_z the magnetic vector potential seen by a
point fixed in the rotor and C(t) a field uniform over the conductor that the end connection sets.
"""

from __future__ import annotations

import math

import numpy as np

MU0 = 4e-7 * math.pi

__all__ = [
    "loss_period_deg", "slot_period_deg", "slot_frequency", "skin_depth", "period_average",
    "shell_loss_from_potential", "shell_loss_from_br", "rotor_frame",
]


# ---------------------------------------------------------------------------------------------------------
# periods and resolution
# ---------------------------------------------------------------------------------------------------------
def loss_period_deg(slots: int, poles: int, phases: int = 3) -> float:
    """Period of the sleeve loss in electrical degrees (rotor frame, steady state).

    For an integer-slot winding (slots divisible by phases x poles) fed with balanced sinusoidal
    currents, advancing the rotor by one phase belt (360 / (2 phases) electrical degrees, q slots)
    maps the stator slots onto themselves and the current pattern onto itself (phase sequence A, -C,
    B, ...), so every field quantity seen from the rotor repeats after 60 electrical degrees
    (three phases). The loss, quadratic in the field, repeats with it (it may also repeat sooner).
    """
    if slots % (phases * poles):
        raise ValueError("loss_period_deg: only integer-slot windings (slots divisible by phases x poles)")
    return 360.0 / (2 * phases)


def slot_period_deg(slots: int, poles: int) -> float:
    """Slot pitch in electrical degrees: the period of slot passing seen from the rotor."""
    return 360.0 * (poles // 2) / slots


def slot_frequency(slots: int, rpm: float) -> float:
    """Slot-passing frequency seen from the rotor, Hz."""
    return slots * rpm / 60.0


def skin_depth(sigma: float, freq: float, mu_r: float = 1.0) -> float:
    """Skin depth sqrt(2 / (omega mu sigma)) in m (inf for sigma = 0 or freq = 0)."""
    if sigma <= 0 or freq <= 0:
        return math.inf
    return math.sqrt(2.0 / (2 * math.pi * freq * MU0 * mu_r * sigma))


# ---------------------------------------------------------------------------------------------------------
# averaging
# ---------------------------------------------------------------------------------------------------------
def _integral(t: np.ndarray, y: np.ndarray, a: float, b: float) -> float:
    """Integral of the piecewise-linear interpolant of (t, y) over [a, b] (inside [t0, t-1])."""
    inside = (t > a) & (t < b)
    tt = np.concatenate([[a], t[inside], [b]])
    yy = np.interp(tt, t, y)
    return float(np.trapezoid(yy, tt))


def period_average(t, y, period: float, n_periods: int = 1, tol: float = 1e-9) -> dict:
    """Time average of a periodic-in-steady-state signal over its last ``n_periods`` full periods,
    and the same average one block earlier as a convergence check.

    ``t`` increasing (s), ``y`` the signal, ``period`` (s). The averages are integrals of the
    piecewise-linear interpolant (trapezoid rule; for samples that fall on the window edges and a
    periodic signal this equals the mean of the samples of one period). Returns ``mean`` (last
    block), ``previous`` (the block before it, or nan when the record is too short), ``rel_change`` =
    |mean - previous| / |mean| (nan if unavailable) and the window ``(t_start, t_end)``.
    Raises ValueError when the record is shorter than one block.
    """
    t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
    if t.ndim != 1 or t.shape != y.shape or len(t) < 2:
        raise ValueError("period_average: t and y must be 1D arrays of the same length >= 2")
    if np.any(np.diff(t) <= 0):
        raise ValueError("period_average: t must be strictly increasing")
    if not period > 0 or n_periods < 1:
        raise ValueError("period_average: need period > 0 and n_periods >= 1")
    block = n_periods * period
    t1 = t[-1]
    t0 = t1 - block
    span = t1 - t[0]
    if t0 < t[0] - tol * max(span, block):
        raise ValueError(f"period_average: record of {span:.6g} s is shorter than {n_periods} period(s) "
                         f"of {period:.6g} s")
    t0 = max(t0, t[0])
    mean = _integral(t, y, t0, t1) / block
    prev, rel = math.nan, math.nan
    if t0 - block >= t[0] - tol * max(span, block):
        prev = _integral(t, y, max(t0 - block, t[0]), t0) / block
        rel = abs(mean - prev) / abs(mean) if mean else (0.0 if prev == 0 else math.inf)
    return dict(mean=mean, previous=prev, rel_change=rel, window=(t0, t1))


# ---------------------------------------------------------------------------------------------------------
# thin-shell estimate of the resistance-limited loss
# ---------------------------------------------------------------------------------------------------------
def rotor_frame(phi, t, field, omega_mech: float, theta):
    """Resample a field sampled at fixed (stator-frame) angles ``phi`` (rad, uniform over one full
    period of the field, endpoint excluded) and times ``t`` into rotor-frame angles ``theta``:
    out[k, j] = field[k](phi = theta[j] + omega_mech t[k]). ``field`` shape (n_t, n_phi); the field
    must be periodic in phi with period phi[-1] - phi[0] + dphi (pass a full ring or one
    period of an even-periodic field). Linear interpolation."""
    phi, t, theta = (np.asarray(x, dtype=float) for x in (phi, t, theta))
    f = np.asarray(field, dtype=float)
    per = phi[-1] - phi[0] + (phi[1] - phi[0])
    out = np.empty((len(t), len(theta)))
    for k in range(len(t)):
        out[k] = np.interp(theta + omega_mech * t[k], phi, f[k], period=per)
    return out


def _dt(a: np.ndarray, t: np.ndarray, method: str) -> np.ndarray:
    """Time derivative along axis 0 (backward: (a_k - a_k-1)/dt at t_k, the first row dropped)."""
    if method == "backward":
        return np.diff(a, axis=0) / np.diff(t)[(slice(None),) + (None,) * (a.ndim - 1)]
    if method == "central":
        return np.gradient(a, t, axis=0)
    raise ValueError("method must be 'backward' or 'central'")


def shell_loss_from_potential(t, theta, a, sigma: float, radii, length: float, net_current: str = "ring",
                              method: str = "backward", sectors: int = 1) -> dict:
    """Resistance-limited eddy loss of a conducting shell from the vector potential of the
    NON-conducting solution, seen in the rotor frame.

    Derivation. If the shell's own eddy currents do not change the field (resistance-limited: shell
    thickness << skin depth and small shell reaction, which the sigma-proportionality of the FE loss
    checks), the rotor-frame field is that of the non-conducting solution. With axial current only,
    E_z = -dA/dt + C(t), J_z = sigma E_z and the loss is P = int sigma E_z^2 dV. The uniform field C(t)
    is set by the end connection. Ideal end return: the axial current filaments are joined at both
    ends by resistance-free paths (the 2D idealisation: the circumferential end-region current of the
    real sleeve and its resistance are left out). There is no external circuit, so the total axial
    current through the cross-section is zero, int J_z dA = 0, which gives C = <dA/dt>, the mean over
    the conducting cross-section. Hence

        P = sigma L int_{r_in}^{r_out} int_0^{2 pi} (dA/dt - <dA/dt>)^2 r dtheta dr   (time-averaged).

    ``net_current``: 'ring' subtracts the mean over the whole ring at every time (zero total current
    in the full ring: the physical case), 'none' takes C = 0 (each filament shorted to an external
    ideal return), 'radius' subtracts the mean at each radius separately (zero current in each thin
    sub-shell), 'sector' subtracts the mean over each of ``sectors`` equal angular sectors starting at
    theta[0] (zero current in each sector, which is what a model of one sector that forces zero net
    current in its conductor would impose). For a field that is anti-periodic over a pole pitch the
    ring mean of dA/dt vanishes, so 'ring' and 'none' agree.

    Arguments: ``t`` (s, uniform, spanning an integer number of periods of the loss, endpoints
    included), ``theta`` (rad, uniform over the full ring, endpoint excluded), ``a`` [Wb/m] of shape
    (n_r, n_t, n_theta) at ``radii`` (m, increasing; the shell is radii[0] .. radii[-1], trapezoid
    rule through the thickness), ``sigma`` S/m, ``length`` m. ``method`` 'backward' takes the
    derivative as Maxwell's backward-Euler time stepping does ((A_k - A_k-1)/dt at t_k) and averages
    over the n_t - 1 intervals; 'central' uses central differences and the trapezoid rule.

    Returns ``loss`` (W, time average), ``loss_t`` (instantaneous, W) and ``t`` (its times).
    """
    t, theta, radii = (np.asarray(x, dtype=float) for x in (t, theta, radii))
    a = np.asarray(a, dtype=float)
    if a.shape != (len(radii), len(t), len(theta)):
        raise ValueError(f"a must have shape (n_r, n_t, n_theta) = {(len(radii), len(t), len(theta))}, got {a.shape}")
    if net_current not in ("ring", "none", "radius", "sector"):
        raise ValueError("net_current must be 'ring', 'none', 'radius' or 'sector'")
    e = -_dt(np.moveaxis(a, 1, 0), t, method)            # (n_t', n_r, n_theta)
    w = _trapz_weights(radii) * radii                    # r dr weights (trapezoid in r)
    if net_current == "radius":
        e = e - e.mean(axis=2, keepdims=True)
    elif net_current == "ring":
        cmean = np.einsum("krj,r->k", e, w) / (w.sum() * e.shape[2])
        e = e - cmean[:, None, None]
    elif net_current == "sector":
        if sectors < 1 or len(theta) % sectors:
            raise ValueError("net_current='sector' needs n_theta divisible by sectors >= 1")
        m = len(theta) // sectors
        es = e.reshape(e.shape[0], e.shape[1], sectors, m)
        cmean = np.einsum("krsj,r->ks", es, w) / (w.sum() * m)
        e = (es - cmean[:, None, :, None]).reshape(e.shape)
    dth = 2 * math.pi / len(theta)
    ring = (e ** 2).sum(axis=2) * dth                    # (n_t', n_r): int E^2 dtheta
    loss_t = sigma * length * (ring * w).sum(axis=1)
    if method == "backward":
        tt = t[1:]
        loss = float(loss_t.mean())
    else:
        tt = t
        loss = float(np.trapezoid(loss_t, t) / (t[-1] - t[0]))
    return dict(loss=loss, loss_t=loss_t, t=tt)


def _trapz_weights(x: np.ndarray) -> np.ndarray:
    if len(x) == 1:
        raise ValueError("need at least two radii through the shell thickness")
    w = np.zeros(len(x))
    d = np.diff(x)
    w[:-1] += d / 2
    w[1:] += d / 2
    return w


def shell_loss_from_br(t, theta, br, sigma: float, r_in: float, r_out: float, length: float) -> dict:
    """Harmonic form of the thin-shell estimate from the radial flux density B_r(theta, t) at the
    mid-radius r_m = (r_in + r_out)/2, rotor frame.

    With B_r = (1/r) dA/dtheta, a rotor-frame harmonic B_r = c_nk exp(i(n theta + w_k t)) has
    A_nk = r c_nk / (i n) and E_nk = -dA/dt = -w_k r c_nk / n, so with the field taken as uniform
    through the thickness d = r_out - r_in::

        P = sigma L d 2 pi r_m sum_{n != 0} sum_k (w_k r_m / n)^2 |c_nk|^2

    (complex coefficients of both signs, so a real travelling wave B cos(n theta - w t) gives
    sigma L d 2 pi r_m (w r_m B / n)^2 / 2). The n = 0 term is absent: B_r has no mean over a
    ring (no net flux) and the mean of E over the ring, which B_r cannot determine, is exactly what
    the zero-total-current condition removes. ``t`` uniform over an integer number of periods with
    the endpoint EXCLUDED, ``theta`` uniform over the full ring with the endpoint excluded, ``br``
    shape (n_t, n_theta), Tesla.

    Each harmonic travels relative to the shell at v = -w_k / n (rad/s, rotor frame; a harmonic at
    rest in the stator moves at -Omega). With J_z = sigma E_z and E_z = v r B_r for a travelling
    field, the tangential force density J_z B_r gives a torque on the shell of P_h / v per harmonic
    (cross terms vanish in the time and ring average), so the eddy currents add the torque
    ``torque`` = sum P_h / v_h (N m, positive in the direction of +theta, the direction of rotation;
    a stator-fixed harmonic gives the drag -P_h / Omega).

    Returns ``loss`` (W), ``torque`` (N m) and ``spectrum`` (n, f_Hz, contribution W) of the largest
    contributions.
    """
    t, theta = np.asarray(t, dtype=float), np.asarray(theta, dtype=float)
    br = np.asarray(br, dtype=float)
    nt, nth = br.shape
    if (nt, nth) != (len(t), len(theta)):
        raise ValueError("br must have shape (n_t, n_theta)")
    c = np.fft.fft2(br) / (nt * nth)                     # c[k, n] for exp(+i(w_k t + n theta)) (numpy sign)
    T = (t[1] - t[0]) * nt
    w = 2 * math.pi * np.fft.fftfreq(nt, d=T / nt)       # rad/s
    n = np.fft.fftfreq(nth, d=1.0 / nth)                  # harmonic order
    rm, d = 0.5 * (r_in + r_out), r_out - r_in
    W, N = np.meshgrid(w, n, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(N != 0, (W * rm / np.where(N != 0, N, 1.0)) ** 2 * np.abs(c) ** 2, 0.0)
    contrib = sigma * length * d * 2 * math.pi * rm * term
    with np.errstate(divide="ignore", invalid="ignore"):
        per_v = np.where(W != 0, -N / np.where(W != 0, W, 1.0), 0.0)       # 1 / v_h
    torque = float((contrib * per_v).sum())
    order = np.argsort(contrib.ravel())[::-1][:12]
    spec = [(int(abs(N.ravel()[i])), abs(W.ravel()[i]) / (2 * math.pi), float(contrib.ravel()[i])) for i in order]
    return dict(loss=float(contrib.sum()), torque=torque, spectrum=spec)
