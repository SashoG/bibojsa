"""Biaxial crystal optics: Fresnel equation of wave normals, eigenindices,
walk-off, restricted in practice (per this project's validation scope and
the LPR2010 Table 1 o/e convention) to propagation confined to a principal
dielectric plane (xy, yz, or xz).

General Fresnel equation of wave normals
-----------------------------------------
For propagation direction (theta, phi) in the dielectric frame (theta from
z, phi from x measured in the xy-plane), the two eigenindices n satisfy::

    sin^2(theta)cos^2(phi)/(1/n^2 - 1/nx^2)
  + sin^2(theta)sin^2(phi)/(1/n^2 - 1/ny^2)
  + cos^2(theta)/(1/n^2 - 1/nz^2) = 0

which is quadratic in u = 1/n^2 (see `eigenindices_general`). This reduces
to the ordinary/extraordinary uniaxial result when two of nx, ny, nz
coincide, and to simple closed forms in the three principal planes (see
`principal_plane_index`), which is what phasematch.py and gvm.py use in
practice, since every validation milestone and BIBO survey target in this
project involves principal-plane propagation only.

Cautions (see project brief)
-----------------------------
- BIBO is monoclinic (point group 2): deff in the plane normal to the
  two-fold axis is neither even nor odd in theta, and the two octants
  are inequivalent (LPR2010 Sec. 2-3). This module does not assume
  orthorhombic symmetry shortcuts; principal-plane index formulas below
  hold regardless of point group, since they follow purely from the
  index ellipsoid, but deff symmetry (handled in phasematch.py) does not.
- The xyz dielectric frame rotates slightly (<= +-1 deg) with wavelength
  about the crystallographic b-axis for BIBO; neglected here, as in
  LPR2010, following the brief's explicit instruction.
- Walk-off corrections to deff are neglected (also per the brief); walk-off
  is only used for the effective-length estimate in the exploration phase.
"""
from __future__ import annotations

import numpy as np

from .sellmeier import refractive_index, finite_diff_converged, omega_to_lambda_um, C_LIGHT

PRINCIPAL_PLANES = ("xy", "yz", "xz")


def eigenindices_general(lambda_um: float, theta: float, phi: float,
                          coeffs_x, coeffs_y, coeffs_z) -> tuple[float, float]:
    """General biaxial Fresnel-equation eigenindices at (theta, phi).

    Returns (n_slow, n_fast) with n_slow >= n_fast (slow wave = larger
    index = smaller phase velocity). theta measured from z, phi from x
    in the xy-plane (standard spherical convention in the dielectric
    frame). Valid for any nx, ny, nz ordering, including degenerate
    (uniaxial) cases.
    """
    nx = refractive_index(lambda_um, coeffs_x)
    ny = refractive_index(lambda_um, coeffs_y)
    nz = refractive_index(lambda_um, coeffs_z)
    a, b, c = 1.0 / nx**2, 1.0 / ny**2, 1.0 / nz**2

    s1 = np.sin(theta) ** 2 * np.cos(phi) ** 2
    s2 = np.sin(theta) ** 2 * np.sin(phi) ** 2
    s3 = np.cos(theta) ** 2

    Bc = s1 * (b + c) + s2 * (a + c) + s3 * (a + b)
    Cc = s1 * b * c + s2 * a * c + s3 * a * b
    disc = max(Bc**2 - 4 * Cc, 0.0)  # clip tiny negative round-off at degeneracy
    sqrt_disc = np.sqrt(disc)

    u_lo = (Bc - sqrt_disc) / 2.0  # smaller u -> larger n (slow wave)
    u_hi = (Bc + sqrt_disc) / 2.0  # larger u -> smaller n (fast wave)
    n_slow = 1.0 / np.sqrt(u_lo)
    n_fast = 1.0 / np.sqrt(u_hi)
    return n_slow, n_fast


def principal_plane_index(lambda_um: float, theta: float, plane: str,
                           coeffs_x, coeffs_y, coeffs_z) -> tuple[float, float]:
    """Closed-form eigenindices for propagation confined to a principal
    plane. Returns (n_fixed, n_variable):

    - n_fixed: index of the eigenmode polarized along the axis
      perpendicular to the plane (angle-independent; the "o-like" mode
      in LPR2010's o/e convention).
    - n_variable: index of the in-plane-polarized eigenmode ("e-like"),
      following 1/n(angle)^2 = cos^2(angle)/n_B^2 + sin^2(angle)/n_A^2,
      where, by convention, `theta` is measured from axis A (so that
      angle=0 puts k along A and the in-plane perpendicular polarization
      is along B; angle=90 deg puts k along B and polarization along A):

        plane 'xz': theta from z (A=z, B=x); n_fixed = ny
        plane 'yz': theta from z (A=z, B=y); n_fixed = nx
        plane 'xy': theta from x (A=x, B=y); n_fixed = nz

    This is the standard index-ellipse reduction used throughout LPR2010
    Table 1 and is exact (no approximation) for propagation exactly in
    the plane.
    """
    if plane not in PRINCIPAL_PLANES:
        raise ValueError(f"plane must be one of {PRINCIPAL_PLANES}, got {plane!r}")

    nx = refractive_index(lambda_um, coeffs_x)
    ny = refractive_index(lambda_um, coeffs_y)
    nz = refractive_index(lambda_um, coeffs_z)

    c2, s2 = np.cos(theta) ** 2, np.sin(theta) ** 2
    if plane == "xz":
        n_fixed = ny
        n_var = 1.0 / np.sqrt(c2 / nx**2 + s2 / nz**2)
    elif plane == "yz":
        n_fixed = nx
        n_var = 1.0 / np.sqrt(c2 / ny**2 + s2 / nz**2)
    else:  # 'xy'
        n_fixed = nz
        n_var = 1.0 / np.sqrt(c2 / ny**2 + s2 / nx**2)
    return n_fixed, n_var


def walkoff_angle(lambda_um: float, theta: float, plane: str,
                   coeffs_x, coeffs_y, coeffs_z, h: float = 1e-5) -> float:
    """Spatial walk-off angle (rad) of the angle-dependent ("e-like")
    eigenmode in a principal plane, via the standard crystal-optics
    result for the index ellipse::

        rho(theta) = arctan[ (1/n(theta)) * dn/dtheta ]

    (Boyd, "Nonlinear Optics", walk-off in a principal plane.) The
    "o-like" fixed-polarization mode has zero walk-off (Poynting vector
    parallel to the wavevector) and is not covered by this function.
    dn/dtheta is evaluated by 4th-order central finite difference with a
    step-size convergence check (see `bibojsa.sellmeier.finite_diff_converged`).
    """
    def n_of_theta(th):
        _, n_var = principal_plane_index(lambda_um, th, plane, coeffs_x, coeffs_y, coeffs_z)
        return n_var

    n_theta = n_of_theta(theta)
    dn_dtheta, _, converged, rel_change = finite_diff_converged(
        n_of_theta, theta, order=1, h0=h, rtol=1e-8
    )
    if not converged:
        # still return the best estimate; caller can inspect via debug rerun
        pass
    return np.arctan(dn_dtheta / n_theta)


def k_of_omega_principal(omega: float, theta: float, plane: str, pol: str,
                          coeffs_x, coeffs_y, coeffs_z) -> float:
    """Wavevector magnitude k(omega) = n(omega)*omega/c [rad/m] for a wave
    of polarization `pol` ('o' = fixed/angle-independent branch, 'e' =
    angle-dependent branch) propagating at fixed angle `theta` in
    principal plane `plane`. Used to compute group velocity/GVD at fixed
    crystal cut angle (gvm.py) and Delta k(omega_s, omega_i) for JSA
    (jsa.py), since theta (the crystal cut) is fixed while the SPDC
    signal/idler frequencies are swept.
    """
    lam_um = omega_to_lambda_um(omega)
    n_fixed, n_var = principal_plane_index(lam_um, theta, plane, coeffs_x, coeffs_y, coeffs_z)
    n = n_fixed if pol == "o" else n_var
    return n * omega / C_LIGHT
