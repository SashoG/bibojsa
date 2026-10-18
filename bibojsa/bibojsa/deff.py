"""Effective nonlinearity deff for BIBO in the principal planes.

Implements LPR2010 Table 1 ("Effective nonlinearity deff of BIBO in the
principal optical planes (Kleinman symmetry assumed)") verbatim, using
the Kleinman-averaged coefficients from LPR2010 Sec. 2 / [73]:

    d11 = 2.53, d12 = d26 = 3.2, d13 = d35 = -1.76,
    d14 = d25 = d36 = 1.66   [pm/V, xyz dielectric frame]

Conventions (LPR2010 Sec. 2-3, followed exactly):

- theta, phi are the standard polar/azimuthal spherical angles in the
  dielectric xyz frame (theta from z, phi from x).
- The o/e labels are per principal plane (o = polarized perpendicular to
  the plane, e = polarized in the plane), used in the sequence
  lambda1 lambda2 lambda3 with lambda1 >= lambda2 > lambda3, i.e.
  (idler, signal, pump) in this project's SPDC naming.
- BIBO is monoclinic, point group 2, with the two-fold axis along the
  dielectric x-axis (this is why exactly the eight coefficients d11, d12,
  d13, d14, d25, d26, d35, d36 survive: all have an odd number of x
  indices). Consequently deff in the y-z plane (the plane normal to the
  two-fold axis) is NEITHER even NOR odd in theta: the two octants are
  inequivalent, and expressions must not be symmetrized (no orthorhombic
  shortcuts -- see the project brief's caution and LPR2010 Sec. 2-3).
- Octant convention (LPR2010 Sec. 3): first octant (0 <= phi <= 90 deg,
  0 <= theta <= 90 deg) plus the adjacent second octant
  (0 <= phi <= 90 deg, 90 deg <= theta <= 180 deg). Table 1's phi = 180
  and 270 deg rows duplicate these by 2/m symmetry and are not
  implemented ("all cases with phi = 180 deg and 270 deg in Table 1 will
  be irrelevant because they simply duplicate others").

Table 1 rows implemented (phi = 0 deg for x-z, phi = +90 deg for y-z;
Omega = optic-axis half-angle from z, wavelength-dependent):

    x-y plane, theta = 90 deg, angle = phi:
        ooe   (type I-)           : -d13 sin(phi)
        eoe = oee (II- = III-)    :  d14 sin(2 phi)

    y-z plane, phi = +90 deg, angle = theta (0..180 deg; octants
    inequivalent -- the sin(2 theta) term changes sign about 90 deg while
    the cos^2/sin^2 terms do not):
        eeo   (type I+)           : -d12 cos^2(theta) - d13 sin^2(theta)
                                     + d14 sin(2 theta)
        oeo = eoo (II+ = III+)    :  0   (exactly; the surviving d-tensor
                                     components cannot couple two x-polarized
                                     o-waves to one in-plane e-wave)

    x-z plane, phi = 0 deg, angle = theta (0..180 deg):
        eeo   (I+),  Omega < theta < 180 deg - Omega : -d14 sin(2 theta)
        oeo = eoo (II+ = III+), same range           :  d12 cos(theta)
        ooe   (I-),  theta < Omega or theta > 180 deg - Omega
                                                     : -d12 cos(theta)
        eoe = oee (II- = III-), same range           : -d14 sin(2 theta)

    (The x-z angle-range annotations mark where each interaction type is
    phase-matchable -- inside vs. outside the optic axes at Omega and
    180 deg - Omega; the algebraic expressions themselves are valid at any
    theta in the plane, so `enforce_range=False` allows evaluating them
    outside the nominal region, e.g. for plotting deff across the whole
    plane.)

Only |deff| is physically significant for the SPDC conversion-efficiency
and brightness-proxy purposes of this project (deff enters as deff^2);
the signed values above follow Table 1's sign conventions so they can be
compared against the paper directly.

Walk-off corrections to deff are neglected (project brief; LPR2010 does
the same), as is the <= +-1 deg wavelength-dependent rotation of the xyz
frame about the crystallographic b-axis.

Reference: V. Petrov, M. Ghotbi, O. Kokabee, A. Esteban-Martin, F. Noack,
M. Ebrahim-Zadeh, "Femtosecond nonlinear frequency conversion based on
BiB3O6," Laser & Photonics Reviews 4, 53-98 (2010) [LPR2010], Table 1 and
Sec. 2-3.
"""
from __future__ import annotations

import numpy as np

from .materials.bibo import BIBO, D_COEFFS_PM_PER_V
from .sellmeier import optic_axis_angle_deg

# Interaction-type equivalences per LPR2010 Table 1 (lambda1 lambda2 lambda3
# = idler, signal, pump sequence): the two type-II/III labellings that
# differ only by which daughter is 'e' share one deff expression.
_EQUIVALENT = {
    "oee": "eoe",
    "eoo": "oeo",
}

# Roman-numeral type labels (LPR2010 Table 1), for reporting convenience.
TYPE_LABELS = {
    ("xy", "ooe"): "I-",
    ("xy", "eoe"): "II-/III-",
    ("yz", "eeo"): "I+",
    ("yz", "oeo"): "II+/III+",
    ("xz", "eeo"): "I+",
    ("xz", "oeo"): "II+/III+",
    ("xz", "ooe"): "I-",
    ("xz", "eoe"): "II-/III-",
    # partner sequences (e/o roles of the two daughters swapped; same deff)
    ("xy", "oee"): "II-/III-",
    ("yz", "eoo"): "II+/III+",
    ("xz", "eoo"): "II+/III+",
    ("xz", "oee"): "II-/III-",
}


def normalize_oe_sequence(pol_idler: str, pol_signal: str, pol_pump: str) -> str:
    """Build the LPR2010 lambda1-lambda2-lambda3 o/e sequence (idler,
    signal, pump) and map it onto Table 1's canonical representative
    (oee -> eoe, eoo -> oeo)."""
    seq = f"{pol_idler}{pol_signal}{pol_pump}"
    if any(c not in ("o", "e") for c in seq):
        raise ValueError(f"polarizations must be 'o'/'e', got {seq!r}")
    return _EQUIVALENT.get(seq, seq)


def deff_bibo_pm_per_v(plane: str, pol_pump: str, pol_signal: str, pol_idler: str,
                        angle_rad: float, lambda_pump_um: float | None = None,
                        phi_sign: int = +1, enforce_range: bool = True,
                        d: dict[str, float] = D_COEFFS_PM_PER_V) -> float:
    """Signed effective nonlinearity deff [pm/V] for BIBO, LPR2010 Table 1.

    Parameters
    ----------
    plane : 'xy', 'yz', or 'xz' principal plane.
    pol_pump, pol_signal, pol_idler : 'o'/'e' per the principal-plane
        convention (same convention as `bibojsa.phasematch.ThreeWaveConfig`).
    angle_rad : the propagation angle *in radians*: phi (from x) for the
        x-y plane (theta = 90 deg implied); theta (from z) for the y-z
        plane (phi = +90 deg implied, see `phi_sign`) and the x-z plane
        (phi = 0 deg implied). For y-z and x-z, theta may run over
        [0, 180 deg] -- both (inequivalent, monoclinic) octants.
    lambda_pump_um : pump vacuum wavelength (um). Only used for the x-z
        plane's optic-axis-range check (Omega is wavelength dependent);
        required there when `enforce_range=True`, ignored otherwise.
    phi_sign : +1 or -1, selects phi = +90 deg or phi = -90 deg for the
        y-z plane (LPR2010 lists 'phi = +-90 deg' with upper/lower signs
        on the d12/d13 terms: deff = -+ d12 cos^2 -+ d13 sin^2 + d14 sin2theta).
        Ignored for the other planes. Default +1 (first/second octant,
        consistent with the project's octant convention).
    enforce_range : if True (default), raise ValueError when the requested
        x-z interaction type is evaluated outside its phase-matchable
        theta range relative to the optic axes (inside Omega..180-Omega
        for eeo/oeo; outside it for ooe/eoe), per Table 1's annotations.
        The y-z and x-y rows carry no such range restriction.
    d : d-coefficient dict (pm/V); defaults to LPR2010's Kleinman set.

    Returns the *signed* deff in pm/V (only |deff| matters for the SPDC
    efficiency/brightness figures in this project).
    """
    if plane not in ("xy", "yz", "xz"):
        raise ValueError(f"plane must be 'xy', 'yz' or 'xz', got {plane!r}")
    if phi_sign not in (+1, -1):
        raise ValueError(f"phi_sign must be +1 or -1, got {phi_sign!r}")

    seq = normalize_oe_sequence(pol_idler, pol_signal, pol_pump)
    a = float(angle_rad)
    d11, d12, d13, d14 = d["d11"], d["d12"], d["d13"], d["d14"]

    if plane == "xy":
        # theta = 90 deg; angle = phi
        if seq == "ooe":
            return -d13 * np.sin(a)
        if seq == "eoe":
            return d14 * np.sin(2.0 * a)
        raise ValueError(
            f"interaction {seq!r} is not a phase-matchable LPR2010 Table 1 "
            f"entry in the x-y plane (available: ooe [I-], eoe/oee [II-/III-])."
        )

    if plane == "yz":
        # phi = +-90 deg; angle = theta. Monoclinic: neither even nor odd
        # in theta -- do NOT symmetrize.
        s = float(phi_sign)
        if seq == "eeo":
            return (-s * d12 * np.cos(a) ** 2
                    - s * d13 * np.sin(a) ** 2
                    + d14 * np.sin(2.0 * a))
        if seq == "oeo":
            # Exactly zero: two o-waves are both x-polarized, and every
            # surviving point-group-2 coefficient coupling (x, x, in-plane)
            # (d15 = d16 = 0) vanishes.
            return 0.0
        raise ValueError(
            f"interaction {seq!r} is not an LPR2010 Table 1 entry in the "
            f"y-z plane (available: eeo [I+], oeo/eoo [II+/III+, deff=0])."
        )

    # plane == 'xz'; phi = 0 deg; angle = theta in [0, 180 deg]
    if seq in ("eeo", "oeo"):
        inside_required = True
    elif seq in ("ooe", "eoe"):
        inside_required = False
    else:
        raise ValueError(
            f"interaction {seq!r} is not an LPR2010 Table 1 entry in the "
            f"x-z plane (available: eeo/oeo [+ types, between the optic "
            f"axes], ooe/eoe [- types, outside them])."
        )

    if enforce_range:
        if lambda_pump_um is None:
            raise ValueError(
                "lambda_pump_um is required for the x-z plane range check "
                "(the optic-axis angle Omega is wavelength dependent); pass "
                "enforce_range=False to skip it."
            )
        omega_deg = optic_axis_angle_deg(lambda_pump_um, BIBO["x"], BIBO["y"], BIBO["z"])
        theta_deg = np.degrees(a)
        inside = omega_deg < theta_deg < 180.0 - omega_deg
        if inside != inside_required:
            region = ("Omega < theta < 180deg - Omega" if inside_required
                      else "theta < Omega or theta > 180deg - Omega")
            raise ValueError(
                f"interaction {seq!r} in the x-z plane is phase-matchable for "
                f"{region} (Omega({lambda_pump_um} um) = {omega_deg:.2f} deg), "
                f"but theta = {theta_deg:.2f} deg was requested. Pass "
                f"enforce_range=False to evaluate the bare Table 1 expression "
                f"anyway (e.g. for plotting)."
            )

    if seq == "eeo":
        return -d14 * np.sin(2.0 * a)
    if seq == "oeo":
        return d12 * np.cos(a)
    if seq == "ooe":
        return -d12 * np.cos(a)
    # seq == 'eoe'
    return -d14 * np.sin(2.0 * a)


def deff_bibo_from_config(config, theta_rad: float, lambda_pump_um: float | None = None,
                           phi_sign: int = +1, enforce_range: bool = True) -> float:
    """Convenience wrapper taking a `bibojsa.phasematch.ThreeWaveConfig`
    (which stores pump/signal/idler polarizations and the principal plane)
    plus the crystal-cut angle, and returning the signed deff [pm/V]."""
    return deff_bibo_pm_per_v(config.plane, config.pol_pump, config.pol_signal,
                               config.pol_idler, theta_rad,
                               lambda_pump_um=lambda_pump_um, phi_sign=phi_sign,
                               enforce_range=enforce_range)
