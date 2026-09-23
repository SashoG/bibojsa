"""KTP (KTiOPO4, potassium titanyl phosphate) Sellmeier data.

Orthorhombic, biaxial, point group mm2; n_x < n_y < n_z (alpha, beta,
gamma axes respectively in Kato's notation). Sellmeier recipe:

    n_x      : Kato & Takaoka (2002) "double-resonance" fit::

                   n^2(lambda_um) = A + B*lambda_um^2/(lambda_um^2-C)
                                      + D/(lambda_um^2-E)

    n_y      : Konig & Wong (2004) "modified Fan" fit, Eq. (3) of that
               paper, obtained from their own extensive three-wave-mixing
               (difference-frequency) measurements around 791/1582 nm
               specifically because "there are no published data on n_y
               near 1600 nm" in the equations they surveyed::

                   n_y^2 = A1 + A2/(1 - A3/lambda_um^2) - A4*lambda_um^2
                         = A1 + A2*lambda_um^2/(lambda_um^2-A3) - A4*lambda_um^2

               (the two forms are algebraically identical; the second is
               used here since it matches this module's
               `TwoPoleLinearSellmeierCoeffs` with B2=0.)

    n_z      : Fradkin et al. (1999) "two-pole + linear" fit. Konig & Wong
               explicitly adopted Fradkin's n_z unmodified ("We fit the
               data in Fig. 1 to a set of Sellmeier equations by using
               Fradkin's equation for n_z and a modified Fan's equation
               for n_y") -- i.e. this project's pre-existing n_z choice
               (from milestone 2's earlier iteration) was already exactly
               what Konig & Wong used, and required no change once their
               paper was located.

MILESTONE-2 UPDATE (superseding the previous hybrid recipe): the
project's earlier n_y was Kato & Takaoka's own fit (chosen as the closest
available proxy for Konig & Wong's, which was believed unsourceable --
paywalled, not reproduced verbatim in any secondary reference found by
search). That proxy reproduced n_y's *absolute* value well but not its
*dispersion slope* near 791 nm closely enough, which shifted the computed
GVM regime for the milestone-2 design point away from the literature's
symmetric/EPM classification (see tests/test_milestone2_ktp.py and
PROGRESS.md for the full quantified discussion of that gap and its
resolution). Now that Konig & Wong, Appl. Phys. Lett. 84, 1644 (2004) has
been located directly (papers/09.), its own Eq. (3) coefficients are used
verbatim instead: A1=2.09930, A2=0.922683, A3=0.0467695, A4=0.0138408
(lambda in um). This closes the milestone-2 gap: GD_p, GD_s, GD_i, D, and
the QPM period all now match Laudenbach et al. (2017) to well under 1%.

Validity range for n_y: restricted to 0.75-1.65 um, matching the range
Konig & Wong actually measured and fit (pump 756-832 nm, probe/DF
1520-1630 nm) -- NOT extrapolated to BIBO/BIBO-like full-transparency
ranges the way n_x's Kato-Takaoka fit nominally is, since Konig & Wong's
own fit was constructed and validated only in this window.

References:
    K. Kato and E. Takaoka, "Sellmeier and thermo-optic dispersion
    formulas for KTP," Appl. Opt. 41, 5040-5044 (2002). Coefficients as
    transcribed by refractiveindex.info (Polyanskiy, Sci. Data 11, 94
    (2024)):
        n_alpha^2 (=nx) = 3.29100 + 0.04140/(lam^2-0.03978) + 9.35522/(lam^2-31.45571)

    F. Konig and F. N. C. Wong, "Extended phase matching of second-
    harmonic generation in periodically poled KTiOPO4 with zero
    group-velocity mismatch," Appl. Phys. Lett. 84, 1644 (2004)
    [papers/09. in this project] -- source of n_y (Eq. 3) used here, and
    of the precise zero-GVM/EPM design point (grating period
    Lambda_e = 46.15 um, pump lambda_e = 791.1-792.6 nm, degenerate
    2*lambda_e ~ 1582-1584 nm, FWHM extended-phase-matching bandwidth
    67-100 nm).

    K. Fradkin, A. Arie, A. Skliar, and G. Rosenman, "Tunable midinfrared
    source by difference frequency generation in bulk periodically poled
    KTiOPO4," Appl. Phys. Lett. 74, 914-916 (1999); coefficients as
    reproduced in D. A. Heberle and Z. H. Levine, J. Res. NIST 118, 375
    (2013), and adopted directly (unmodified) by Konig & Wong (2004) for
    n_z:
        n_z^2 = 2.12725 + 1.18431*lam^2/(lam^2-0.0514852)
                         + 0.6603*lam^2/(lam^2-100.00507)
                         - 0.00968956*lam^2

    F. Laudenbach, R.-B. Jin, C. Greganti, M. Hentschel, P. Walther, and
    H. Hubel, "Numerical Investigation of Photon-Pair Generation in
    Periodically Poled MTiOXO4 (M=K,Rb,Cs; X=P,As)," Phys. Rev. Applied
    8, 024035 (2017) [arXiv:1706.04559] -- source of the precise
    milestone-2 validation numbers (group delays, dispersion parameter D,
    purity, QPM period) for 791 nm(o) -> 1582 nm(o)+1582 nm(e).

Used for validation milestone 2 (this project): type-II periodically-poled
KTP (PPKTP), collinear propagation along the crystallographic x-axis
(walk-off-free principal axis), degenerate near 1582-1584 nm from a
~791-793 nm pump, quasi-phase-matched with grating period Lambda ~ 46.1-46.2 um.
"""
from ..sellmeier import ZernikeSellmeierCoeffs, TwoPoleLinearSellmeierCoeffs

_LMIN, _LMAX = 0.43, 3.54
_SOURCE_X = "Kato & Takaoka, Appl. Opt. 41, 5040 (2002), via refractiveindex.info"
_SOURCE_Y = "Konig & Wong, Appl. Phys. Lett. 84, 1644 (2004), Eq. (3) [papers/09.]"
_SOURCE_Z = "Fradkin et al., Appl. Phys. Lett. 74, 914 (1999), adopted unmodified by Konig & Wong (2004)"

KTP_X = ZernikeSellmeierCoeffs(
    A=3.29100, B=0.04140, C=0.03978, D=9.35522, E=31.45571,
    lambda_min=_LMIN, lambda_max=_LMAX, axis="x", source=_SOURCE_X,
)

# Kato & Takaoka's own n_y, retained only for reference/comparison -- NOT
# used in the KTP dict below (superseded by Konig & Wong's own n_y, see
# module docstring: their fit's dispersion *slope* near 791 nm differed
# enough from ours to misclassify the GVM regime at the milestone-2
# design point).
KTP_Y_KATO_TAKAOKA = ZernikeSellmeierCoeffs(
    A=3.45018, B=0.04341, C=0.04597, D=16.98825, E=39.43799,
    lambda_min=_LMIN, lambda_max=_LMAX, axis="y", source=_SOURCE_X,
)
KTP_Y = TwoPoleLinearSellmeierCoeffs(
    A=2.09930, B1=0.922683, C1=0.0467695, B2=0.0, C2=1.0e6, D=0.0138408,
    lambda_min=0.75, lambda_max=1.65, axis="y", source=_SOURCE_Y,
)
KTP_Z = TwoPoleLinearSellmeierCoeffs(
    A=2.12725, B1=1.18431, C1=0.0514852, B2=0.6603, C2=100.00507, D=0.00968956,
    lambda_min=_LMIN, lambda_max=_LMAX, axis="z", source=_SOURCE_Z,
)

KTP = {"x": KTP_X, "y": KTP_Y, "z": KTP_Z}
