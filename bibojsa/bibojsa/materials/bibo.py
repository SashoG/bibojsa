"""BIBO (BiB3O6, bismuth triborate) Sellmeier data.

Miyata-Umemura-Kato fit as reproduced in LPR2010 Table 2, valid
0.32-3.09 um. Dielectric-frame convention n_x < n_y < n_z (monoclinic
point group 2; the xyz frame rotates slightly (<= +-1 deg) with
wavelength about the crystallographic b-axis -- neglected here, as in
LPR2010, and noted as an approximation, see bibojsa.biaxial docstring).

n^2(lambda) = A + B/(lambda^2 - C) - D*lambda^2, lambda in vacuum um.

Reference: V. Petrov et al., "Application of the nonlinear crystal
BiB3O6...", Laser & Photonics Reviews 4, 53 (2010) [LPR2010], Table 2.
Validation numbers (Table 3): n_x=1.77668, n_y=1.80641, n_z=1.94582 at
632.8 nm; n_x=1.75752, n_y=1.78400, n_z=1.91711 at 1064 nm; optic axis
2*Omega = 52.6 deg (matches 632.8 nm computation, see tests/ and note in
test_sellmeier_bibo.py -- the brief's "at 1064 nm" label for this number
looks like a transcription slip; flagged for confirmation against the
original LPR2010 table).
"""
from ..sellmeier import SellmeierCoeffs

_LMIN, _LMAX = 0.32, 3.09

BIBO = {
    "x": SellmeierCoeffs(A=3.0759, B=0.03169, C=0.03323, D=0.01402,
                          lambda_min=_LMIN, lambda_max=_LMAX, axis="x",
                          source="LPR2010 Table 2 (Miyata-Umemura-Kato fit)"),
    "y": SellmeierCoeffs(A=3.1698, B=0.03666, C=0.03599, D=0.01819,
                          lambda_min=_LMIN, lambda_max=_LMAX, axis="y",
                          source="LPR2010 Table 2 (Miyata-Umemura-Kato fit)"),
    "z": SellmeierCoeffs(A=3.6546, B=0.05116, C=0.03713, D=0.02299,
                          lambda_min=_LMIN, lambda_max=_LMAX, axis="z",
                          source="LPR2010 Table 2 (Miyata-Umemura-Kato fit)"),
}

# Kleinman-symmetry d-tensor coefficients, xyz dielectric frame, pm/V
# (LPR2010 Table 1 and associated text).
D_COEFFS_PM_PER_V = {
    "d11": 2.53,
    "d12": 3.2,   # = d26
    "d26": 3.2,
    "d13": -1.76,  # = d35
    "d35": -1.76,
    "d14": 1.66,   # = d25 = d36
    "d25": 1.66,
    "d36": 1.66,
}


# --------------------------------------------------------------------
# Alternative BIBO dispersion fits, used only for the sensitivity check
# in the manuscript (Sec. 2.5): the GVM conditions depend on Sellmeier
# *slopes*, which are far less well constrained than absolute indices,
# so the crossing wavelengths must be quoted with the spread over the
# published fits and not to the precision of a single one.
# --------------------------------------------------------------------

# H. Hellwig, J. Liebertz, L. Bohaty, J. Appl. Phys. 88, 240 (2000);
# direct index measurements, 0.35-2.4 um.
BIBO_HELLWIG = {
    "x": SellmeierCoeffs(A=3.0740, B=0.0323, C=0.0316, D=0.01337,
                         lambda_min=0.35, lambda_max=2.40, axis="x",
                         source="Hellwig, Liebertz & Bohaty, JAP 88, 240 (2000)"),
    "y": SellmeierCoeffs(A=3.1685, B=0.0373, C=0.0346, D=0.01750,
                         lambda_min=0.35, lambda_max=2.40, axis="y",
                         source="Hellwig, Liebertz & Bohaty, JAP 88, 240 (2000)"),
    "z": SellmeierCoeffs(A=3.6545, B=0.0511, C=0.0371, D=0.0226,
                         lambda_min=0.35, lambda_max=2.40, axis="z",
                         source="Hellwig, Liebertz & Bohaty, JAP 88, 240 (2000)"),
}

# N. Umemura, K. Miyata, K. Kato, Opt. Mater. 30, 532 (2007); OPO-based
# refinement, 0.48-3.1 um.  Only n_y is materially changed from Hellwig.
BIBO_UMEMURA = {
    "x": SellmeierCoeffs(A=3.07403, B=0.03231, C=0.03163, D=0.013376,
                         lambda_min=0.48, lambda_max=3.10, axis="x",
                         source="Umemura, Miyata & Kato, Opt. Mater. 30, 532 (2007)"),
    "y": SellmeierCoeffs(A=3.16940, B=0.03717, C=0.03483, D=0.01827,
                         lambda_min=0.48, lambda_max=3.10, axis="y",
                         source="Umemura, Miyata & Kato, Opt. Mater. 30, 532 (2007)"),
    "z": SellmeierCoeffs(A=3.6545, B=0.05112, C=0.03713, D=0.02261,
                         lambda_min=0.48, lambda_max=3.10, axis="z",
                         source="Umemura, Miyata & Kato, Opt. Mater. 30, 532 (2007)"),
}

ALTERNATIVE_FITS = {"hellwig2000": BIBO_HELLWIG, "umemura2007": BIBO_UMEMURA}
