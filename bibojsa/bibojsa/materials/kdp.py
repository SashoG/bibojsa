"""KDP (KH2PO4, potassium dihydrogen phosphate) Sellmeier data.

Tetragonal, negative uniaxial (point group -42m); optic axis along the
crystallographic c-axis. Zernike (1964) dispersion fit, "double-resonance"
Sellmeier form::

    n^2(lambda_um) = A + B*lambda_um^2/(lambda_um^2 - C) + D/(lambda_um^2 - E)

valid 0.2138-1.529 um at 24.8 degC.

Reference: F. Zernike, Jr., "Refractive indices of ammonium dihydrogen
phosphate and potassium dihydrogen phosphate between 2000 A and 1.5 um,"
J. Opt. Soc. Am. 54, 1215-1219 (1964). Coefficients as transcribed by
refractiveindex.info (Polyanskiy, Sci. Data 11, 94 (2024)):
    n_o^2 = 2.259276 + 13.00522 lam^2/(lam^2-400) + 0.01008956/(lam^2-0.0129426)
    n_e^2 = 2.132668 + 3.2279924 lam^2/(lam^2-400) + 0.008637494/(lam^2-0.0122810)

Used for validation milestone 1 (this project): type-II collinear SPDC in
KDP pumped near 415 nm producing degenerate 830 nm photon pairs with
group-velocity matching between the pump (e-polarized) and the o-polarized
daughter photon -- Mosley et al., Phys. Rev. Lett. 100, 133601 (2008).
"""
from ..sellmeier import ZernikeSellmeierCoeffs

_LMIN, _LMAX = 0.2138, 1.529
_SOURCE = "Zernike, J. Opt. Soc. Am. 54, 1215 (1964), via refractiveindex.info"

KDP_O = ZernikeSellmeierCoeffs(
    A=2.259276, B=13.00522, C=400.0, D=0.01008956, E=0.0129426,
    lambda_min=_LMIN, lambda_max=_LMAX, axis="o", source=_SOURCE,
)
KDP_E = ZernikeSellmeierCoeffs(
    A=2.132668, B=3.2279924, C=400.0, D=0.008637494, E=0.0122810,
    lambda_min=_LMIN, lambda_max=_LMAX, axis="e", source=_SOURCE,
)

# Uniform (x, y, z) interface for the biaxial/principal-plane machinery:
# uniaxial crystal, optic axis along z (crystallographic c-axis), so
# nx = ny = n_o and nz = n_e. Negative uniaxial: n_e < n_o.
KDP = {"x": KDP_O, "y": KDP_O, "z": KDP_E, "o": KDP_O, "e": KDP_E}
