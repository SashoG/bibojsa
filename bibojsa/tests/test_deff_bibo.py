"""Validate the BIBO deff implementation (LPR2010 Table 1, Kleinman).

Three independent lines of validation:

1. Analytic spot checks of the Table 1 expressions at hand-computable
   angles (values worked out directly from the printed expressions).

2. Independent first-principles cross-check: contract the full point-
   group-2 d-tensor (two-fold axis || dielectric x; surviving elements
   d11, d12, d13, d14, d25=d14, d26=d12, d35=d13, d36=d14 under Kleinman)
   with explicit unit polarization vectors for each principal-plane o/e
   combination, and compare |deff| against the Table 1 expressions over a
   dense angle grid. This catches any transcription error in the table
   entries, since the contraction uses only the raw d_il values and
   geometry (sign conventions for the e/o unit vectors are fixed only up
   to +-1 per wave, so magnitudes are compared).

3. Literature benchmark: the well-known BIBO type-I SHG configuration in
   the y-z plane (1064 nm -> 532 nm at theta ~ 168.9 deg, i.e. the
   *second* octant) has |deff| ~ 3.6-3.7 pm/V, roughly 1.5x the value in
   the first octant at the mirror angle (~2.4 pm/V) -- the monoclinic
   octant inequivalence LPR2010 Sec. 2-3 emphasizes and the project brief
   explicitly warns about. Also the milestone-3 SPDC design point (x-z,
   ooe, theta ~ 11 deg) should give |deff| = d12 cos(theta) ~ 3.14 pm/V.
"""
import numpy as np
import pytest

from bibojsa.deff import deff_bibo_pm_per_v, deff_bibo_from_config, normalize_oe_sequence
from bibojsa.materials.bibo import BIBO, D_COEFFS_PM_PER_V
from bibojsa.phasematch import ThreeWaveConfig
from bibojsa.sellmeier import optic_axis_angle_deg

D11 = D_COEFFS_PM_PER_V["d11"]
D12 = D_COEFFS_PM_PER_V["d12"]
D13 = D_COEFFS_PM_PER_V["d13"]
D14 = D_COEFFS_PM_PER_V["d14"]


# ----------------------------------------------------------------------
# 1. Analytic spot checks of the Table 1 expressions
# ----------------------------------------------------------------------

def test_xy_ooe_expression():
    phi = np.radians(30.0)
    val = deff_bibo_pm_per_v("xy", "e", "o", "o", phi)  # pump e, daughters o
    assert val == pytest.approx(-D13 * np.sin(phi), rel=1e-12)
    assert abs(val) == pytest.approx(1.76 * 0.5, abs=1e-9)


def test_xy_eoe_expression():
    phi = np.radians(30.0)
    val = deff_bibo_pm_per_v("xy", "e", "o", "e", phi)  # eoe: idler e, signal o, pump e
    assert val == pytest.approx(D14 * np.sin(2 * phi), rel=1e-12)


def test_yz_eeo_expression_upper_sign():
    theta = np.radians(40.0)
    val = deff_bibo_pm_per_v("yz", "o", "e", "e", theta, phi_sign=+1)
    expected = (-D12 * np.cos(theta) ** 2 - D13 * np.sin(theta) ** 2
                + D14 * np.sin(2 * theta))
    assert val == pytest.approx(expected, rel=1e-12)


def test_yz_oeo_is_exactly_zero():
    for theta_deg in (10.0, 45.0, 100.0, 170.0):
        assert deff_bibo_pm_per_v("yz", "o", "e", "o", np.radians(theta_deg)) == 0.0


def test_xz_expressions():
    lam_p = 0.8
    th_in = np.radians(90.0)   # between the optic axes
    th_out = np.radians(11.0)  # outside them (Omega(0.8um) ~ 25.9 deg)
    assert deff_bibo_pm_per_v("xz", "o", "e", "e", th_in, lam_p) == pytest.approx(
        -D14 * np.sin(2 * th_in), abs=1e-12)                       # eeo
    assert deff_bibo_pm_per_v("xz", "o", "e", "o", th_in, lam_p) == pytest.approx(
        D12 * np.cos(th_in), abs=1e-12)                            # oeo
    assert deff_bibo_pm_per_v("xz", "e", "o", "o", th_out, lam_p) == pytest.approx(
        -D12 * np.cos(th_out), rel=1e-12)                          # ooe
    assert deff_bibo_pm_per_v("xz", "e", "o", "e", th_out, lam_p) == pytest.approx(
        -D14 * np.sin(2 * th_out), rel=1e-12)                      # eoe


def test_equivalent_sequences_share_expression():
    """LPR2010 Table 1: eoe = oee and oeo = eoo."""
    phi = np.radians(25.0)
    assert (deff_bibo_pm_per_v("xy", "e", "o", "e", phi)
            == deff_bibo_pm_per_v("xy", "e", "e", "o", phi))  # eoe vs oee
    theta = np.radians(70.0)
    assert (deff_bibo_pm_per_v("xz", "o", "e", "o", theta, 0.8)
            == deff_bibo_pm_per_v("xz", "o", "o", "e", theta, 0.8))  # oeo vs eoo
    assert normalize_oe_sequence("o", "e", "e") == "eoe"  # oee -> eoe
    assert normalize_oe_sequence("e", "o", "e") == "eoe"
    assert normalize_oe_sequence("e", "o", "o") == "oeo"


def test_invalid_combinations_raise():
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("xy", "e", "e", "e", 0.5)   # eee not in Table 1
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("xy", "o", "e", "e", 0.5)   # eeo not in x-y plane
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("yz", "e", "o", "o", 0.5)   # ooe not in y-z plane
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("bad", "e", "o", "o", 0.5)


def test_xz_optic_axis_range_enforcement():
    """Table 1's x-z angle ranges: '+' types (eeo, oeo) live between the
    optic axes (Omega < theta < 180 - Omega), '-' types (ooe, eoe)
    outside. Omega(800 nm) ~ 25.9 deg."""
    lam_p = 0.8
    omega = optic_axis_angle_deg(lam_p, BIBO["x"], BIBO["y"], BIBO["z"])
    assert 25.0 < omega < 27.0  # sanity on the 800 nm optic-axis half-angle

    # '+' type outside its range -> raises; with enforce_range=False -> value
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("xz", "o", "e", "e", np.radians(11.0), lam_p)
    bare = deff_bibo_pm_per_v("xz", "o", "e", "e", np.radians(11.0), lam_p,
                               enforce_range=False)
    assert bare == pytest.approx(-D14 * np.sin(np.radians(22.0)), rel=1e-12)

    # '-' type inside the between-axes region -> raises
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("xz", "e", "o", "o", np.radians(90.0), lam_p)

    # missing pump wavelength with enforcement on -> raises
    with pytest.raises(ValueError):
        deff_bibo_pm_per_v("xz", "e", "o", "o", np.radians(11.0))


# ----------------------------------------------------------------------
# 2. First-principles d-tensor contraction cross-check
# ----------------------------------------------------------------------

# 3x6 d-matrix (Voigt), point group 2 with two-fold axis || x, Kleinman:
#   row x: [d11, d12, d13, d14,   0,   0]
#   row y: [  0,   0,   0,   0, d25, d26]   (d25 = d14, d26 = d12)
#   row z: [  0,   0,   0,   0, d35, d36]   (d35 = d13, d36 = d14)
_D_MATRIX = np.array([
    [D11, D12, D13, D14, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, D14, D12],
    [0.0, 0.0, 0.0, 0.0, D13, D14],
])


def _contract(e_pump, e_sig, e_idl):
    """|deff| = |e_p . d : (e_s e_i)| via the Voigt contracted product."""
    V = np.array([
        e_sig[0] * e_idl[0],
        e_sig[1] * e_idl[1],
        e_sig[2] * e_idl[2],
        e_sig[1] * e_idl[2] + e_sig[2] * e_idl[1],
        e_sig[0] * e_idl[2] + e_sig[2] * e_idl[0],
        e_sig[0] * e_idl[1] + e_sig[1] * e_idl[0],
    ])
    return abs(float(e_pump @ (_D_MATRIX @ V)))


def _unit_vectors(plane, angle):
    """(o_hat, e_hat) for propagation at `angle` in `plane` (phi=0 for xz,
    phi=+90deg for yz, theta=90deg for xy), one consistent sign choice."""
    c, s = np.cos(angle), np.sin(angle)
    if plane == "xz":       # k = (s, 0, c)
        return np.array([0.0, 1.0, 0.0]), np.array([c, 0.0, -s])
    if plane == "yz":       # k = (0, s, c)
        return np.array([1.0, 0.0, 0.0]), np.array([0.0, c, -s])
    # 'xy': k = (cos phi, sin phi, 0), angle = phi
    return np.array([0.0, 0.0, 1.0]), np.array([-s, c, 0.0])


_ALL_ROWS = [
    ("xy", "ooe"), ("xy", "eoe"),
    ("yz", "eeo"), ("yz", "oeo"),
    ("xz", "eeo"), ("xz", "oeo"), ("xz", "ooe"), ("xz", "eoe"),
]


@pytest.mark.parametrize("plane,seq", _ALL_ROWS)
def test_table1_matches_tensor_contraction(plane, seq):
    """|deff| from the Table 1 expressions must equal the direct d-tensor
    contraction for every row, over a dense angle grid spanning both
    octants (0..180 deg for yz/xz thetas, 0..360 deg for xy phi).
    Magnitudes compared (e/o unit-vector sign conventions are only fixed
    up to +-1 per wave)."""
    pol_idl, pol_sig, pol_pump = seq[0], seq[1], seq[2]
    angles = np.radians(np.linspace(1.0, 359.0 if plane == "xy" else 179.0, 73))
    for a in angles:
        o_hat, e_hat = _unit_vectors(plane, a)
        vecs = {"o": o_hat, "e": e_hat}
        expected = _contract(vecs[pol_pump], vecs[pol_sig], vecs[pol_idl])
        got = abs(deff_bibo_pm_per_v(plane, pol_pump, pol_sig, pol_idl, a,
                                      enforce_range=False))
        assert got == pytest.approx(expected, abs=1e-10), (
            f"{plane}/{seq} at {np.degrees(a):.1f} deg: table={got:.6f}, "
            f"contraction={expected:.6f} pm/V"
        )


# ----------------------------------------------------------------------
# 3. Literature benchmarks / monoclinic octant inequivalence
# ----------------------------------------------------------------------

def test_yz_shg_benchmark_and_octant_inequivalence():
    """BIBO's flagship type-I SHG configuration (1064 -> 532 nm) lives in
    the y-z plane at theta ~ 168.9 deg (second octant) with |deff| ~ 3.6-3.7
    pm/V; the mirror angle in the first octant (~11.1 deg) gives only
    ~2.4 pm/V. The ratio ~1.5 is the monoclinic octant inequivalence."""
    th2 = np.radians(168.9)
    th1 = np.radians(180.0 - 168.9)
    d2 = abs(deff_bibo_pm_per_v("yz", "o", "e", "e", th2))
    d1 = abs(deff_bibo_pm_per_v("yz", "o", "e", "e", th1))
    assert d2 == pytest.approx(3.644, abs=0.01)
    assert d1 == pytest.approx(2.389, abs=0.01)
    assert d2 / d1 > 1.4  # clearly inequivalent octants


def test_yz_phi_sign_mirror_symmetry():
    """2/m symmetry: phi = -90 deg at theta is equivalent (in |deff|) to
    phi = +90 deg at 180 deg - theta."""
    for theta_deg in (11.1, 40.0, 75.0, 130.0):
        a = abs(deff_bibo_pm_per_v("yz", "o", "e", "e", np.radians(theta_deg),
                                    phi_sign=-1))
        b = abs(deff_bibo_pm_per_v("yz", "o", "e", "e",
                                    np.radians(180.0 - theta_deg), phi_sign=+1))
        assert a == pytest.approx(b, rel=1e-12)


def test_milestone3_spdc_design_point_deff():
    """The milestone-3 BIBO SPDC configuration (x-z plane, ooe, pump 800 nm,
    theta ~ 10.96 deg at the lambda_s = 1.34 um EPM candidate) should give
    |deff| = d12 cos(theta) ~ 3.14 pm/V -- BIBO's key practical advantage
    over e.g. BBO (~2 pm/V) for bright pair generation."""
    config = ThreeWaveConfig(plane="xz", pol_pump="e", pol_signal="o", pol_idler="o")
    theta = np.radians(10.96)
    val = deff_bibo_from_config(config, theta, lambda_pump_um=0.8)
    assert abs(val) == pytest.approx(D12 * np.cos(theta), rel=1e-12)
    assert abs(val) == pytest.approx(3.14, abs=0.01)
