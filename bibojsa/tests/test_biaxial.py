"""Validate the general Fresnel-equation eigenindex solver against:
(1) the uniaxial closed form (KDP: one eigenmode = n_o exactly, for any
    (theta, phi)), and
(2) the principal-plane closed form (BIBO, phi=0 -> xz plane, phi=pi/2
    -> yz-like cut within phi sweep, etc.), which is what phasematch.py
    and gvm.py actually use.
Also sanity-checks walk-off: zero along principal axes, nonzero between.
"""
import numpy as np
import pytest

from bibojsa.materials.kdp import KDP
from bibojsa.materials.bibo import BIBO
from bibojsa.biaxial import eigenindices_general, principal_plane_index, walkoff_angle

LAM_KDP = 0.6328
LAM_BIBO = 0.8


@pytest.mark.parametrize("theta,phi", [(0.3, 0.7), (1.1, 0.2), (0.5, 1.4), (0.9, 0.9)])
def test_uniaxial_limit_gives_ordinary_index_exactly(theta, phi):
    """For KDP (nx=ny=no, nz=ne), one eigenindex must equal no exactly,
    independent of (theta, phi) -- the defining property of the ordinary
    ray in a uniaxial crystal."""
    n_slow, n_fast = eigenindices_general(LAM_KDP, theta, phi, KDP["x"], KDP["y"], KDP["z"])
    no = KDP["x"].n(LAM_KDP)
    ne = KDP["z"].n(LAM_KDP)
    assert no > ne  # KDP is negative uniaxial
    # no should show up as whichever of (slow, fast) matches, since no > ne
    assert min(abs(n_slow - no), abs(n_fast - no)) < 1e-10


@pytest.mark.parametrize("plane,theta", [("xz", 0.2), ("xz", 0.9), ("xz", 1.3)])
def test_general_solver_matches_principal_plane_xz(plane, theta):
    n_slow, n_fast = eigenindices_general(LAM_BIBO, theta, 0.0, BIBO["x"], BIBO["y"], BIBO["z"])
    n_fixed, n_var = principal_plane_index(LAM_BIBO, theta, plane, BIBO["x"], BIBO["y"], BIBO["z"])
    got = sorted([n_slow, n_fast])
    expected = sorted([n_fixed, n_var])
    assert got == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("theta", [0.2, 0.8, 1.2])
def test_general_solver_matches_principal_plane_yz(theta):
    # yz-plane: direction (0, sin theta, cos theta) -> phi = pi/2
    n_slow, n_fast = eigenindices_general(LAM_BIBO, theta, np.pi / 2, BIBO["x"], BIBO["y"], BIBO["z"])
    n_fixed, n_var = principal_plane_index(LAM_BIBO, theta, "yz", BIBO["x"], BIBO["y"], BIBO["z"])
    got = sorted([n_slow, n_fast])
    expected = sorted([n_fixed, n_var])
    assert got == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("phi", [0.2, 0.8, 1.2])
def test_general_solver_matches_principal_plane_xy(phi):
    # xy-plane: direction (sin(theta=pi/2)cos(phi), sin(theta=pi/2)sin(phi), 0)
    n_slow, n_fast = eigenindices_general(LAM_BIBO, np.pi / 2, phi, BIBO["x"], BIBO["y"], BIBO["z"])
    n_fixed, n_var = principal_plane_index(LAM_BIBO, phi, "xy", BIBO["x"], BIBO["y"], BIBO["z"])
    got = sorted([n_slow, n_fast])
    expected = sorted([n_fixed, n_var])
    assert got == pytest.approx(expected, abs=1e-9)


def test_walkoff_vanishes_along_principal_axes():
    rho0 = walkoff_angle(LAM_BIBO, 1e-4, "xz", BIBO["x"], BIBO["y"], BIBO["z"])
    rho90 = walkoff_angle(LAM_BIBO, np.pi / 2 - 1e-4, "xz", BIBO["x"], BIBO["y"], BIBO["z"])
    rho_mid = walkoff_angle(LAM_BIBO, np.pi / 4, "xz", BIBO["x"], BIBO["y"], BIBO["z"])
    assert abs(np.degrees(rho0)) < 0.1
    assert abs(np.degrees(rho90)) < 0.1
    assert abs(np.degrees(rho_mid)) > 1.0  # clearly nonzero in between


def test_principal_plane_rejects_bad_plane_name():
    with pytest.raises(ValueError):
        principal_plane_index(LAM_BIBO, 0.5, "bad", BIBO["x"], BIBO["y"], BIBO["z"])
