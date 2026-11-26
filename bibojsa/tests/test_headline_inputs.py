"""Consistency of the featured-point inputs with the released tables (v1.2.0).

The rows of survey/headline_points.csv were computed from the values that
`landscape.trace_gvm_curve` returns at full precision.  bibojsa v1.1.0 listed
the nondegenerate members G, H, J, K, M in examples/landscape_headline.py
rounded to six digits; at H that was enough to move the adaptive collection
window from 6.15 to 5.90 pump sigma and the unfiltered purity from 0.95124
to 0.95152.  These tests pin the inputs to the traced values, the released
purity of H to its window, and two statements of the manuscript's first
revision: the fold of the x-z k_p' = k_e' curve at 722.6 nm and the absence
of any phase-matched x-y ooe configuration.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "examples"))

from landscape_headline import FEATURED, as_pt  # noqa: E402
from landscape_points import cfg_for, run_case  # noqa: E402

from bibojsa import landscape as L  # noqa: E402
from bibojsa.jsa import (fwhm_omega_from_fwhm_wavelength,  # noqa: E402
                         gaussian_pump_sigma)
from bibojsa.materials.bibo import BIBO  # noqa: E402
from bibojsa.phasematch import delta_k  # noqa: E402
from bibojsa.sellmeier import lambda_um_to_omega  # noqa: E402

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
BY_TAG = {f[0]: f for f in FEATURED}
TRACED = ("G", "H", "J", "K", "M")      # full-precision members


def test_featured_list_is_in_tag_order():
    tags = [f[0] for f in FEATURED]
    assert tags == sorted(tags)
    assert tags == list("ABCDEFGHJKMR")


def test_traced_members_conserve_energy_and_are_phase_matched():
    for tag in TRACED:
        _, plane, seq, _, pnm, ls, li, th, _, L_mm, _, _ = BY_TAG[tag]
        lam_p = pnm / 1e3
        assert abs(1.0 / lam_p - 1.0 / ls - 1.0 / li) < 1e-9, tag
        wp, ws = lambda_um_to_omega(lam_p), lambda_um_to_omega(ls)
        half_phase = delta_k(np.radians(th), wp, ws, wp - ws,
                             cfg_for(plane, seq), X, Y, Z) * L_mm * 1e-3 / 2
        assert abs(float(half_phase)) < 1e-3, (tag, float(half_phase))


def test_traced_members_are_the_traced_roots():
    for tag in TRACED:
        _, plane, seq, target, pnm, ls, li, th, *_ = BY_TAG[tag]
        pts = L.trace_gvm_curve(plane, seq, target, X, Y, Z, [pnm / 1e3],
                                n_ls=320, n_theta=1201)
        best = min(pts, key=lambda p: abs(p.lam_s_um - ls))
        assert abs(best.lam_s_um - ls) < 1e-9, tag
        assert abs(best.lam_i_um - li) < 1e-9, tag
        assert abs(best.theta_deg - th) < 1e-7, tag


def test_point_H_released_purity_and_window():
    f = BY_TAG["H"]
    res, edge, _, wsa, wia, *_ = run_case(as_pt(f), f[9], f[10], n=200)
    sigma = gaussian_pump_sigma(
        fwhm_omega_from_fwhm_wavelength(f[10] * 1e-3, f[4] / 1e3))
    assert abs(res.purity - 0.95124) < 2e-5
    assert abs(0.5 * (wsa[-1] - wsa[0]) / sigma - 6.15) < 0.01
    assert abs(0.5 * (wia[-1] - wia[0]) / sigma - 6.15) < 0.01
    assert abs(edge - 0.01112) < 1e-4


def test_kp_eq_ke_fold_is_at_722p6_nm():
    def n_points(lam_p):
        return sum(len(L.trace_gvm_curve("xz", "oeo", 0, X, Y, Z, [lam_p],
                                         n_ls=320, branch=b, n_theta=1201))
                   for b in ("first", "second"))
    assert n_points(0.7220) == 2      # two arms below the fold
    assert n_points(0.7226) == 2
    assert n_points(0.7227) == 0      # none above it


def test_xy_ooe_is_never_phase_matched():
    for lam_p in np.arange(0.34, 1.52, 0.02):
        ls = np.linspace(max(lam_p * 1.0005, L.LAM_MIN),
                         min(2 * lam_p * 0.999999, L.LAM_MAX), 80)
        for branch in ("first", "second"):
            th, _ = L.solve_theta_grid(np.full_like(ls, lam_p), ls, "xy",
                                       "ooe", X, Y, Z, n_theta=601,
                                       branch=branch)
            assert not np.isfinite(th).any(), lam_p
