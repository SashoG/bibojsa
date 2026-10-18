"""The second half of each type-II phase-matching surface (17 Oct 2026).

`solve_theta_grid` orders lambda_i >= lambda_s, so "oeo" (idler o, signal
e) covers only configurations with the e daughter shorter; "eoo" is the
other half.  These tests pin the continuity of the level curves through
the degenerate point and two reference points of the previously untraced
half, found independently in physics_check/ (bibo_ref.py):

  x-z eoo, pump 640 nm: o-wave 1.130 um, e-wave 1.476 um, theta = 50.25 deg,
      k_p' = k_e' (the continuation of point A's k_p' = k_e' curve);
  x-y oee, pump 462 nm: e-wave 0.618 um, o(z)-wave 1.827 um, phi = 46.99 deg,
      k_p' = k_o' (asymmetric GVM with |d_eff| = d14 sin 2phi = 1.66 pm/V).
"""
import numpy as np

from bibojsa import landscape as L
from bibojsa.deff import deff_bibo_pm_per_v
from bibojsa.materials.bibo import BIBO

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]


def test_region_ok_partner_sequences():
    th = np.array([10.0, 46.0, 80.0])
    assert L.region_ok("xz", "eoo", th, 0.64, X, Y, Z).tolist() == [False, True, True]
    assert L.region_ok("xz", "oee", th, 0.64, X, Y, Z).tolist() == [True, False, False]


def test_xz_eoo_continues_the_kp_eq_ke_curve():
    # signal = o (shorter), idler = e (longer)
    th, li = L.solve_theta_grid(0.640, 1.130, "xz", "eoo", X, Y, Z)
    assert abs(np.degrees(th) - 50.25) < 0.05
    assert abs(li - 1.476) < 0.002
    g = L.gvm_functions(0.640, 1.130, li, th, "xz", "eoo", X, Y, Z)
    assert abs(g["g90"]) * 1e12 < 1.0            # k_p' = k_i' = k_e', fs/mm
    assert abs(g["g0"]) * 1e12 > 150.0
    assert abs(abs(g["theta_pm"]) - 90.0) < 0.5      # +-90 deg are the same orientation
    # same |d_eff| expression as oeo
    d = deff_bibo_pm_per_v("xz", "o", "o", "e", th, lambda_pump_um=0.640)
    assert abs(abs(d) - 3.2 * np.cos(th)) < 1e-12


def test_xz_eoo_traced_curve_passes_through_point_A():
    # the eoo 90-deg curve and the oeo 0-deg curve share the degenerate member
    p_oeo = L.degenerate_gvm_point("xz", "oeo", 0, X, Y, Z, (0.60, 0.80))
    p_eoo = L.degenerate_gvm_point("xz", "eoo", 90, X, Y, Z, (0.60, 0.80))
    assert p_oeo is not None and p_eoo is not None
    assert abs(p_oeo.lam_p_um - p_eoo.lam_p_um) < 1e-6
    assert abs(p_oeo.theta_deg - p_eoo.theta_deg) < 1e-4   # root-finder tolerance
    # and the eoo curve exists below 682 nm, where the oeo curve does not
    pts = L.trace_gvm_curve("xz", "eoo", 90, X, Y, Z, [0.640], n_ls=200)
    assert pts, "no eoo 90-deg crossing at 640 nm"
    assert abs(pts[0].lam_s_um - 1.130) < 0.01 and abs(pts[0].lam_i_um - 1.476) < 0.01
    assert not L.trace_gvm_curve("xz", "oeo", 0, X, Y, Z, [0.640], n_ls=200)


def test_xy_oee_asymmetric_point():
    # idler = o (z-polarized, longer), signal = e (in-plane, shorter)
    pts = L.trace_gvm_curve("xy", "oee", 90, X, Y, Z, [0.462], n_ls=200)
    assert pts, "no x-y oee 90-deg crossing at 462 nm"
    p = pts[0]
    assert abs(p.lam_s_um - 0.618) < 0.005 and abs(p.lam_i_um - 1.827) < 0.01
    assert abs(p.theta_deg - 46.99) < 0.1
    d = deff_bibo_pm_per_v("xy", "e", "e", "o", np.radians(p.theta_deg))
    assert abs(abs(d) - 1.66) < 0.01
    # the author's half of x-y eoe has no GVM crossing here
    assert not L.trace_gvm_curve("xy", "eoe", 90, X, Y, Z, [0.462], n_ls=200)
