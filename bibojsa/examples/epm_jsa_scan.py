"""JSA / Schmidt / heralded-purity study at the BIBO x-z ooe symmetric-GVM
(theta_PM = 45 deg) operating point: pump 800 nm (e), signal 1.3409 um (o),
idler 1.9831 um (o), theta_cut = 10.9577 deg.

Because this is a *triple* GVM point (all three pairwise group-delay
differences are only a few fs/mm), the first-order phase-matching
bandwidth is enormous and the sinc's true extent is governed by
second-order dispersion. Delta k is therefore evaluated exactly (full
Sellmeier, no Taylor expansion) and grid windows are chosen adaptively
from the actual decay of |f|, clipped to the Sellmeier validity range
(0.32-3.09 um) -- any support truncated by the validity limit is
reported per point (`edge` in the JSON output = max |f| on the grid
boundary relative to max |f|).

Outputs:
  examples/epm_records.jsonl         per-point cache (resumable runs)
  examples/epm_scan_results.json     purity vs (L, pump bandwidth) records
  figures/epm_jsa_maps.png           |JSA| at representative (L, bw)
  figures/epm_purity_scans.png       purity vs bandwidth and vs L
  figures/epm_schmidt.png            Schmidt spectrum at the best point
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig
from bibojsa.jsa import build_jsa, fwhm_omega_from_fwhm_wavelength
from bibojsa.sellmeier import lambda_um_to_omega, omega_to_lambda_um

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
CONFIG = ThreeWaveConfig(plane="xz", pol_pump="e", pol_signal="o", pol_idler="o")

LAMBDA_P = 0.800
LAMBDA_S = 1.3409
LAMBDA_I = 1.9831
THETA = np.radians(10.9577)

WP0 = lambda_um_to_omega(LAMBDA_P)
WS0 = lambda_um_to_omega(LAMBDA_S)
WI0 = lambda_um_to_omega(LAMBDA_I)

# Sellmeier validity clip (small safety margin inside 0.32-3.09 um)
W_MAX = lambda_um_to_omega(0.33)   # highest frequency allowed
W_MIN = lambda_um_to_omega(3.05)   # lowest frequency allowed


def adaptive_axes(L_m, fwhm_nm, n):
    """Choose (omega_s, omega_i) axes wide enough to contain the JSA.

    The JSA support is the pump stripe (width ~ pump sigma along
    omega_s + omega_i) extended along the antidiagonal until the exact
    sinc kills it. Probe |sinc(Delta k L/2)| along the antidiagonal
    through the design point to find that extent, then clip to Sellmeier
    validity.
    """
    from bibojsa.phasematch import delta_k

    def first_ext(deltas):
        """Extent along the antidiagonal: first delta (from 0 outward)
        where |sinc| < 0.05 or the Sellmeier validity limit is hit;
        vectorized over the whole probe array."""
        ws, wi = WS0 + deltas, WI0 - deltas
        valid = (ws <= W_MAX) & (ws >= W_MIN) & (wi <= W_MAX) & (wi >= W_MIN)
        x = np.where(valid,
                     delta_k(THETA, np.where(valid, ws + wi, WP0),
                             np.where(valid, ws, WS0), np.where(valid, wi, WI0),
                             CONFIG, X, Y, Z) * L_m / 2.0,
                     0.0)
        dead = (~valid) | (np.abs(np.sinc(x / np.pi)) < 0.05)
        idx = np.argmax(dead)
        return abs(deltas[idx]) if dead.any() else abs(deltas[-1])

    probe = np.linspace(0, 6e14, 400)
    ext_plus = first_ext(probe)
    ext_minus = first_ext(-probe)

    sigma = fwhm_omega_from_fwhm_wavelength(fwhm_nm * 1e-3, LAMBDA_P)
    pad = 3.0 * sigma  # pump-stripe thickness margin
    ws_lo = max(WS0 - 1.15 * ext_minus - pad, W_MIN)
    ws_hi = min(WS0 + 1.15 * ext_plus + pad, W_MAX)
    wi_lo = max(WI0 - 1.15 * ext_plus - pad, W_MIN)
    wi_hi = min(WI0 + 1.15 * ext_minus + pad, W_MAX)
    return np.linspace(ws_lo, ws_hi, n), np.linspace(wi_lo, wi_hi, n)


def run_point(L_mm, fwhm_nm, n=140):
    L_m = L_mm * 1e-3
    ws_axis, wi_axis = adaptive_axes(L_m, fwhm_nm, n)
    fw = fwhm_omega_from_fwhm_wavelength(fwhm_nm * 1e-3, LAMBDA_P)
    res = build_jsa(ws_axis, wi_axis, WP0, THETA, CONFIG, X, Y, Z,
                    L_m, fw, pump_shape="gaussian", include_phase=True)
    A = np.abs(res.f)
    m = A.max()
    edge = max(A[0, :].max(), A[-1, :].max(), A[:, 0].max(), A[:, -1].max()) / m if m > 0 else 0.0
    return res, edge


_CACHE_PATH = "examples/epm_records.jsonl"


def cached_point(L_mm, fwhm_nm, n=140):
    """run_point with a JSONL cache so interrupted runs resume cheaply."""
    key = {"L_mm": float(L_mm), "fwhm_nm": float(fwhm_nm), "n": n}
    if os.path.exists(_CACHE_PATH):
        with open(_CACHE_PATH) as fh:
            for line in fh:
                rec = json.loads(line)
                if all(rec.get(k) == v for k, v in key.items()):
                    return rec
    res, edge = run_point(L_mm, fwhm_nm, n)
    rec = dict(key, purity=res.purity, K=res.schmidt_number, edge=edge)
    with open(_CACHE_PATH, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return rec


def main():
    os.makedirs("figures", exist_ok=True)
    records = []

    # --- purity vs pump bandwidth at three lengths -----------------------
    bw_grid = [1, 2, 4, 8, 16, 24, 32, 40]        # nm FWHM at 800 nm
    L_show = [2.0, 5.0, 10.0]                      # mm
    for L in L_show:
        for bw in bw_grid:
            rec = cached_point(L, bw)
            records.append(rec)
            print(f"L={L:5.1f} mm  bw={bw:5.1f} nm  P={rec['purity']:.4f}  "
                  f"K={rec['K']:6.2f}  edge={rec['edge']:.3f}", flush=True)

    # --- purity vs L at the per-L best bandwidth from the grid ----------
    L_grid = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    best_bw = {}
    purity_vs_L = []
    for L in L_grid:
        best = None
        for bw in bw_grid:
            rec = cached_point(L, bw)
            records.append(rec)
            if best is None or rec["purity"] > best[1]:
                best = (bw, rec["purity"], rec["edge"])
        best_bw[L] = best
        purity_vs_L.append((L, best[1]))
        print(f"L={L:5.1f} mm  best bw={best[0]:5.1f} nm  P={best[1]:.4f}  "
              f"edge={best[2]:.3f}", flush=True)

    # --- convergence check at the overall best point ---------------------
    L_best = max(best_bw, key=lambda L: best_bw[L][1])
    bw_best = best_bw[L_best][0]
    conv = []
    for n in (100, 140, 200, 260):
        rec = cached_point(L_best, bw_best, n=n)
        conv.append({"n": n, "purity": rec["purity"], "edge": rec["edge"]})
        print(f"conv n={n}: P={rec['purity']:.5f} edge={rec['edge']:.3f}", flush=True)

    with open("examples/epm_scan_results.json", "w") as fh:
        json.dump({"design_point": {"lambda_p_um": LAMBDA_P, "lambda_s_um": LAMBDA_S,
                                     "lambda_i_um": LAMBDA_I,
                                     "theta_deg": float(np.degrees(THETA))},
                   "records": records, "purity_vs_L_best_bw": purity_vs_L,
                   "best": {"L_mm": L_best, "fwhm_nm": bw_best},
                   "convergence": conv}, fh, indent=1)

    # --- figures ----------------------------------------------------------
    # JSA maps at three (L, best bw) combos
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, L in zip(axes, L_show):
        bw = max((r for r in records if r["L_mm"] == L),
                 key=lambda r: r["purity"])["fwhm_nm"]
        res, _ = run_point(L, bw, n=200)
        ls = omega_to_lambda_um(res.omega_s)
        li = omega_to_lambda_um(res.omega_i)
        ax.pcolormesh(li, ls, np.abs(res.f), shading="auto", rasterized=True)
        ax.set_xlabel("idler wavelength (um)")
        ax.set_ylabel("signal wavelength (um)")
        ax.set_title(f"L = {L:g} mm, pump FWHM = {bw:g} nm\nP = {res.purity:.3f}")
    fig.suptitle("|JSA|, BIBO x-z ooe, 800 nm pump, theta = 10.96 deg (theta_PM = 45 deg)")
    fig.savefig("figures/epm_jsa_maps.png", dpi=160)
    plt.close(fig)

    # purity scans
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for L in L_show:
        rs = [r for r in records if r["L_mm"] == L][:len(bw_grid)]
        ax1.plot([r["fwhm_nm"] for r in rs], [r["purity"] for r in rs],
                 marker="o", label=f"L = {L:g} mm")
    ax1.set_xlabel("pump intensity FWHM (nm)")
    ax1.set_ylabel("unfiltered heralded purity")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    Ls, Ps = zip(*purity_vs_L)
    ax2.semilogx(Ls, Ps, marker="s")
    ax2.set_xlabel("crystal length (mm)")
    ax2.set_ylabel("purity at best scanned bandwidth")
    ax2.grid(True, alpha=0.3, which="both")
    fig.suptitle("Heralded purity at the 1.34 um triple-GVM point (no heralding filter)")
    fig.savefig("figures/epm_purity_scans.png", dpi=160)
    plt.close(fig)

    # Schmidt spectrum at the best point
    res, _ = run_point(L_best, bw_best, n=200)
    fig, ax = plt.subplots(figsize=(5, 3.5), constrained_layout=True)
    ax.bar(range(1, 21), res.schmidt_lambdas[:20])
    ax.set_xlabel("Schmidt mode index n")
    ax.set_ylabel("lambda_n")
    ax.set_title(f"Schmidt spectrum, L = {L_best:g} mm, FWHM = {bw_best:g} nm\n"
                 f"P = {res.purity:.3f}, K = {res.schmidt_number:.2f}")
    fig.savefig("figures/epm_schmidt.png", dpi=160)
    plt.close(fig)

    print("done", flush=True)


if __name__ == "__main__":
    main()
