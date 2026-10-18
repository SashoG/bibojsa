# Prompt for Claude Code — BIBO SPDC / JSA computational skeleton

Copy everything below the line into Claude Code (or attach this file and say "execute this brief").

---

## Project brief

You are building the computational core for a theory paper with the working title *"Group-velocity-matched photon-pair generation in the monoclinic crystal BiB₃O₆: joint spectral engineering and heralded-photon purity."* I am a physicist (doctoral background in quantum electronics / nonlinear optics) returning to research; the physics conventions below are authoritative — follow them exactly and cite them in docstrings.

Create a Python package named `bibojsa` with tests, example notebooks, and plots. Prioritize correctness and validation against known literature values over features.

## Physics scope

We compute, for spontaneous parametric down-conversion (SPDC) in bulk birefringent crystals:

1. **Dispersion**: refractive indices from Sellmeier equations; group indices k′(ω) = ∂k/∂ω; GVD k″(ω).
2. **Biaxial optics**: for propagation direction (θ, φ) in the dielectric frame xyz with nx < ny < nz, solve the Fresnel equation of wave normals for the two eigenindices (slow/fast); compute spatial walk-off angles for each eigenwave.
3. **Phase matching**: collinear SPDC ωp = ωs + ωi with Δk = kp − ks − ki = 0. In the principal planes, classify interactions by the o/e convention of Petrov et al., Laser Photon. Rev. 4, 53 (2010) [henceforth LPR2010], Table 1, ordering λ1 ≥ λ2 > λ3 (idler, signal, pump). Implement deff per that table (Kleinman symmetry) with coefficients d11 = 2.53, d12 = d26 = 3.2, d13 = d35 = −1.76, d14 = d25 = d36 = 1.66 pm/V (xyz dielectric frame).
4. **GVM classification**: phase-matching-function orientation tan θPM = −(k′p − k′s)/(k′p − k′i). Find loci of θPM = 0° and 90° (asymmetric GVM, k′p = k′s or k′p = k′i) and θPM = 45° (symmetric GVM / extended phase matching, k′p = (k′s + k′i)/2).
5. **Joint spectral amplitude**: f(ωs, ωi) = α(ωs + ωi) · φ(ωs, ωi), with pump envelope α Gaussian or sech² (parametrized by FWHM bandwidth or pulse duration) and φ = sinc(ΔkL/2)·exp(iΔkL/2). Grid f on a configurable (Ns × Ni) window, SVD → Schmidt coefficients λn, heralded purity P = Σλn², Schmidt number K = 1/P. Optional Gaussian spectral filter on the heralding arm before decomposition. Support both amplitude conventions (with/without the exp(iΔkL/2) phase) via a flag, defaulting to including it.

## Crystal data (hard-code, with unit tests)

**BIBO** (Miyata–Umemura–Kato fit, reproduced as Table 2 of LPR2010; λ in µm, valid 0.32–3.09 µm), n² = A + B/(λ² − C) − Dλ²:

| axis | A | B [µm²] | C [µm²] | D [µm⁻²] |
|---|---|---|---|---|
| nx | 3.0759 | 0.03169 | 0.03323 | 0.01402 |
| ny | 3.1698 | 0.03666 | 0.03599 | 0.01819 |
| nz | 3.6546 | 0.05116 | 0.03713 | 0.02299 |

Unit tests against LPR2010 Table 3: at 632.8 nm nx = 1.77668, ny = 1.80641, nz = 1.94582; at 1064 nm nx = 1.75752, ny = 1.78400, nz = 1.91711 (tolerance 2×10⁻⁵). Also test the optic-axis angle: 2Ω = 52.6° at 1064 nm, with Ω = arcsin[(nz/ny)·√((ny²−nx²)/(nz²−nx²))].

**Validation crystals** (add with literature Sellmeier fits, cite source in docstring): BBO, KDP, KTP (bulk, Fradkin/Kato-type fits including ny and nz), and optionally a QPM-KTP mode where Δk includes 2π/Λ.

## Validation milestones (do these before touching BIBO results)

1. **KDP asymmetric GVM**: reproduce that type-II KDP pumped near 415 nm produces degenerate photons at 830 nm with θPM ≈ 0°/90°, i.e. the Mosley et al. PRL 100, 133601 (2008) configuration; heralded purity of the broadband photon high without filtering.
2. **KTP symmetric GVM**: reproduce extended phase matching for type-II KTP with degenerate photons near 1582 nm (θPM ≈ 45°), unfiltered purity ≈ 0.80–0.85 for a bandwidth-matched Gaussian pump (sinc-sidelobe-limited) — König & Wong APL 84, 1644 (2004); Evans et al. PRL 105, 253601 (2010); Graffitti et al. PRA 98, 053811 (2018).
3. **Sanity vs LPR2010**: for BIBO type-I ooe in the x–z plane pumped at 800 nm, reproduce Fig. 10b qualitatively: GVM parameter Δ21 → small near λsignal ≈ 1340 nm (signal and idler co-propagating with the pump). This point is our first EPM candidate.

## Exploration targets (after validation)

For pump wavelengths 405, 515, 775, 780, 800, 810 nm: scan all principal-plane interaction types of BIBO (Table 1 of LPR2010, both relevant octants: 0 ≤ θ ≤ 180° where applicable), find phase-matched (θ, λs) solutions, compute θPM, |deff|, walk-off angles for all three waves, and produce:

- Maps of θPM over the phase-matching loci, highlighting crossings of 0°, 45°, 90°.
- For each highlighted operating point: JSA plots, Schmidt spectra, purity vs crystal length (0.1–10 mm log grid) and vs pump bandwidth, unfiltered and with heralding-arm filtering.
- Walk-off-limited effective length estimate: L_eff at which transverse walk-off equals the beam waist, for waists 20–200 µm; brightness proxy deff²·min(L, L_eff)².

## Engineering requirements

- Python ≥3.11, numpy/scipy/matplotlib; no exotic dependencies. Type hints, docstrings with the governing equations and literature references.
- Frequencies internally in rad/s; user-facing API in vacuum wavelength (µm). Derivatives k′, k″ by ≥4th-order central finite differences on k(ω) with step-size convergence check (or use automatic differentiation if you prefer — but then still test against finite differences).
- Package layout: `bibojsa/{sellmeier.py, biaxial.py, phasematch.py, gvm.py, jsa.py, materials/}`, `tests/` (pytest), `notebooks/01_validation_kdp_ktp.ipynb`, `notebooks/02_bibo_survey.ipynb`.
- Every physics function unit-tested against the numbers above. CI-style `pytest` must pass before you present results.
- Plots publication-lean: label axes with units, no styling gimmicks.

## Cautions

- BIBO is monoclinic point group 2: deff in the plane normal to the two-fold axis is neither even nor odd in θ; two octants are inequivalent (LPR2010 Sec. 2–3). Do not assume orthorhombic symmetry shortcuts.
- The xyz dielectric frame rotates slightly with wavelength about the b-axis (≤ ±1°); we neglect this (as LPR2010 does for deff) — note it in the docs as an approximation.
- Neglect walk-off corrections to deff (note as approximation); walk-off enters only via L_eff.
- Near degeneracy for type-I, first- and second-order GVM terms can vanish simultaneously — JSA grids must be wide enough to capture higher-order curvature of the phase-matching band; add a convergence check on purity vs grid size and window.

## Working style

Proceed milestone by milestone; after each validation milestone, stop and show me the comparison table (computed vs literature) before continuing. Ask before adding any physics not specified here.
