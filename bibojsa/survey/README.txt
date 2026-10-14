BIBO GVM/JSA survey data
========================

Data underlying "Mapping group-velocity matching in monoclinic BiB3O6
for high-purity photon pairs" (Data File 1).  All files are plain CSV
with a header row.  Regenerated 18 September 2026 (v0.6 of the
manuscript).

Conventions
-----------
* Interaction labels are the LPR2010 o/e sequence of
  (lambda_1, lambda_2, lambda_3) with lambda_1 >= lambda_2 > lambda_3,
  i.e. idler-signal-pump.  An interaction whose two daughters carry the
  same label (ooe, eeo) is type I; one whose daughters differ (oeo,
  eoe) is type II.
* theta is the crystal-cut angle (deg): from z in the y-z and x-z
  planes, from x in the x-y plane.  In the x-z plane the mirror
  partner at 180-theta is equivalent and only theta <= 90 is listed.
* theta_PM is the JSA ridge orientation in the (omega_s, omega_i)
  plane measured from the SIGNAL axis,
      tan(theta_PM) = -(k_p' - k_s') / (k_p' - k_i'),
  so theta_PM = 0 means k_p' = k_s', 90 means k_p' = k_i', +45 is
  symmetric GVM / extended phase matching, and -45 (= 135) is the
  ordinary anticorrelated orientation.
* d31 = k_p' - k_i' and d32 = k_p' - k_s', in fs/mm.
* lambda_i_in_transparency flags the idler against BIBO's 1 cm^-1
  infrared absorption onset at 2.60 um (Ghotbi et al. 2006; LPR2010
  Table 4), which is not the same as the 3.09 um validity range of
  the Sellmeier fit.  Rows with the flag 0 are inside the dispersion
  model but their idler is absorbed in centimetre-long crystals
  (manuscript Sec. 5.1); they are reported for completeness.  (Data
  File 1 of earlier drafts used 2.50 um for this flag.)
* Purities are Schmidt purities P = sum(lambda_n^2) of the complex JSA
  including the exp(i dk L/2) phase; "edge" is the boundary amplitude
  max_boundary|f| / max|f|, and values with edge >= 0.15 are
  truncation-limited rather than converged.
* Every filtered purity is accompanied by the fraction of pairs the
  filter transmits.  A purity quoted without its transmission is not
  meaningful.

Landscape (continuous-pump survey, examples/landscape_survey.py)
---------------------------------------------------------------
landscape_curves.csv      Every point on every traced theta_PM =
                          0/45/90 deg level curve, for all eight
                          principal-plane interactions, pump
                          0.34-1.52 um in 2 nm steps.  1160 rows.
landscape_degenerate.csv  The frequency-degenerate member of each
                          family (lambda_s = lambda_i = 2 lambda_p).
landscape_field.csv       theta_PM sampled on a coarse
                          (lambda_p, lambda_s) grid, for contouring.
landscape_coverage.csv    Per-interaction, per-target record of what
                          was found, including the negative results.

Operating points (examples/landscape_points.py)
-----------------------------------------------
points2_summary.csv       One row per selected operating point.
points2/point<id>_purity.csv       purity and filter transmission over
                                   the (L, pump FWHM) grid
points2/point<id>_convergence.csv  purity vs window scale and grid size
points2/point<id>_jsa.csv          |f| in angular-frequency detuning
points2/point<id>_schmidt.csv      Schmidt spectrum
points2/point<id>_brightness.csv   aperture length and brightness proxy
                                   vs waist

Featured points (examples/landscape_headline.py)
------------------------------------------------
headline_points.csv       Table 2 of the manuscript (points A-F, R).
                          purity_T90/T70: heralding filter on the arm
                          named in herald_arm (the narrower marginal:
                          idler at theta_PM = 0, signal at 90; idler by
                          convention at symmetric points) -- since
                          14 Oct 2026; before, always the idler.
headline_filter.csv       Purity vs transmitted heralding flux (same
                          arm convention, column herald_arm).
headline_convergence.csv  Window/grid convergence at the quoted
                          configurations.
headline_jsa_<tag>.csv    |f| in angular-frequency detuning.
family_purity.csv         Best window-converged unfiltered purity along
                          each GVM family, vs pump wavelength.
sellmeier_sensitivity.csv The three degenerate GVM points and the
                          purity at the 0 deg point recomputed under
                          each of the three published BIBO Sellmeier
                          fits (Miyata-Umemura-Kato, Hellwig 2000,
                          Umemura 2007).  The crossing wavelengths move
                          by 1.0-1.3% between fits; the purity moves by
                          5e-4.  The Umemura row reproduces the
                          published values of Jin et al. (2020).

Collection-window and length dependence (examples/window_tail.py)
-----------------------------------------------------------------
window_tail.csv           Unfiltered purity of each featured point on
                          square collection windows of half-width
                          6/8/12/16/24 pump sigma (sigma = 1/e
                          half-width of the pump INTENSITY in angular
                          frequency) at fixed 1.0 rad/ps pixel;
                          least-squares fit P(W) = P_inf + C/W (P_inf,
                          C_per_sigma, residual, fit_points); the
                          manuscript's adaptive window in sigma and nm;
                          the T = 90 % / 70 % filtered purities at 6 and
                          16 sigma with the filter FWHM in nm.  Windows
                          that hit the Sellmeier limit are flagged
                          clipped_* = 1 and excluded from the fit.
length_extension.csv      Points A, B, C at L = 10/15/20/30 mm over pump
                          bandwidths 2-16 nm: P on the 6-sigma window
                          and P_inf.

Superseded
----------
survey_loci.csv, survey_crossings.csv, survey_summary.csv,
convergence_all.csv and points/ are the earlier fixed-pump-list survey
(405/515/775/780/800/810 nm).  They are retained because the
manuscript's point R is drawn from them, but the landscape files above
supersede them: the fixed pump list samples the GVM level curves at
arbitrary places and misses all three degenerate members.

Reproducing
-----------
    cd bibojsa
    python examples/landscape_survey.py all       # curves + field
    python examples/landscape_points.py           # per-point JSA/purity
    python examples/landscape_headline.py all     # tables and figures
    python examples/window_tail.py                # window/length tables
    python -m pytest tests                        # 104 tests
