BIBO GVM/JSA survey data
========================

Data underlying "Mapping group-velocity matching in monoclinic BiB3O6
for high-purity photon pairs" (Data File 1).  All files are plain CSV
with a header row.  Regenerated 17 October 2026 (v0.11 of the
manuscript): the survey now covers BOTH halves of every type-II
phase-matching surface (see "Interaction labels" below); the
eight-interaction run of 18 September 2026 is not part of the public
repository (its rows are a strict subset of these; available from the
author).  landscape_field.csv, points2/ and the
superseded fixed-pump files were not regenerated and are unchanged.

Conventions
-----------
* Interaction labels are the LPR2010 o/e sequence of
  (lambda_1, lambda_2, lambda_3) with lambda_1 >= lambda_2 > lambda_3,
  i.e. idler-signal-pump.  An interaction whose two daughters carry the
  same label (ooe, eeo) is type I; one whose daughters differ (oeo,
  eoe) is type II.
* Because lambda_i >= lambda_s is enforced, a type-II sequence fixes
  which daughter is the e wave: "oeo" is the half of the surface on
  which the e daughter is the SHORTER one, "eoo" the half on which it
  is the LONGER one (LPR2010 writes "oeo = eoo" because the d_eff
  expression is shared; the phase matching is not).  Likewise "eoe"
  (o daughter shorter) and "oee" (o daughter longer); in the x-y plane
  the o wave is the z-polarized one.  The GVM level curves are
  continuous across the two halves, with the target label flipping:
  in x-z the physical condition k_p' = k_e' is oeo/target 0 together
  with eoo/target 90, EPM is oeo/45 + eoo/45, and k_p' = k_o' is
  oeo/90 + eoo/0.  Before 17 October 2026 only the first eight
  sequences of LPR2010 Table 1 were traced (oeo, eoe, eeo, ooe per
  plane); the partner sequences eoo and oee were added then.
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
                          0/45/90 deg level curve, for all twelve
                          principal-plane interactions (the eight of
                          LPR2010 Table 1 plus their type-II partners
                          x-z eoo, x-z oee, y-z eoo, x-y oee), pump
                          0.34-1.52 um in 2 nm steps.  2193 rows
                          (1160 before 17 Oct 2026; those rows are
                          unchanged).
landscape_degenerate.csv  The frequency-degenerate member of each
                          family (lambda_s = lambda_i = 2 lambda_p).
landscape_coverage.csv    Per-interaction, per-target record of what
                          was found, including the negative results.
landscape_field.csv       theta_PM field over the phase-matching
                          surfaces (18 Sep 2026 run; the eight LPR2010
                          sequences only -- not regenerated on 17 Oct).

Operating points (examples/landscape_points.py)
-----------------------------------------------
points2_summary.csv       One row per selected operating point
                          (18 Sep 2026 run on the eight LPR2010
                          sequences; unchanged on 17 Oct 2026 -- the
                          new-half points G, H, J, K, M are in the
                          headline files below).
points2/point<id>_purity.csv       purity and filter transmission over
                                   the (L, pump FWHM) grid
points2/point<id>_convergence.csv  purity vs window scale and grid size
points2/point<id>_jsa.csv          |f| in angular-frequency detuning
points2/point<id>_schmidt.csv      Schmidt spectrum
points2/point<id>_brightness.csv   aperture length and brightness proxy
                                   vs waist

Featured points (examples/landscape_headline.py)
------------------------------------------------
headline_points.csv       Table 2 of the manuscript (A, B, C, D, E,
                          G, H, K) plus F, R (featured until v0.10) and
                          J, M (x-y EPM member at 532 nm; continuation
                          of C at 1150 nm).  G, H, K, J, M lie on the
                          half of the type-II surfaces added on
                          17 Oct 2026 (sequences eoo / oee).
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
                          frequency) at a fixed pixel (1.0 rad/ps for
                          A-F, R; 0.25 rad/ps for G-M, whose narrow-arm
                          phase-matching width is 2 rad/ps at 10 mm);
                          least-squares fit P(W) = P_inf + C/W (P_inf,
                          C_per_sigma, residual, fit_points); the
                          manuscript's adaptive window in sigma and nm;
                          the T = 90 % / 70 % filtered purities at 6 and
                          16 sigma with the filter FWHM in nm.  Windows
                          that hit the Sellmeier limit are flagged
                          clipped_* = 1 and excluded from the fit.
length_extension.csv      Points A, B, C (and G, H at 0.4 rad/ps pixel)
                          at L = 10/15/20/30 mm over pump bandwidths
                          2-16 nm: P on the 6-sigma window and P_inf.
                          20 Oct 2026: 1-nm rows for G and H added
                          (examples/add_1nm_rows.py, 0.25 rad/ps pixel)
                          so that H's 30-mm optimum, 0.968 at 1 nm, is
                          in the file; the bandwidth grid of
                          window_tail.py now starts at 1 nm.

Checks added for the first revision (examples/revision_checks.py, v1.2.0)
-------------------------------------------------------------------------
resolution_check.csv      Folds and ends of the five type-II GVM curves
                          and the two weak x-z eoe features, retraced
                          on a 0.1-nm pump grid at the survey's node
                          counts (320 signal nodes, 1201 angles) and at
                          1280 / 2401.  Folds: 722.6 nm (k_p' = k_e')
                          and 882.4 nm (EPM); curve ends move by at
                          most 4 nm against the 2-nm survey grid; the
                          eoe theta_PM = 45 deg feature spans 776-786 nm.
negative_results_check.csv  x-y ooe, x-y eoe, x-z oee and y-z eeo
                          re-scanned at 1 nm / 640 signal nodes / 2401
                          angles (survey field: 5 nm / 220 / 1201), both
                          roots: nodes tested, nodes phase matched, and
                          the theta_PM range, mapped to (-90, 90] deg.
xy_altfit_curves.csv      The two x-y type-II GVM curves traced with the
                          Miyata-Umemura-Kato fit (the survey's), the
                          Hellwig (2000) fit and the Umemura (2007) fit.
                          lambda_e = in-plane daughter, lambda_o =
                          z-polarized daughter.  inside_fit_range = 1
                          when pump and idler lie inside the fit's own
                          stated range (the fits are evaluated outside
                          it where the curve requires).
xy_altfit_H.csv           The member of the x-y k_p' = k_o' curve pumped
                          at 464 nm under each fit, with its purity
                          (10 mm, 2 nm) on 6- and 12-sigma windows.
window_prescription.csv   Adaptive window against the prescription
                          R = 10 x (mean FWHM of the signal and idler
                          central peaks) of Cai, Baek, Jin & Kaneda,
                          Adv. Phys. Res. 4, e00167 (2025), at A, B, H
                          and D; FWHMs from cuts through the JSI peak
                          and from the marginal spectra.
window_tail_dense.csv     P(W) at A, B, C, D, G, H on every window
                          half-width listed in revision_checks.DENSE
                          that the Sellmeier range allows.
window_tail_fits.csv      (P_inf, C) of P = P_inf + C/W fitted over
                          sub-ranges of those windows, including the
                          three windows {6, 8, 12} sigma on which the
                          B and C entries of window_tail.csv rest.
window_tail_conditions.csv  Tail coefficient of the EPM point D at
                          10 mm (2, 2.4, 3 nm) against its tabulated
                          2 mm, 12 nm.  Recomputed here on five windows;
                          agrees with window_tail.csv to about 1e-4 in
                          P and 1e-3 in C (frequency grid and windows
                          fitted; compare window_tail_fits.csv).
deff_benchmarks.csv       |d_eff| at A, D, C with the Kleinman average
                          (3.2 pm/V) and with d26 = 3.48 pm/V, against
                          Jin et al. (2020); the last columns are a
                          reconstruction (d26 cos(theta - rho) x Miller
                          factor), an inference and not a statement of
                          that paper.

headline_points.csv and the other headline_*.csv are unchanged from
v1.1.0; examples/landscape_headline.py now carries the traced inputs of
G, H, J, K, M at full precision and regenerates them byte for byte
(v1.1.0 listed rounded inputs, which gave H on a 5.90- instead of a
6.15-sigma window: 0.95152 instead of the tabulated 0.95124).

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
    python examples/full_landscape_preview.py     # merged curves + preview
    python examples/revision_checks.py            # the revision checks above
    python -m pytest tests                        # 114 tests
