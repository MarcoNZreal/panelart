# Simulation studies

Scripts and results behind the statements about simulation behaviour in the package documentation.
This folder is not part of the built package (see `.Rbuildignore`).

## Layout

* `*.R` in this folder: the study scripts, grouped below.
* `results/`: the CSV tables the scripts wrote, and `results/simulation_results.md`, which collects
  the tables with short readings (numbered sections).

## Running a script

Install the package, set the working directory to this folder and run, for example:

```r
source("kink_calib2.R")                       # cheap: exact finite-sample variance ratios
Sys.setenv(RHO = "0.5", PHI = "0.7", NONSEP = "1", NREP = "20", MODES = "pc")
source("kink_cross1.R")                       # configured through environment variables
```

Scripts write their CSV tables into the working directory; the copies in `results/` are the
tables as produced when the studies were run. Several scripts read their settings from
environment variables (`NREP`, `PHI`, `RHO`, `MODES`, ...): read the top of each script.

**Provenance caveats.**

* The studies were run one at a time while the package was being developed. The interface and
  some defaults changed in between, so a script may need small changes to run against the current
  version, and re-running will not reproduce a table exactly where a default changed afterwards.
  `results/simulation_results.md` says in the relevant sections which setting a table used.
* Replication counts are small (typically 40 to 400) because the studies were run in a single-core
  sandbox with a time limit per run. Treat them as smoke-test-scale results with sampling error of
  a few percentage points, or larger for power at 40 to 60 replications.
* Scripts were syntax-checked; only a few of the cheap ones were re-run on the final version.

## Scripts by topic

| Topic | Scripts |
| --- | --- |
| Mean breaks (size, power, AR(1), bandwidth, localization, scaling, heterogeneity, factors) | `run_sims.R`, `run_sims2.R` (and `inst/scripts/simulation_study.R` in the package) |
| Trend breaks and kink rates and intervals | `trend_study.R`, `trend_study2.R`, `kink_rate.R`, `kink_ci.R`, `kink_ci2.R`, `kink_perseries.R` |
| Kink search, thinning, extra knots, thresholds | `kink_search.R`, `kink_thin.R`, `kink_thin2.R`, `kink_thin3.R`, `kink_theory_checks.R`, `kink_trap.R`, `kink_extras.R`, `kink_extras2.R`, `kink_alphathin.R`, `kink_power2N.R`, `kink_power2N_b.R`, `kink_fixedalpha.R`, `kink_pair1.R`, `kink_pair2.R`, `kink_pair3.R`, `kink_pair_diag.R` |
| Serial dependence | `kink_lemma48.R`, `kink_dep1.R` to `kink_dep7.R`, `kink_dep_diag.R`, `lrv_ar_proto.R`, `lrv_ar_test.R`, `lrv_ar_test2.R` |
| Cross-sectional dependence | `kink_cross1.R`, `kink_cross2.R`, `kink_cross_oracle.R`, `lrv_pc_proto.R` |
| Calibration of the penalty | `kink_calib1.R`, `kink_calib2.R`, `kink_calib3.R` |

## Summary of what was checked

`tests/test-panelart.R` runs with `R CMD check` and passes. It checks the criterion against brute
force (difference about 3e-14), simple break recovery, no false positives on a null panel, the
long-run covariance, and that each main function runs.

Monte Carlo studies (sandbox, 50 to 300 replications per cell; treat them as smoke-test-scale
results with sampling error of a few percentage points) are in `simulation_results.md`, with the
CSV tables and the scripts `run_sims.R`, `run_sims2.R`. Headlines:

* **Limit law.** Simulated quantiles of |xi| (1.52, 7.68, 11.00 at 50, 90, 95%) match the
  closed-form density (1.505, 7.687, 11.033).
* **Independent errors.** The Monte Carlo penalty `"mc"` has 4 to 6% false positives on null data and
  more power than the union-bound penalty; pooling gain is large (N = 10, jump 0.4: exact break
  count 94% against 3% for classic BIC).
* **AR(1) errors.** The iid-calibrated penalty fails from phi = 0.3. `"hac_mc"` keeps size
  (at most 2% on the null) and at phi = 0.6 finds the exact count in 78% of runs against 64% for
  classic BIC, which has 16% false positives there. At phi = 0.8 nothing has power at this signal.
* **Bandwidth.** Too short a bandwidth gives false positives; too long a bandwidth, applied to data that
  still contain the breaks, inflates the threshold above the signal. The first pass of the iterated
  scheme therefore uses a short bandwidth, and the automatic rule is used from the second pass on.
* **Scaling.** At fixed pooled signal, power and the localization error as a fraction of T do not
  depend on T. At the same pooled signal localization is worse at larger N because the per-series
  signal is smaller, as the notes predict.
* **Strong factor.** For jumps orthogonal to the factors, factor projection finds the exact
  count in 95 to 97% of runs against 13 to 70% for `"hac_mc"` on the raw panel. For jumps inside
  the factor space every method is weak.
* **Trend breaks** (`type = "trend"`, T = 200, N = 10, 60 replications; `study_trend.csv`). With a
  common trend and no break, mean-based methods find breaks in every run (about 3 on average for
  `art_fit`), while the trend version has a 7% false-positive rate. For a slope break of 0.02 or 0.04
  per series it finds exactly one break in 95% of runs, against 3 to 5 steps for the mean-based
  methods. Dates of a kink are less precise than dates of a jump (mean error 7 to 12 dates, with
  25 to 43% of dates within 5 of the truth). For a pure level jump it still detects the break (97%)
  with somewhat lower precision than the mean version.
* **Kinks** (`art_kink`, same design plus two kinks; `study_trend2.csv`, `study_kink_rate.csv`). For
  slope changes the join-point estimator puts 73 to 100% of dates within 5 of the truth against 28
  to 52% for `type = "trend"`, with mean errors of 2 to 4 dates against 6 to 12, but it mostly
  misses level jumps (30% detection). The date error as a fraction of T falls like T^-0.57 for
  join points and T^-0.36 for separate lines in a 30-replication check (heuristic predictions
  -1/2 and -1/3). Classic ART with line fits finds one break in the two-kink design, in the wrong
  place. Theory notes: section 10 of the notes document.
* **Kink intervals** (`art_ci` on `art_kink` fits; `study_kink_ci.csv`, `study_kink_perseries.csv`,
  `study_kink_ci2.csv`). The predicted standard deviation of the knot date, 2 sigma / (norm(gamma)
  sqrt(T tau (1 - tau))), matches simulation (2.24 against 2.38, 4.47 against 5.12, 2.67 against
  3.04, 4.22 against 4.61 dates), and a variance inflation factor 1 + N sigma^2 / (norm(gamma)^2
  norm(h)^2) matches the robust spread when the per-series signal is small. Coverage of 95%
  intervals is 85 to 94% (heavy tails from occasional mislocated knots), against 60% at N = 500
  without the inflation factor.
* **Kink search** (`kink_search.R`, `study_kink_repartition.csv`, `study_kink_tent.csv`). Repartition from
  starts offset by up to 40 dates (13% of T) reaches the same knots every time. Offsetting kinks (a
  tent) are no harder than same-direction ones at spacing T/3. The search adds an occasional
  duplicate knot next to a true one at high signal (exact count 73% at slope 0.04 against 90% at
  0.02); the thinning step described in the next bullet was added after this study.
* **Kink thinning** (`kink_thin2.R`, `kink_thin3.R`). In two-kink designs thinning raises the exact-count
  rate from 0.835 to 0.892 (tent, slope 0.04) and from 0.875 to 0.905 (same direction, slope 0.04) at
  400 replications, and lowers the share of runs with a duplicate knot from 4 to 9% to 3 to 4%.
  One-kink fits are unchanged. Remaining extra knots are mostly ordinary false positives.
* **Penalty size for kinks** (`kink_theory_checks.R`, `study_kink_union_vs_mc.csv`). With the larger union-bound
  penalty the search had no false positives and no far knots at slope 0.04 (exact count 1.00), but at
  slope 0.03 it ended with one knot about 50 dates from both true knots in 54% of runs, because a
  second knot adds less than the penalty. The Monte Carlo penalty has about 5 to 6% false positives
  and finds both knots. The removal cost of a knot equals norm(gamma_hat)^2 times norm(h)^2 (checked to
  1e-14), and the sufficient signal in the covering theorem is about 3000 times the simulated one.
* **Extra knots after thinning** (`kink_extras.R`, `kink_extras2.R`). With the union-bound penalty no extra
  knot survived in 900 runs across three designs (two or three kinks) and no true knot was missed;
  with the Monte Carlo penalty extras survive at about its 5 to 7% false-positive level.
* **Thinning threshold** (`kink_alphathin.R`, `kink_power2N_b.R`, `kink_fixedalpha.R`). `art_kink(alpha_thin = "2N")`
  acts as a stricter overall level for extra knots: in 200 two-kink panels it removed all 10 isolated false-positive
  knots and all 8 duplicate pairs left by thinning at the penalty (182 to 200 exact fits). A fixed penalty at the
  same level gives the same 200 exact fits. Both lose weak true knots next to a strong one (at population gain 33.6
  the weak knot was found in 9% of runs with 2N thinning and 2% with the fixed penalty, against 66% with
  thinning at the penalty), so it is an option and not the default.
* **Pair test** (`kink_pair_diag.R`, `kink_pair1.R`, `kink_pair2.R`). The experimental `pair_test = TRUE` merged
  11 of 19 duplicate pairs, the same as thinning at the penalty, and no more; it costs no power on weak knots.
* **Serially correlated errors for kinks** (`kink_dep1.R` to `kink_dep5.R`, `kink_lemma48.R`).
  `art_kink(penalty = "hac_mc")` multiplies the Monte Carlo penalty by an iterated long-run variance estimate
  and `art_ci` then uses it. Under AR(1) errors at T = 300, N = 10 the i.i.d.-calibrated penalty finds knots in
  98 to 100% of null panels, `hac_mc` in 5 to 6% (phi = 0.3, 0.5, 0.7). Single-kink localization sd matches the
  long-run-variance prediction (1.59 against 1.58 at phi = 0.5) and 95% intervals cover 93 to 95%, against 45 to 77%
  with the residual variance. Limits: power falls at phi = 0.7 (both kinks found 62% at slope 0.04); the pooled
  version over-rejects when series differ in persistence (11% at phi ~ U(0, 0.8); `scale = "series"` gives 7.5%
  but is worse at phi ~ U(0.5, 0.9)); at phi >= 0.8 the long-run variance is underestimated and size fails
  (21% at phi = 0.8, 83% at phi = 0.9 for T = 300).
  The iteration can settle at too few knots when two same-sign kinks look like one bend (7 of 40 runs at
  phi = 0.7, slope 0.06). `lrv_extra = 2` removes those misses but raises the extra-knot rate when a real
  kink is present (single-kink runs with exactly one knot 0.975 to 0.855 at phi = 0.7), so it is off by default.
* **Knot-free first pass** (`lrv_ar_proto.R`, `lrv_ar_test.R`, `lrv_ar_test2.R`, `kink_dep7.R`). The default
  `art_kink(penalty = "hac_mc", lrv = "ar_start")` takes its first long-run variance from AR models fitted to the
  autocovariances of the first differences, which the kinks barely affect, and uses the residual estimate
  afterwards. It removes the same-sign trap (both kinks found 0.975 against 0.825 at phi = 0.7, slope 0.06;
  0.90 against 0.275 at phi = 0.8; 1.0 against 0.1 for AR(2) errors) and the high-persistence size failure
  (null false positives 3% and 0% at phi = 0.8 and 0.9 against 21% and 83%). It matches `hac` on every null case
  tried, including ARMA and MA errors. `lrv = "ar"` (knot-free in every pass) is not recommended: it
  over-rejects under ARMA/MA errors (14%, 16%, 58%).
* **Cross-sectional dependence** (`kink_cross1.R`, `kink_cross2.R`). `art_kink(penalty = "hac_mc", cross = ...)`:
  `"separable"` (eigenvalue-weighted penalty from the first-difference covariance, sandwich intervals through
  `fit$Omega`), `"whiten"` (whitened panel), `"full"` (N x N long-run covariance from residuals).
  With equicorrelation 0.5 and N = 10, ignoring the dependence gives 34% null false positives; `separable` and
  `whiten` give 5%. Jumps along the common factor are located about 2.3 times less precisely than jumps
  orthogonal to it, which the sandwich variance reproduces (2.14 against 2.14 empirical for aligned jumps);
  whitening finds weak orthogonal jumps (93% against 18% at slope 0.015) but loses aligned ones (62.5% against 95%).
  When the common factor is more persistent than the idiosyncratic part, `separable` fails (31% false positives)
  and `full` improves it (13%) without reaching nominal size.
* **Principal-direction-wise long-run covariance** (`lrv_pc_proto.R`, `kink_cross_oracle.R`). `cross = "pc"` finds the
  factor directions from the first-difference covariance, estimates each factor direction's long-run variance
  with the knot-free AR estimate of section 19, pools the rest, and floors each factor variance at the separable
  estimate (`pc_floor = TRUE`). Null false positives: 5% when the common factor is more persistent than the
  noise (13% `full`, 31% `separable`), 2% when the dependence is separable, 5% for uncorrelated series. It is
  conservative: power at slope 0.05 is 50% (aligned) and 40% (orthogonal), against 70% and 52.5% for an oracle
  with the true covariance (itself conservative, 0% null false positives) and 80% and 70% for `full`, whose
  higher power comes from under-penalizing. Shrinking the factor variances is not a remedy (x 0.5 gives 25%).
* **Calibration** (`kink_calib1.R`, `kink_calib2.R`, `kink_calib3.R`). The true 95% quantile of the null maximum under
  dependent errors is below the i.i.d.-in-time Monte Carlo quantile for the same long-run covariance, by 5% to 25%
  for AR(1) 0.3 to 0.9 at T = 300 and 14% for the factor design; the shortfall is the finite-sample variance of a
  step-like increment (Lemma 48), computed exactly from the autocovariances. The knot-free AR estimate now returns
  that effective variance (omega^2 times the step ratio) instead of omega^2, which moves `ar_start`/`ar` null rates
  to nominal at phi <= 0.8 (`ar` is liberal at 0.9: 13%). Residual-based estimates carry the same downward bias, so
  the two cancel. Use `nsim >= 2000` with the weighted penalty (`cross != "none"`): at 500 draws its simulation
  error is 6% (sd) for factor-type eigenvalues.
* **First-pass overshoot** (`pc_lrv = "resid"`, now the default for `cross = "pc"`, with `lrv_extra` defaulting to 2).
  The factor variances come from the knot-free estimate in the first pass only; later passes re-estimate them from
  residuals after the fitted knots plus two unpenalized extra knots. In the non-separable factor design (slope
  0.05) both kinks are found in 75% of aligned and 78% of orthogonal runs, against 57.5% and 42.5% with the AR
  estimate in every pass (oracle at the true quantile: 83% and 92%), with null false positives of 3% and no
  loss in coverage; a second design (factor AR(1) 0.5) gives 92% against 80%. Evidence is from two designs with
  N = 10, T = 300 and one factor.
* **Intervals.** 95% limit-law intervals cover 94 to 98% in four of five designs and 91% in the
  worst (N = 50, jump 0.2).

Not yet done: more than one factor or weak factors, non-Gaussian errors, comparison with the R
packages `tree` and `strucchange`, and any real data.
