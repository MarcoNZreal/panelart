# panelart

<!-- badges: start -->
[![R-CMD-check](https://github.com/YOUR-GITHUB-USERNAME/panelart/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YOUR-GITHUB-USERNAME/panelart/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package for **Atheoretical Regression Trees (ART)** applied to common breaks in a panel of
time series: breaks in the mean, and continuous piecewise-linear (join-point) breaks, with
penalties that allow for serial and cross-sectional dependence in the errors. It extends the ART
method of Rea, Reale, Cappelli and Brown (2010) from one series to a panel and adds a search
algorithm, confidence intervals and tools for Monte Carlo studies.

> **Status: pre-release research software.**
>
> * Not on CRAN. The interface and several defaults changed during development and may change again.
> * The evidence is simulation only: Gaussian errors, mostly N = 10 series and T = 300 dates.
>   There is no real-data application yet, and no comparison yet with the R packages `tree` and
>   `strucchange`.
> * Several help pages refer to "the notes": theory notes with proofs, each result carrying a
>   status label (proved, by inspection, heuristic, unproven). They accompany a paper in
>   preparation and are not part of this repository yet.
> * The dependence options for join points (`penalty = "hac_mc"`, `cross = ...`) are the least
>   tested part. Check their size in your own design before relying on them.

Data convention: `Y` is a `T x N` matrix (rows = dates, columns = series). A break date `d` means
the mean changes between row `d` and row `d + 1`. For join points, `d` is the knot.

## Install

```r
# install.packages("remotes")
remotes::install_github("YOUR-GITHUB-USERNAME/panelart")
```

Only base R packages (`stats`, `graphics`, `grDevices`, `parallel`) are needed.

## Quick start

```r
library(panelart)

## Breaks in the mean
set.seed(1)
sim <- art_sim_data(T = 200, N = 10, breaks = c(0.4, 0.7), delta_bar = 0.6)
fit <- art_fit(sim$Y, penalty = "mc")   # independent errors, calibrated penalty
fit
summary(fit)                            # segments, penalty, limit-law intervals
plot(fit)
art_classic(sim$Y, pruning = "bic")$dates   # classic ART for comparison

## Join points (continuous piecewise-linear breaks)
set.seed(1)
simk <- art_sim_data(T = 200, N = 10, breaks = 0.5, delta_bar = 0, slope_bar = 0.03)
fk <- art_kink(simk$Y)                  # independent errors
fk$dates
art_ci(fk)                              # normal-theory intervals for the knots

## Serially correlated errors
set.seed(2)
simd <- art_sim_data(T = 300, N = 10, breaks = 0.5, delta_bar = 0, slope_bar = 0.04,
                     errors = "ar1", phi = 0.5)
art_kink(simd$Y, penalty = "hac_mc")$dates
```

## Which option for which errors

| Errors | Mean breaks (`art_fit`) | Join points (`art_kink`) |
| --- | --- | --- |
| Independent over time and across series | `penalty = "mc"` | `penalty = "mc"` (default) |
| Serially correlated | `penalty = "hac_mc"` | `penalty = "hac_mc"`. In simulations the null size held for AR(1) errors up to 0.9 at T = 300, but not for an ARMA(1,1) case with 0.8 and -0.4 (31% false positives); power falls as persistence rises |
| Correlated across series | `whiten`, or `art_fit_factor()` | `cross = "separable"` or `"whiten"` |
| A common factor more persistent than the noise | `art_fit_factor()` | `cross = "pc"` |

These recommendations rest on the simulations summarized in `simulations/README.md`; they are
not guarantees. Use `nsim = 2000` (the default) or more for the Monte Carlo penalties.

## Main functions

| Function | Purpose |
| --- | --- |
| `art_fit()` | Extended estimator for common mean breaks. Penalties: `"iid"` (union bound, as in the proofs), `"mc"` (Monte Carlo calibrated, independent errors), `"hac"` / `"hac_mc"` (long-run covariance, iterated), `"fixed"`. Options for whitening (`whiten = "diff"`, `"lrv"`, or a known `Omega`). |
| `art_fit(type = "trend")` | Breaks in intercept and/or slope: a separate line is fitted in each segment for every series (discontinuous piecewise-linear fit). Checked by simulation only. No confidence intervals. Not available in `art_fit_factor()`. |
| `art_kink()` | Continuous piecewise-linear (join-point) breaks with common knots. Penalties `"mc"`, `"iid"`, `"fixed"`, `"hac_mc"`; dependence options `lrv`, `scale`, `lrv_extra`, `cross`; thinning of duplicate knots (`thin`, on by default). Much better dates for pure kinks than `type = "trend"`; it mostly misses level jumps. |
| `art_fit_factor()` | Factor projection: whitened factor pipeline plus orthogonal pipeline, dates merged. |
| `art_classic()` | Classic greedy regression tree on time, with BIC / cost-complexity / no pruning. `type = "trend"` fits a line in each node. |
| `art_ci()` | Confidence intervals for break dates: limit law for mean breaks, Gauss-Newton or sandwich variance for knots. |
| `art_lrcov()` | Long-run covariance from segment residuals (Parzen or Bartlett kernel, Andrews bandwidth). |
| `art_penalty*()` | Penalties: `art_penalty`, `art_penalty_iid`, `art_penalty_mc`, `art_penalty_mc_cov`, `art_penalty_kink_mc`, `art_penalty_kink_mc_w`. |
| `art_sim_data()` | Simulate panels: iid, AR(1), equicorrelated or factor errors; jumps random, inside or orthogonal to the factor space. |
| `art_sim_study()`, `art_sim_summary()`, `art_eval()`, `art_methods()` | Monte Carlo study over a grid of scenarios and methods. |

`inst/scripts/simulation_study.R` contains five ready-to-run studies (size and power, AR(1),
strong factor, interval coverage, bandwidth window).

## What is implemented, and how it differs from the theory

* The criterion is the pooled reduction in the sum of squares, computed from cumulative sums.
* **Penalty constants.** The proofs need `C0 >= 100` in the union-bound penalty. That is very
  conservative: for `N = 1`, `T = 200` the union-bound penalty is about `63 sigma^2` against about
  `10 sigma^2` for a Monte Carlo null quantile. `"iid"` and `"hac"` implement the proof form (with
  `C0 = 10` by default, not 100). `"mc"` and `"hac_mc"` calibrate by simulation and are the
  practical defaults for simulation work. Check their size in your design.
* After whitening the penalty is the iid form for `"iid"`/`"hac"` and the Monte Carlo quantile
  with variance one for `"mc"`/`"hac_mc"`. The union-bound form is very conservative after
  whitening.
* **`"hac_mc"`** draws null panels from `N(0, Omega_hat)` independently over time. This captures
  the long-run covariance, not the finite-sample serial structure.
* **Long-run covariance.** Parzen kernel by default (the notes use Bartlett), Andrews AR(1)
  plug-in bandwidth, first-order demeaning bias correction. It is still biased downward in finite
  samples (about 7–9% for AR(1) with `phi = 0.5`, `T = 300`).
* For `N > 200` the covariance is approximated by its diagonal (a warning-free fallback;
  cross-sectional dependence is then ignored).
* `art_fit_factor()` estimates loadings by principal components from the same data and does not
  model the dependence this creates. `r` is chosen by the Bai–Ng IC_p2 criterion plus `extra = 1`.
* `art_classic()` mimics the R package `tree` stopping rules (`minsize`, `mincut`, `mindev`) but is
  not a reimplementation of it. Its BIC rule is a BIC-type approximation.
* The theory is for Gaussian errors, a fixed number of breaks, and regimes of length at least
  `lambda * T`. See the notes for the assumptions behind each result.

## Testing

`tests/test-panelart.R` runs with `R CMD check` (about 15 seconds). It checks the criterion against
brute force, break recovery, size on null panels, the long-run covariance estimates, the
dependence options and a number of regression cases found during development. GitHub Actions runs
`R CMD check` on Linux, macOS and Windows.

## Citation

Please cite the package (`citation("panelart")`) and the original ART method:

Rea, W., Reale, M., Cappelli, C. and Brown, J. (2010). Identification of changes in mean with
regression trees: an application to market research. *Econometric Reviews* 29(5-6), 754-777.
doi:10.1080/07474938.2010.482001

A reference for the paper accompanying this package will be added when it is available.

## License

GPL (>= 2). See `LICENSE.md`.

## Issues

Bug reports and questions are welcome through the issue tracker.
