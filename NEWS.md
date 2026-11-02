# panelart 0.1.0 (development version)

First version. Research software accompanying a paper in preparation; the
interface and some defaults may still change before a CRAN release.

## Mean breaks

* `art_classic()`: classic greedy regression tree on time (BIC, cost-complexity
  or no pruning); `type = "trend"` fits a line in each node.
* `art_fit()`: extended estimator with a forward-repartition-eliminate search and
  hull-merge thinning. Penalties `"iid"` (union bound), `"mc"` (Monte Carlo),
  `"hac"` and `"hac_mc"` (iterated long-run covariance) and `"fixed"`; whitening
  with `whiten = "diff"`, `"lrv"` or a known `Omega`; `type = "trend"` fits a
  separate line in each segment.
* `art_fit_factor()`: factor projection (whitened factor pipeline plus an
  orthogonal pipeline), with `art_nfactors()` for the number of factors.
* `art_ci()`: limit-law confidence intervals for break dates
  (`art_xi_density()`, `art_xi_quantile()`).

## Join points (continuous piecewise-linear breaks)

* `art_kink()`: join-point estimator with common knots, search by forward steps
  and repartition, and thinning by re-optimized removal costs.
* Serial dependence: `penalty = "hac_mc"` with an iterated long-run variance;
  a knot-free first pass (`lrv = "ar_start"`, the default for `"hac_mc"`) that
  uses first differences; options `scale`, `lrv_extra`.
* Cross-sectional dependence: `cross = "separable"`, `"whiten"`, `"full"` and
  `"pc"` (principal-direction-wise long-run covariance, the most recent and least
  tested), with `art_penalty_kink_mc_w()` for the weighted Monte Carlo penalty.
  `art_ci()` uses a sandwich variance when a long-run covariance is available.
* Experimental, off by default: `alpha_thin = "2N"`, `pair_test = TRUE`.

## Tools

* `art_lrcov()`, `art_penalty*()`, `art_whiten_matrix()`, `art_bandwidth_check()`.
* `art_sim_data()`, `art_sim_study()`, `art_sim_summary()`, `art_eval()` and
  `art_methods()` for Monte Carlo studies; `inst/scripts/simulation_study.R`.

## Known limitations

* Evidence is from simulations with Gaussian errors, mostly N = 10 and T = 300;
  there is no real-data application yet.
* Not yet compared with the R packages `tree` and `strucchange`.
* Defaults for the dependence options were chosen on one or two designs.
