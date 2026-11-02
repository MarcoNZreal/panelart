# panelart simulation results (sandbox runs)

All runs: panelart 0.1.0 (current source), one core, Gaussian errors, sigma = 1. Replications are 50 to 300 per cell, so sampling error is about 2 to 7 percentage points on proportions. These are smoke-test-scale results, not a full study. Scripts: `run_sims.R` (studies 1 to 5) and `run_sims2.R` (studies 6 to 8). Files ending `_before_fix.csv` hold runs made before the first-pass bandwidth fix described below.

Definitions: *exact* = estimated number of breaks equals the true number; *spurious* = at least one estimated date farther than max(5, 0.02 T) dates from every true break (this also counts poorly localized weak breaks); *found* = fraction of true breaks with an estimate within that tolerance. Methods: `classic_bic` = classic ART with BIC-type pruning; `panel_iid` = proof-form union-bound penalty (C0 = 10); `panel_mc` = Monte Carlo penalty for independent errors; `panel_hac` = iterated long-run covariance with Monte Carlo penalty; `panel_factor` = factor projection with whitened factor series; `panel_factor_nowhite` = the same without whitening.

**Bandwidth fix.** In the iterated long-run covariance scheme, the first pass leaves the breaks in the data. An automatic bandwidth computed there reads the steps as strong autocorrelation and picks a very long window, which inflates the threshold above the signal. The first pass now uses floor(T^(1/3)) and the automatic rule is used from the second pass on. Studies 2 and 8 below were rerun after the fix. Studies 1, 3, 4, 5, 6, 7 are unaffected (they do not use the automatic bandwidth in a way the fix changes; study 7 `panel_hac` results were identical before and after).

## 1. Independent errors, T = 200, two breaks (null when jump = 0), 200 replications

| N | jump per series | method | breaks found (mean) | exact | spurious | found | mean abs error (dates) |
|---|---|---|---|---|---|---|---|
| 1 | 0 | classic_bic | 0.05 | 0.95 | 0.05 | – | – |
| 1 | 0 | panel_iid | 0.00 | 1.00 | 0.00 | – | – |
| 1 | 0 | panel_mc | 0.07 | 0.94 | 0.07 | – | – |
| 1 | 0.2 | classic_bic | 0.27 | 0.02 | 0.20 | 0.02 | 40.1 |
| 1 | 0.2 | panel_iid | 0.00 | 0.00 | 0.00 | 0.00 | – |
| 1 | 0.2 | panel_mc | 0.28 | 0.01 | 0.23 | 0.02 | 42.8 |
| 1 | 0.4 | classic_bic | 0.64 | 0.08 | 0.37 | 0.12 | 34.8 |
| 1 | 0.4 | panel_iid | 0.00 | 0.00 | 0.00 | 0.00 | – |
| 1 | 0.4 | panel_mc | 0.69 | 0.09 | 0.41 | 0.12 | 35.8 |
| 1 | 0.8 | classic_bic | 1.89 | 0.79 | 0.57 | 0.57 | 10.4 |
| 1 | 0.8 | panel_iid | 0.34 | 0.00 | 0.13 | 0.11 | 33.8 |
| 1 | 0.8 | panel_mc | 1.73 | 0.77 | 0.46 | 0.59 | 8.7 |
| 10 | 0 | classic_bic | 0.00 | 1.00 | 0.00 | – | – |
| 10 | 0 | panel_iid | 0.00 | 1.00 | 0.00 | – | – |
| 10 | 0 | panel_mc | 0.06 | 0.95 | 0.05 | – | – |
| 10 | 0.2 | classic_bic | 0.00 | 0.00 | 0.00 | 0.00 | – |
| 10 | 0.2 | panel_iid | 0.00 | 0.00 | 0.00 | 0.00 | – |
| 10 | 0.2 | panel_mc | 1.08 | 0.20 | 0.49 | 0.26 | 31.8 |
| 10 | 0.4 | classic_bic | 0.78 | 0.03 | 0.14 | 0.32 | 32.5 |
| 10 | 0.4 | panel_iid | 0.68 | 0.01 | 0.12 | 0.28 | 33.4 |
| 10 | 0.4 | panel_mc | 2.04 | 0.94 | 0.24 | 0.88 | 2.3 |
| 10 | 0.8 | classic_bic | 2.00 | 0.99 | 0.02 | 0.99 | 0.5 |
| 10 | 0.8 | panel_iid | 2.00 | 1.00 | 0.00 | 1.00 | 0.3 |
| 10 | 0.8 | panel_mc | 2.10 | 0.91 | 0.10 | 1.00 | 0.4 |
| 50 | 0 | classic_bic | 0.00 | 1.00 | 0.00 | – | – |
| 50 | 0 | panel_iid | 0.00 | 1.00 | 0.00 | – | – |
| 50 | 0 | panel_mc | 0.06 | 0.95 | 0.04 | – | – |
| 50 | 0.2 | classic_bic | 0.00 | 0.00 | 0.00 | 0.00 | – |
| 50 | 0.2 | panel_iid | 0.56 | 0.00 | 0.10 | 0.23 | 33.5 |
| 50 | 0.2 | panel_mc | 1.98 | 0.85 | 0.34 | 0.78 | 6.0 |
| 50 | 0.4 | classic_bic | 0.59 | 0.00 | 0.01 | 0.29 | 33.1 |
| 50 | 0.4 | panel_iid | 2.00 | 1.00 | 0.00 | 1.00 | 0.3 |
| 50 | 0.4 | panel_mc | 2.02 | 0.97 | 0.03 | 1.00 | 0.3 |
| 50 | 0.8 | classic_bic | 2.00 | 1.00 | 0.00 | 1.00 | 0.0 |
| 50 | 0.8 | panel_iid | 2.00 | 1.00 | 0.00 | 1.00 | 0.0 |
| 50 | 0.8 | panel_mc | 2.06 | 0.94 | 0.06 | 1.00 | 0.0 |

Reading: with jump 0 the rows show false-positive rates; `panel_mc` is 4 to 6% (nominal 5% per segment). Pooling gain: at N = 10, jump 0.4, `panel_mc` has the exact count in 94% of runs against 3% for classic BIC. At N = 50, jump 0.2 the figures are 86% and 0%. The union-bound penalty is safe but loses power at small signal.

## 2. AR(1) errors, T = 300, N = 10, two breaks, 50 replications (after the bandwidth fix)

| phi | jump | method | breaks found (mean) | exact | spurious | found | exact before fix |
|---|---|---|---|---|---|---|---|
| 0 | 0 | classic_bic | 0.00 | 1.00 | 0.00 | – | 1.00 |
| 0 | 0 | panel_hac | 0.04 | 0.96 | 0.04 | – | 0.95 |
| 0 | 0 | panel_mc | 0.06 | 0.94 | 0.06 | – | 0.94 |
| 0 | 0.7 | classic_bic | 2.00 | 1.00 | 0.00 | 1.00 | 1.00 |
| 0 | 0.7 | panel_hac | 2.04 | 0.96 | 0.04 | 1.00 | 0.96 |
| 0 | 0.7 | panel_mc | 2.06 | 0.94 | 0.06 | 1.00 | 0.95 |
| 0.3 | 0 | classic_bic | 0.00 | 1.00 | 0.00 | – | 1.00 |
| 0.3 | 0 | panel_hac | 0.02 | 0.98 | 0.02 | – | 1.00 |
| 0.3 | 0 | panel_mc | 4.44 | 0.02 | 0.98 | – | 0.01 |
| 0.3 | 0.7 | classic_bic | 2.00 | 1.00 | 0.12 | 0.94 | 1.00 |
| 0.3 | 0.7 | panel_hac | 2.00 | 1.00 | 0.10 | 0.95 | 0.99 |
| 0.3 | 0.7 | panel_mc | 5.42 | 0.00 | 1.00 | 0.97 | 0.00 |
| 0.6 | 0 | classic_bic | 0.20 | 0.84 | 0.16 | – | 0.85 |
| 0.6 | 0 | panel_hac | 0.00 | 1.00 | 0.00 | – | 1.00 |
| 0.6 | 0 | panel_mc | 1.02 | 0.00 | 1.00 | – | 0.00 |
| 0.6 | 0.7 | classic_bic | 2.50 | 0.64 | 0.60 | 0.80 | 0.69 |
| 0.6 | 0.7 | panel_hac | 1.78 | 0.78 | 0.28 | 0.74 | 0.30 |
| 0.6 | 0.7 | panel_mc | 1.10 | 0.10 | 0.32 | 0.39 | 0.09 |
| 0.8 | 0 | classic_bic | 14.74 | 0.00 | 1.00 | – | 0.00 |
| 0.8 | 0 | panel_hac | 0.04 | 0.98 | 0.02 | – | 1.00 |
| 0.8 | 0 | panel_mc | 1.00 | 0.00 | 1.00 | – | 0.00 |
| 0.8 | 0.7 | classic_bic | 14.50 | 0.00 | 1.00 | 0.81 | 0.00 |
| 0.8 | 0.7 | panel_hac | 0.34 | 0.02 | 0.18 | 0.08 | 0.00 |
| 0.8 | 0.7 | panel_mc | 1.00 | 0.00 | 0.60 | 0.20 | 0.00 |

Reading: the iid-calibrated `panel_mc` fails as soon as phi is 0.3 (false positives in about 99% of null runs). `panel_hac` keeps size at every phi (at most 2% false positives on the null). After the fix its power at phi = 0.6 is an exact count in 78% of runs with breaks (30% before), against 64% for classic BIC, which has 16% false positives on the null at phi = 0.6 and breaks down at phi = 0.8 (about 14 false breaks). At phi = 0.8 `panel_hac` has almost no power at this signal (the long-run variance is 25).

## 3. Bandwidth window, AR(1) phi = 0.5, T = 300, N = 10, `panel_hac` with a fixed bandwidth, 50 replications

| L | P(any break on null) | exact count with two breaks | fraction found |
|---|---|---|---|
| 2 | 0.72 | 0.20 | 0.88 |
| 5 | 0.14 | 0.86 | 0.90 |
| 10 | 0.02 | 0.96 | 0.91 |
| 20 | 0.00 | 0.92 | 0.86 |
| 40 | 0.00 | 0.00 | 0.00 |
| 80 | 0.00 | 0.00 | 0.00 |

Reading: both ends of the window matter. L = 2 underestimates the long-run variance (72% false positives). L of 10 to 20 works. L of 40 or more finds nothing, because breaks left in the raw covariance inflate the first-pass threshold above the signal, which is the condition L log T much smaller than T in the notes. This is the same mechanism as the automatic-bandwidth problem fixed above. A user-supplied L is used in the first pass as given.

## 4. Localization with fixed pooled jump, T = 300, one break, Delta = 3, 200 replications

| N | per-series delta^2 T | found | mean abs error | 90th pct abs error |
|---|---|---|---|---|
| 3 | 300 | 0.96 | 0.91 | 3 |
| 10 | 90 | 0.94 | 0.96 | 3 |
| 30 | 30 | 0.96 | 1.14 | 3 |
| 100 | 9 | 0.96 | 1.73 | 5 |
| 300 | 3 | 0.94 | 1.92 | 5 |

Reading: the limit law predicts a mean absolute error of E|xi| sigma^2 / Delta = 3 x 0.33 = 1.0 date. The observed 0.91 to 1.14 for N = 3 to 30 agrees. For N = 100 and 300, where the per-series signal falls to 9 and 3, the error grows to 1.7 and 1.9, the direction the per-series condition (H) predicts.

## 5. Coverage of 95% limit-law intervals (with noise-bias correction), T = 400, one break, 300 replications

| N | jump per series | Delta | one break found | coverage | mean width (dates) |
|---|---|---|---|---|---|
| 10 | 0.3 | 0.9 | 0.94 | 0.965 | 25.0 |
| 10 | 0.5 | 2.5 | 0.95 | 0.947 | 8.9 |
| 10 | 1 | 10.0 | 0.96 | 0.976 | 2.2 |
| 50 | 0.2 | 2.0 | 0.94 | 0.907 | 11.0 |
| 50 | 0.4 | 8.0 | 0.94 | 0.936 | 2.8 |

Reading: coverage is near 95% in four of five designs. At N = 50, jump 0.2 it is 0.907 (0.879 before the bias correction), about two standard errors low. Coverage is conditional on the fit finding exactly one break.

## 6. Scaling in T and N at fixed pooled signal s = Delta T / sigma^2, `panel_mc`, 100 replications

Two breaks (null when s = 0); per-series jump squared = s / (N T). If the theory is right, power depends on s and N but not on T, and the error as a fraction of T is constant in T.

| s | N | T | exact | spurious | found | mean abs error (dates) | error / T |
|---|---|---|---|---|---|---|---|
| 0 | 5 | 100 | 0.92 | 0.08 | – | – | – |
| 0 | 5 | 400 | 0.98 | 0.02 | – | – | – |
| 0 | 5 | 1000 | 0.95 | 0.05 | – | – | – |
| 0 | 50 | 100 | 0.90 | 0.10 | – | – | – |
| 0 | 50 | 400 | 0.98 | 0.02 | – | – | – |
| 0 | 50 | 1000 | 0.96 | 0.04 | – | – | – |
| 300 | 5 | 100 | 0.98 | 0.13 | 0.94 | 1.2 | 0.0124 |
| 300 | 5 | 400 | 0.95 | 0.30 | 0.83 | 4.6 | 0.0116 |
| 300 | 5 | 1000 | 0.95 | 0.42 | 0.79 | 13.6 | 0.0136 |
| 300 | 50 | 100 | 0.70 | 0.25 | 0.76 | 6.8 | 0.0685 |
| 300 | 50 | 400 | 0.66 | 0.57 | 0.52 | 30.3 | 0.0757 |
| 300 | 50 | 1000 | 0.66 | 0.46 | 0.63 | 62.9 | 0.0629 |
| 600 | 5 | 100 | 0.95 | 0.05 | 1.00 | 0.4 | 0.0039 |
| 600 | 5 | 400 | 0.97 | 0.13 | 0.95 | 2.4 | 0.0060 |
| 600 | 5 | 1000 | 0.96 | 0.10 | 0.97 | 5.6 | 0.0056 |
| 600 | 50 | 100 | 0.92 | 0.12 | 0.96 | 1.0 | 0.0097 |
| 600 | 50 | 400 | 0.96 | 0.24 | 0.89 | 3.2 | 0.0081 |
| 600 | 50 | 1000 | 0.92 | 0.31 | 0.84 | 9.7 | 0.0097 |

Reading: false positives are 2 to 10% on the null at every T and N. At fixed s, the exact-count rate is stable in T (N = 5, s = 300: 0.98, 0.95, 0.95). The localization error as a fraction of T is also stable in T. At the same pooled signal the error is about 6 times larger at N = 50 than at N = 5 for s = 300 (0.07 against 0.012), and about 2 times for s = 600, because the per-series signal is smaller. That is the per-series condition (H).

## 7. Heterogeneous jump sizes, N = 20, T = 300, three breaks at 0.25, 0.5, 0.75, 60 replications

Jump sizes are per-series root-mean-square values for the large, medium and small break. Columns give the fraction of runs in which each break is found.

| errors | sizes | method | large | medium | small | breaks found (mean) | spurious |
|---|---|---|---|---|---|---|---|
| ar1 | 0.5/0.3/0.17 | classic_bic | 0.87 | 0.05 | 0.00 | 1.02 | 0.10 |
| ar1 | 0.5/0.3/0.17 | panel_hac | 0.97 | 0.58 | 0.05 | 1.92 | 0.28 |
| ar1 | 0.5/0.3/0.17 | panel_mc | 0.93 | 0.62 | 0.22 | 3.47 | 0.80 |
| ar1 | 1/0.6/0.35 | classic_bic | 1.00 | 1.00 | 0.00 | 2.00 | 0.00 |
| ar1 | 1/0.6/0.35 | panel_hac | 1.00 | 0.98 | 0.58 | 2.77 | 0.20 |
| ar1 | 1/0.6/0.35 | panel_mc | 1.00 | 0.95 | 0.53 | 3.62 | 0.73 |
| ar1 | 0.5/0.15/0.07 | classic_bic | 0.95 | 0.00 | 0.00 | 1.00 | 0.05 |
| ar1 | 0.5/0.15/0.07 | panel_hac | 0.95 | 0.02 | 0.02 | 1.12 | 0.13 |
| ar1 | 0.5/0.15/0.07 | panel_mc | 0.97 | 0.25 | 0.15 | 3.62 | 0.92 |
| ar1 | 1/0.3/0.15 | classic_bic | 1.00 | 0.03 | 0.00 | 1.03 | 0.00 |
| ar1 | 1/0.3/0.15 | panel_hac | 1.00 | 0.58 | 0.00 | 1.90 | 0.32 |
| ar1 | 1/0.3/0.15 | panel_mc | 1.00 | 0.57 | 0.08 | 3.68 | 0.93 |
| iid | 0.5/0.3/0.17 | classic_bic | 0.97 | 0.03 | 0.00 | 1.02 | 0.02 |
| iid | 0.5/0.3/0.17 | panel_hac | 1.00 | 0.88 | 0.17 | 2.32 | 0.25 |
| iid | 0.5/0.3/0.17 | panel_mc | 1.00 | 0.90 | 0.27 | 2.48 | 0.27 |
| iid | 1/0.6/0.35 | classic_bic | 1.00 | 1.00 | 0.00 | 2.00 | 0.00 |
| iid | 1/0.6/0.35 | panel_hac | 1.00 | 1.00 | 0.95 | 3.02 | 0.07 |
| iid | 1/0.6/0.35 | panel_mc | 1.00 | 1.00 | 0.95 | 3.03 | 0.08 |
| iid | 0.5/0.15/0.07 | classic_bic | 1.00 | 0.00 | 0.00 | 1.00 | 0.00 |
| iid | 0.5/0.15/0.07 | panel_hac | 1.00 | 0.20 | 0.00 | 1.45 | 0.23 |
| iid | 0.5/0.15/0.07 | panel_mc | 1.00 | 0.25 | 0.00 | 1.60 | 0.33 |
| iid | 1/0.3/0.15 | classic_bic | 1.00 | 0.00 | 0.00 | 1.00 | 0.00 |
| iid | 1/0.3/0.15 | panel_hac | 1.00 | 0.90 | 0.05 | 2.23 | 0.28 |
| iid | 1/0.3/0.15 | panel_mc | 1.00 | 0.92 | 0.05 | 2.25 | 0.27 |

Reading: pooled methods find medium and small breaks much better than classic BIC. For sizes 1/0.6/0.35 with iid errors the small break is found in 95% of runs by `panel_mc` and `panel_hac` and in 0% by classic BIC. When the small jump is too weak (pooled signal around 135) no method finds it. `panel_hac` loses nothing to `panel_mc` under iid errors despite the large jump (medium break found in 90% against 92% of runs for sizes 1/0.3/0.15). With AR(1) errors (phi = 0.3) the iid-calibrated `panel_mc` has false positives in most runs. Results for `panel_hac` were identical before and after the bandwidth fix.

## 8. Strong factor, T = 300, factor variance 2, same-direction jumps (staircase), 60 replications (after the bandwidth fix)

| design | method | breaks found (mean) | exact | spurious | found | exact before fix |
|---|---|---|---|---|---|---|
| N=20, 1 factor, jumps in factor space | classic_bic | 1.20 | 0.20 | 0.65 | 0.27 | 0.20 |
| N=20, 1 factor, jumps in factor space | panel_factor | 1.38 | 0.30 | 0.68 | 0.30 | 0.30 |
| N=20, 1 factor, jumps in factor space | panel_factor_nowhite | 1.63 | 0.47 | 0.77 | 0.33 | 0.47 |
| N=20, 1 factor, jumps in factor space | panel_hac | 1.43 | 0.37 | 0.72 | 0.30 | 0.37 |
| N=20, 1 factor, jumps in factor space | panel_mc | 3.50 | 0.22 | 0.93 | 0.50 | 0.22 |
| N=20, 1 factor, jumps orthogonal to factor | classic_bic | 1.45 | 0.42 | 0.18 | 0.62 | 0.42 |
| N=20, 1 factor, jumps orthogonal to factor | panel_factor | 2.05 | 0.95 | 0.05 | 1.00 | 0.20 |
| N=20, 1 factor, jumps orthogonal to factor | panel_factor_nowhite | 1.82 | 0.68 | 0.23 | 0.77 | 0.02 |
| N=20, 1 factor, jumps orthogonal to factor | panel_hac | 1.87 | 0.70 | 0.23 | 0.80 | 0.67 |
| N=20, 1 factor, jumps orthogonal to factor | panel_mc | 3.82 | 0.13 | 0.85 | 0.88 | 0.13 |
| N=50, 1 factor, jumps in factor space | classic_bic | 0.93 | 0.00 | 0.62 | 0.16 | 0.00 |
| N=50, 1 factor, jumps in factor space | panel_factor | 1.02 | 0.05 | 0.68 | 0.16 | 0.05 |
| N=50, 1 factor, jumps in factor space | panel_factor_nowhite | 1.08 | 0.08 | 0.68 | 0.18 | 0.08 |
| N=50, 1 factor, jumps in factor space | panel_hac | 1.03 | 0.03 | 0.67 | 0.18 | 0.03 |
| N=50, 1 factor, jumps in factor space | panel_mc | 4.22 | 0.05 | 0.98 | 0.32 | 0.05 |
| N=50, 1 factor, jumps orthogonal | classic_bic | 0.83 | 0.00 | 0.37 | 0.23 | 0.00 |
| N=50, 1 factor, jumps orthogonal | panel_factor | 2.03 | 0.97 | 0.05 | 0.98 | 0.22 |
| N=50, 1 factor, jumps orthogonal | panel_factor_nowhite | 1.15 | 0.12 | 0.48 | 0.30 | 0.00 |
| N=50, 1 factor, jumps orthogonal | panel_hac | 1.13 | 0.13 | 0.50 | 0.30 | 0.13 |
| N=50, 1 factor, jumps orthogonal | panel_mc | 4.75 | 0.07 | 0.95 | 0.47 | 0.07 |
| N=20, 2 factors, jumps in factor space | classic_bic | 1.02 | 0.02 | 0.42 | 0.30 | 0.02 |
| N=20, 2 factors, jumps in factor space | panel_factor | 1.32 | 0.32 | 0.65 | 0.28 | 0.32 |
| N=20, 2 factors, jumps in factor space | panel_factor_nowhite | 1.38 | 0.38 | 0.60 | 0.33 | 0.38 |
| N=20, 2 factors, jumps in factor space | panel_hac | 1.30 | 0.30 | 0.53 | 0.33 | 0.30 |
| N=20, 2 factors, jumps in factor space | panel_mc | 3.25 | 0.22 | 0.97 | 0.42 | 0.22 |
| N=20, 2 factors, jumps orthogonal | classic_bic | 1.02 | 0.02 | 0.13 | 0.44 | 0.02 |
| N=20, 2 factors, jumps orthogonal | panel_factor | 2.03 | 0.97 | 0.03 | 1.00 | 0.28 |
| N=20, 2 factors, jumps orthogonal | panel_factor_nowhite | 1.53 | 0.43 | 0.27 | 0.62 | 0.00 |
| N=20, 2 factors, jumps orthogonal | panel_hac | 1.48 | 0.38 | 0.27 | 0.60 | 0.38 |
| N=20, 2 factors, jumps orthogonal | panel_mc | 4.20 | 0.08 | 0.93 | 0.77 | 0.08 |
| N=20, 1 factor, null | classic_bic | 0.02 | 0.98 | 0.02 | – | 0.98 |
| N=20, 1 factor, null | panel_factor | 0.10 | 0.92 | 0.08 | – | 0.92 |
| N=20, 1 factor, null | panel_factor_nowhite | 0.15 | 0.87 | 0.13 | – | 0.88 |
| N=20, 1 factor, null | panel_hac | 0.03 | 0.97 | 0.03 | – | 0.97 |
| N=20, 1 factor, null | panel_mc | 1.93 | 0.17 | 0.83 | – | 0.17 |
| N=50, 1 factor, null | classic_bic | 0.00 | 1.00 | 0.00 | – | 1.00 |
| N=50, 1 factor, null | panel_factor | 0.00 | 1.00 | 0.00 | – | 1.00 |
| N=50, 1 factor, null | panel_factor_nowhite | 0.03 | 0.97 | 0.03 | – | 0.97 |
| N=50, 1 factor, null | panel_hac | 0.00 | 1.00 | 0.00 | – | 1.00 |
| N=50, 1 factor, null | panel_mc | 3.13 | 0.05 | 0.95 | – | 0.05 |

Reading: for jumps orthogonal to the factor space, factor projection finds the exact count in 95%, 97% and 97% of runs in the three designs, against 70%, 13% and 38% for `panel_hac` on the raw panel. Before the first-pass bandwidth fix it was 20%, 22% and 28%. Whitening of the factor series matters: without it the exact-count rate is 68%, 12% and 43%, with 23 to 48% false positives. For jumps inside the factor space every method is weak (exact count at most 47% at N = 20 and at most 8% at N = 50): the effective problem is one noisy factor series with a borderline signal. On the null, `panel_factor` has 8% false positives at N = 20 and 0% at N = 50; `panel_hac` has 3% and 0%; the iid-calibrated `panel_mc` fails (83% and 95%).

## Not covered

Heavy-tailed or non-Gaussian errors, time-varying covariance, sparse breaks, strong dependence with larger T, larger factor structures (more than two factors, weak factors), comparison with `tree` and `strucchange`, and real data.

## 9. Trend breaks (`type = "trend"`), T = 200, N = 10, 60 replications

Scenarios add a common baseline slope of 0.01 per series (random sign) where stated; slope breaks are continuous kinks with per-series slope change 0.02, 0.04 or 0.03 at date 100; the level jump is 0.6 per series. `mean_mc` is `art_fit` with the Monte Carlo penalty; `trend_mc` is the same with `type = "trend"`. *found* uses a 5-date tolerance.

| scenario | method | breaks found (mean) | exact | any break | found | mean abs error (dates) |
|---|---|---|---|---|---|---|
| null, no trend | classic_bic | 0.00 | 1.00 | 0.00 | – | – |
| null, no trend | mean_mc | 0.02 | 0.98 | 0.02 | – | – |
| null, no trend | trend_mc | 0.07 | 0.93 | 0.07 | – | – |
| null, common trend 0.01 | classic_bic | 1.53 | 0.00 | 1.00 | – | – |
| null, common trend 0.01 | mean_mc | 3.25 | 0.00 | 1.00 | – | – |
| null, common trend 0.01 | trend_mc | 0.07 | 0.93 | 0.07 | – | – |
| slope break 0.02 | classic_bic | 1.45 | 0.55 | 1.00 | 0.02 | 31.9 |
| slope break 0.02 | mean_mc | 3.03 | 0.00 | 1.00 | 0.10 | 15.2 |
| slope break 0.02 | trend_mc | 1.05 | 0.95 | 1.00 | 0.25 | 11.7 |
| slope break 0.04 | classic_bic | 3.00 | 0.00 | 1.00 | 0.15 | 10.8 |
| slope break 0.04 | mean_mc | 4.45 | 0.00 | 1.00 | 0.25 | 9.1 |
| slope break 0.04 | trend_mc | 1.05 | 0.95 | 1.00 | 0.43 | 7.0 |
| slope break 0.03 + trend | classic_bic | 3.27 | 0.00 | 1.00 | 0.23 | 14.8 |
| slope break 0.03 + trend | mean_mc | 5.42 | 0.00 | 1.00 | 0.23 | 10.2 |
| slope break 0.03 + trend | trend_mc | 1.08 | 0.92 | 1.00 | 0.40 | 8.9 |
| level jump 0.6 (no slope) | classic_bic | 1.00 | 1.00 | 1.00 | 0.98 | 0.8 |
| level jump 0.6 (no slope) | mean_mc | 1.03 | 0.97 | 1.00 | 0.98 | 0.8 |
| level jump 0.6 (no slope) | trend_mc | 1.12 | 0.83 | 0.97 | 0.93 | 1.1 |

Reading: a common trend with no break makes every mean-based method find breaks (100% of runs), while `trend_mc` has 7% false positives (2% for `mean_mc` when there is no trend at all). For slope breaks `trend_mc` finds exactly one break in 92 to 95% of runs, while `mean_mc` returns 3 to 5 steps. Localizing a kink is hard: only 25 to 43% of `trend_mc` dates fall within 5 dates of the truth. For a pure level jump `trend_mc` detects the break in 97% of runs but with a lower exact-count rate (83% against 97%) and a larger error. The trend criterion is not covered by the theory in the notes; these are simulation results only.

## 10. Trend breaks with classic and join-point methods, T = 200, N = 10, 60 replications

Methods: `classic_bic` = classic ART (mean); `classic_trend` = classic ART with line fits; `mean_mc`, `trend_mc` = `art_fit` with Monte Carlo penalty (type mean, trend); `kink_mc` = `art_kink`. Slope changes are continuous kinks (per-series RMS slope change as stated), with a common baseline slope of 0.01 per series where stated. Two kinks are at dates 67 and 133.

| scenario | method | breaks found (mean) | exact | any break | found | mean abs error (dates) |
|---|---|---|---|---|---|---|
| null, no trend | classic_bic | 0.00 | 1.00 | 0.00 | – | – |
| null, no trend | classic_trend | 0.00 | 1.00 | 0.00 | – | – |
| null, no trend | mean_mc | 0.02 | 0.98 | 0.02 | – | – |
| null, no trend | trend_mc | 0.02 | 0.98 | 0.02 | – | – |
| null, no trend | kink_mc | 0.05 | 0.95 | 0.05 | – | – |
| null, common trend 0.01 | classic_bic | 1.62 | 0.00 | 1.00 | – | – |
| null, common trend 0.01 | classic_trend | 0.00 | 1.00 | 0.00 | – | – |
| null, common trend 0.01 | mean_mc | 3.45 | 0.00 | 1.00 | – | – |
| null, common trend 0.01 | trend_mc | 0.05 | 0.95 | 0.05 | – | – |
| null, common trend 0.01 | kink_mc | 0.05 | 0.97 | 0.03 | – | – |
| one kink 0.02 | classic_bic | 1.53 | 0.47 | 1.00 | 0.00 | 32.5 |
| one kink 0.02 | classic_trend | 0.90 | 0.90 | 0.90 | 0.25 | 10.9 |
| one kink 0.02 | mean_mc | 2.82 | 0.00 | 1.00 | 0.08 | 15.3 |
| one kink 0.02 | trend_mc | 1.00 | 1.00 | 1.00 | 0.28 | 11.5 |
| one kink 0.02 | kink_mc | 1.03 | 0.97 | 1.00 | 0.78 | 3.8 |
| one kink 0.04 | classic_bic | 3.03 | 0.00 | 1.00 | 0.07 | 12.0 |
| one kink 0.04 | classic_trend | 1.00 | 1.00 | 1.00 | 0.52 | 5.7 |
| one kink 0.04 | mean_mc | 4.47 | 0.00 | 1.00 | 0.27 | 9.6 |
| one kink 0.04 | trend_mc | 1.08 | 0.92 | 1.00 | 0.52 | 5.8 |
| one kink 0.04 | kink_mc | 1.12 | 0.88 | 1.00 | 1.00 | 1.9 |
| two kinks 0.03 + trend | classic_bic | 3.67 | 0.00 | 1.00 | 0.25 | 12.1 |
| two kinks 0.03 + trend | classic_trend | 1.00 | 0.00 | 1.00 | 0.00 | 33.0 |
| two kinks 0.03 + trend | mean_mc | 6.83 | 0.00 | 1.00 | 0.42 | 8.1 |
| two kinks 0.03 + trend | trend_mc | 2.02 | 0.95 | 1.00 | 0.35 | 9.6 |
| two kinks 0.03 + trend | kink_mc | 2.12 | 0.88 | 1.00 | 0.73 | 4.0 |
| level jump 0.6 (no slope) | classic_bic | 1.00 | 1.00 | 1.00 | 1.00 | 0.6 |
| level jump 0.6 (no slope) | classic_trend | 0.00 | 0.00 | 0.00 | 0.00 | – |
| level jump 0.6 (no slope) | mean_mc | 1.07 | 0.93 | 1.00 | 1.00 | 0.6 |
| level jump 0.6 (no slope) | trend_mc | 0.97 | 0.93 | 0.95 | 0.87 | 1.8 |
| level jump 0.6 (no slope) | kink_mc | 0.50 | 0.13 | 0.30 | 0.02 | 39.9 |

Reading: with a common trend and no break every mean-based method (classic ART, `mean_mc`) finds breaks in every run; the trend and kink methods have 3 to 5% false positives, and classic ART with line fits has none. For kinks, `kink_mc` locates the dates best (73 to 100% within 5 dates, mean error 2 to 4 dates), `trend_mc` finds the right count but imprecise dates (28 to 52% within 5 dates), and the mean-based methods return 3 to 7 steps. For a level jump `mean_mc` is best, `trend_mc` is close, `kink_mc` mostly misses it (30% detection), and classic ART with line fits finds nothing. Classic ART with line fits finds a single break in the two-kink design, about 33 dates from the nearest kink on average. The trend and kink criteria are not covered by the theory in the main notes beyond the exact null law and signal scale in section 10 of the notes document.

## 11. Date error against T for a kink, N = 10, one kink in the middle, slope change g / T with g = 4, 30 replications

| T | found (separate lines) | found (join points) | mean abs error separate lines | mean abs error join points |
|---|---|---|---|---|
| 100 | 0.87 | 0.87 | 6.9 | 3.0 |
| 200 | 0.97 | 0.93 | 9.7 | 3.7 |
| 400 | 0.93 | 0.93 | 19.2 | 5.0 |
| 800 | 1.00 | 0.93 | 24.0 | 7.4 |

Reading: as a fraction of T the mean error falls like T^-0.36 for separate lines and T^-0.57 for join points (log-log slopes of 0.64 and 0.43 for the error in dates). The heuristic predictions are -1/3 and -1/2. Four values of T and 30 replications make this rough agreement only.

## 12. Kink intervals: predicted against simulated standard deviation and coverage

Gaussian errors, sigma = 1, one knot in the middle, conditional on exactly one knot found (about 93% of runs). `art_ci()` on `art_kink` fits; 1000 replications. Predicted sd = 2 / (norm(gamma) sqrt(T tau (1 - tau))) dates.

| T | N | slope change per series | predicted sd | empirical sd | mean estimated sd | coverage 95% | mean width |
|---|---|---|---|---|---|---|---|
| 200 | 10 | 0.04 | 2.24 | 2.38 | 2.26 | 0.931 | 8.8 |
| 200 | 10 | 0.02 | 4.47 | 5.12 | 4.62 | 0.915 | 18.1 |
| 200 | 50 | 0.015 | 2.67 | 3.04 | 2.81 | 0.926 | 11.0 |
| 400 | 10 | 0.015 | 4.22 | 4.61 | 4.25 | 0.928 | 16.7 |

## 13. Kink intervals when the per-series signal is small, T = 200, pooled norm(gamma)^2 T^3 = 32000, 300 replications

| N | per-series gamma^2 T^3 | found one knot | plain predicted sd | predicted inflation | robust sd | empirical sd | mean estimated sd | coverage 95% |
|---|---|---|---|---|---|---|---|---|
| 10 | 3200 | 0.92 | 4.47 | 1.03 | 4.45 | 5.56 | 4.63 | 0.888 |
| 50 | 640 | 0.92 | 4.47 | 1.14 | 5.19 | 6.06 | 5.14 | 0.917 |
| 200 | 160 | 0.92 | 4.47 | 1.48 | 7.41 | 8.46 | 6.61 | 0.888 |
| 500 | 64 | 0.91 | 4.47 | 2.00 | 8.90 | 13.29 | 8.29 | 0.850 |

## 14. Two kinks, T = 300, N = 10, slope change 0.03 at dates 100 and 200, 300 replications

Exactly two knots found in 0.91 of runs.

| knot | empirical sd | robust sd | mean estimated sd | coverage 95% |
|---|---|---|---|---|
| 1 | 3.02 | 2.97 | 2.94 | 0.945 |
| 2 | 3.40 | 2.97 | 2.93 | 0.908 |

Reading (sections 12 to 14): the predicted standard deviation matches the central part of the error distribution, and the correction factor matters when the per-series signal is small (at N = 500 the robust sd is 8.9 against a corrected prediction of 8.9, and 4.47 uncorrected). Normal intervals undercover by 1 to 10 points because a few runs mislocate the knot badly; the empirical sd exceeds the robust sd. The theory is a sketch for one knot and heuristic for several; see section 11 of the notes document.

## 15. Kink search: effective signal of an uncovered knot (geometry behind the detectability lemmas)

T = 2000, true knot in the middle, N = 1, slope change 1. ‖h‖² is the squared distance of the hinge at the true knot from the span of 1, t and the hinges at the existing knots. With no other knot, ‖h‖²/T³ = 0.005208 = 1/192 exactly. With one neighbour at distance D the prediction is D² T c₂ with c₂ = 1/16; with a knot on each side at distance D the lower bound of Lemma 25 is D³/24.

| D | one neighbour: ‖h‖² / (D² T c₂) | knots on both sides: ‖h‖² / D³ | both sides: ‖h‖² / (D³/24) |
|---|---|---|---|
| 10 | 0.987 | 0.163 | 3.90 |
| 20 | 0.974 | 0.157 | 3.78 |
| 40 | 0.948 | 0.149 | 3.57 |
| 80 | 0.900 | 0.134 | 3.23 |
| 160 | 0.812 | 0.113 | 2.70 |

Reading: one neighbour leaves a quadratic effective signal, neighbours on both sides a cubic one, and the window bound of Lemma 25 holds in every row (by a factor of 2.7 to 3.9).

## 16. Kink search: repartition from perturbed starts, T = 300, N = 10, slope change 0.04 at 100 and 200, 30 replications, no insertions (cap 2)

| pattern | start offset (dates) | two knots at the end | mean abs error | sd of error |
|---|---|---|---|---|
| toward | 0 | 1.00 | 1.67 | 2.35 |
| toward | 10 | 1.00 | 1.67 | 2.35 |
| toward | 20 | 1.00 | 1.67 | 2.35 |
| toward | 40 | 1.00 | 1.67 | 2.35 |
| same | 0 | 1.00 | 1.67 | 2.35 |
| same | 10 | 1.00 | 1.67 | 2.35 |
| same | 20 | 1.00 | 1.67 | 2.35 |
| same | 40 | 1.00 | 1.67 | 2.35 |

Reading: all eight configurations end at the same knots, so repartition does not depend on the start for offsets up to 40 dates. The predicted sd is about 1.9 dates.

## 17. Same-direction against offsetting kinks, from scratch, T = 300, N = 10, slope changes of size s at dates 100 and 200, 30 replications

| pattern | s | mean knots found | exactly two | found 100 | found 200 | max knots |
|---|---|---|---|---|---|---|
| same direction (staircase of slopes) | 0.02 | 2.10 | 0.90 | 0.67 | 0.73 | 3 |
| same direction (staircase of slopes) | 0.03 | 2.13 | 0.87 | 0.87 | 0.93 | 3 |
| same direction (staircase of slopes) | 0.04 | 2.27 | 0.73 | 0.83 | 1.00 | 3 |
| opposite (tent) | 0.02 | 2.13 | 0.87 | 0.70 | 0.63 | 3 |
| opposite (tent) | 0.03 | 2.10 | 0.90 | 0.90 | 0.93 | 3 |
| opposite (tent) | 0.04 | 2.27 | 0.73 | 0.90 | 0.97 | 3 |

Reading: the tent is no harder than the same-direction case at this spacing. The exact-count rate falls from 0.90 to 0.73 as the signal grows because of occasional duplicate knots next to true ones.

## 18. Thinning for kinks: two-kink designs, T = 300, N = 10, slope changes of size s at dates 100 and 200, 400 replications

`art_kink` with the Monte Carlo penalty (200 calibration simulations). A duplicate run has two knots within 30 dates of the same true knot; a far-knot run has a knot more than 30 dates from both true knots. *Before* is the result of the search, *after* the result of the thinning step.

| pattern | s | exact count before | exact count after | duplicate runs before | duplicate runs after | far-knot runs before | far-knot runs after |
|---|---|---|---|---|---|---|---|
| same direction | 0.03 | 0.892 | 0.912 | 0.043 | 0.030 | 0.072 | 0.065 |
| same direction | 0.04 | 0.875 | 0.905 | 0.065 | 0.037 | 0.070 | 0.065 |
| tent | 0.03 | 0.870 | 0.895 | 0.060 | 0.037 | 0.075 | 0.072 |
| tent | 0.04 | 0.835 | 0.892 | 0.087 | 0.040 | 0.085 | 0.075 |

## 19. Thinning for kinks: one-kink designs, 400 replications

| T | N | slope change | exact count before | exact count after | mean abs error before | mean abs error after |
|---|---|---|---|---|---|---|
| 200 | 10 | 0.04 | 0.912 | 0.912 | 2.07 | 2.07 |
| 200 | 10 | 0.02 | 0.920 | 0.920 | 4.56 | 4.56 |
| 200 | 50 | 0.015 | 0.922 | 0.922 | 2.29 | 2.29 |

## 20. Null false-positive rate of `art_kink` against the number of calibration simulations, T = 200, N = 10, 400 replications

| calibration simulations | penalty (sigma = 1) | P(any knot on a null panel) |
|---|---|---|
| 200 | 23.38 | 0.095 |
| 500 | 23.87 | 0.073 |
| 2000 | 24.59 | 0.053 |

Reading (sections 18 to 20): thinning removes about half of the duplicate knots and raises the exact-count rate by 2 to 6 points in two-kink designs; one-kink fits are unchanged because thinning only acts with two or more knots. The null rate depends on the calibration size: with 2000 simulations it is 5.3% at the 5% level, so the default is now 2000. The far-knot rates in section 18 come from runs with 200 simulations and are above the nominal level for that reason.

## 21. Kink search theory checks (section 14 of the notes)

Removal-cost identity: over 83 knots in simulated two-kink fits (T = 300, N = 10, slope change 0.04), the removal cost of a knot equals the squared norm of its estimated slope change times the squared norm of its hinge residualized on the other regressors, with maximum relative error 1.4e-14. The local bound L^2 R^2 / (3 (L + R)) on the squared norm of the residualized hinge divided by the actual value has minimum 2.01 and median 3.2, so it is valid and conservative by a factor of about three.

Covering condition: for this design (lambda = 1/3, three knots allowed, C0' = 64) the sufficient condition of Theorem 38 requires norm(gamma)^2 T^3 of at least 1.3e9; the design has 4.3e5, a factor of about 3000 short. Observed covering distances are a few dates.

## 22. Union-bound penalty against Monte Carlo penalty for kinks, T = 300, N = 10, 200 replications

Two kinks at dates 100 and 200 with the same sign pattern, slope change s per series (s = 0 is a null panel). A far-knot run has a knot more than 30 dates from both true knots. The union-bound penalty is N + 4 sqrt(N log T) + 10 log T = 102.3; the Monte Carlo penalty is 26.5.

| slope change | penalty | alpha | exact count | far-knot runs | both knots found |
|---|---|---|---|---|---|
| 0 | Monte Carlo | 26.5 | 0.940 | 0.060 | – |
| 0 | union bound | 102.3 | 1.000 | 0.000 | – |
| 0.03 | Monte Carlo | 26.5 | 0.925 | 0.050 | 0.80 |
| 0.03 | union bound | 102.3 | 0.460 | 0.540 | 0.35 |
| 0.04 | Monte Carlo | 26.5 | 0.920 | 0.055 | 0.94 |
| 0.04 | union bound | 102.3 | 1.000 | 0.000 | 0.96 |

Reading: the large penalty has no false positives and no far knots at the strong signal, but at s = 0.03 it ends with a single knot between the two true ones in 54% of runs. The two kinks have the same sign, so the mean path is close to one bend; a second knot lowers the sum of squares by less than the penalty. This is a local minimum of the penalized criterion, not a search failure. The Monte Carlo penalty adds the second knot, at the cost of about 5 to 6% false positives on null panels.

## 23. Do extra knots survive thinning? Monte Carlo against union-bound penalty, N = 10, 100 replications per cell

Knots share the same sign pattern. An extra run ends with more knots than true knots; a far knot is more than 30 dates from every true knot; a duplicate is a second knot within 30 dates of a true knot; a missed window has no knot within 30 dates of some true knot.

| design | penalty | alpha | exact count | extra runs | far knot | duplicate | missed window |
|---|---|---|---|---|---|---|---|
| 2 kinks, T=300, s=0.04 | Monte Carlo | 25.6 | 0.95 | 0.05 | 0.04 | 0.02 | 0.00 |
| 2 kinks, T=300, s=0.04 | union bound | 98.8 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 2 kinks, T=300, s=0.06 | Monte Carlo | 25.6 | 0.94 | 0.06 | 0.04 | 0.03 | 0.00 |
| 2 kinks, T=300, s=0.06 | union bound | 98.8 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 3 kinks, T=400, s=0.05 | Monte Carlo | 24.6 | 0.93 | 0.07 | 0.06 | 0.01 | 0.00 |
| 3 kinks, T=400, s=0.05 | union bound | 101.4 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 |

## 24. Union-bound penalty, 300 replications per design: extras and misfit

| design | extra runs | missed windows | median misfit | max misfit | bound rho0 (Lemma 41) |
|---|---|---|---|---|---|
| 2 kinks, T=300, s=0.04 | 0.00 | 0.00 | 1.28 | 12.4 | 18820 |
| 2 kinks, T=300, s=0.06 | 0.00 | 0.00 | 1.87 | 11.9 | 28230 |
| 3 kinks, T=400, s=0.05 | 0.00 | 0.00 | 2.80 | 23.8 | 55670 |

Reading (sections 23 and 24): with the union-bound penalty no extra knot survived in 900 runs and no true knot was missed; the misfit at the output is of order sigma^2 per knot. The Monte Carlo penalty leaves extras at about its nominal false-positive level, as expected if an extra knot without signal behaves like a null test (Lemma 44 of the notes). The provable misfit bound is valid but more than a thousand times larger than observed.

## 25. Thinning at the penalty against thinning at the 2N noise level (alpha_thin), N = 10, 120 replications

| design | thinning | exact count | extra runs | missed window |
|---|---|---|---|---|
| 2 kinks, T=300, s=0.025 | alpha | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.025 | 2N level | 1.000 | 0.000 | 0.00 |
| 2 kinks, T=300, s=0.03 | alpha | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.03 | 2N level | 1.000 | 0.000 | 0.00 |
| 2 kinks, T=300, s=0.04 | alpha | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.04 | 2N level | 1.000 | 0.000 | 0.00 |
| 3 kinks, T=400, s=0.05 | alpha | 0.892 | 0.108 | 0.00 |
| 3 kinks, T=400, s=0.05 | 2N level | 1.000 | 0.000 | 0.00 |
| null, T=300 | alpha | 0.967 | 0.033 | 0.00 |
| null, T=300 | 2N level | 0.967 | 0.033 | 0.00 |

## 26. Cost of thinning at the 2N level: a weak true knot beside a strong one (T = 300, N = 10, strong slope 0.04 at 100, weak knot at 200, 150 replications)

| weak slope | weak population gain | thinning | weak knot found | strong knot found | exact two knots | extra |
|---|---|---|---|---|---|---|
| 0.004 | 8.4 | alpha | 0.06 | 1.00 | 0.11 | 0.007 |
| 0.004 | 8.4 | 2N level | 0.00 | 1.00 | 0.00 | 0.000 |
| 0.005 | 13.1 | alpha | 0.13 | 1.00 | 0.20 | 0.020 |
| 0.005 | 13.1 | 2N level | 0.00 | 1.00 | 0.00 | 0.000 |
| 0.006 | 18.9 | alpha | 0.29 | 1.00 | 0.34 | 0.013 |
| 0.006 | 18.9 | 2N level | 0.01 | 1.00 | 0.01 | 0.000 |
| 0.008 | 33.6 | alpha | 0.65 | 1.00 | 0.72 | 0.040 |
| 0.008 | 33.6 | 2N level | 0.09 | 1.00 | 0.09 | 0.000 |
| 0.012 | 75.6 | alpha | 0.98 | 1.00 | 0.97 | 0.027 |
| 0.012 | 75.6 | 2N level | 0.83 | 1.00 | 0.83 | 0.000 |

Thinning at the 2N level removes every duplicate in the alternative designs and leaves the null rate unchanged, but it deletes weak true knots (found 0.09 against 0.65 at weak gain 33.6). A lone knot is never thinned, so single-kink power checks do not exercise this trade-off.


## 27. Surviving duplicates with thinning off (200 runs, two kinks, slope 0.04, T = 300, N = 10): 19 runs ended with a duplicate pair. Re-optimized removal cost 13.7 to 39.0 (median 22.8); below alpha in 11 of 19, below the 2N level 52.8 in all 19. See study_kink_traps.csv.

## 28. Surviving extras by kind, 200 two-kink panels (T = 300, N = 10, slope 0.04, Monte Carlo penalty with 500 draws)

Search output without thinning: 171 exact, 10 isolated far knots, 19 duplicate pairs. After thinning at the penalty: 182 exact, 10 far, 8 duplicates. Thinning at the 2N level and a fixed penalty of 52.8 sigma^2: 200 exact each. So the extras at the Monte Carlo level are about half isolated knots and half duplicates, and the 2N level works as a stricter overall level rather than a duplicate-specific repair.

## 29. Pair test against thinning at the penalty

On the 19 duplicate runs the pair test (level 0.01) merged 11 and kept 8, the same counts as thinning at the penalty; bootstrap quantiles 24.3 to 28.6 against a penalty of about 25. The 8 kept pairs are rejections of the single-knot null in 8 of 200 runs at nominal level 0.01 (see study_kink_pair_diag.csv). In the four designs of study_kink_pair1_a.csv its exact counts equal those of thinning at the penalty.

### Pair test and thinning in four designs (120 replications; 3-kink design not run)

| design | method | exact count | extra runs | missed window |
|---|---|---|---|---|
| 2 kinks, T=300, s=0.025 | alpha | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.025 | twoN | 1.000 | 0.000 | 0.00 |
| 2 kinks, T=300, s=0.025 | pair05 | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.025 | pair01 | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.03 | alpha | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.03 | twoN | 1.000 | 0.000 | 0.00 |
| 2 kinks, T=300, s=0.03 | pair05 | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.03 | pair01 | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.04 | alpha | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.04 | twoN | 1.000 | 0.000 | 0.00 |
| 2 kinks, T=300, s=0.04 | pair05 | 0.933 | 0.067 | 0.00 |
| 2 kinks, T=300, s=0.04 | pair01 | 0.933 | 0.067 | 0.00 |
| null, T=300 | alpha | 0.967 | 0.033 | 0.00 |
| null, T=300 | twoN | 0.967 | 0.033 | 0.00 |
| null, T=300 | pair05 | 0.967 | 0.033 | 0.00 |
| null, T=300 | pair01 | 0.967 | 0.033 | 0.00 |

### Weak knot beside a strong one (100 replications): thinning at the penalty, at the 2N level, and the pair test

| weak slope | weak gain | method | weak knot found | exact two | extra |
|---|---|---|---|---|---|
| 0.006 | 18.9 | alpha | 0.28 | 0.33 | 0.01 |
| 0.006 | 18.9 | twoN | 0.02 | 0.02 | 0.00 |
| 0.006 | 18.9 | pair01 | 0.28 | 0.33 | 0.01 |
| 0.008 | 33.6 | alpha | 0.66 | 0.72 | 0.04 |
| 0.008 | 33.6 | twoN | 0.09 | 0.10 | 0.00 |
| 0.008 | 33.6 | pair01 | 0.66 | 0.72 | 0.04 |
| 0.012 | 75.6 | alpha | 0.97 | 0.96 | 0.03 |
| 0.012 | 75.6 | twoN | 0.85 | 0.86 | 0.00 |
| 0.012 | 75.6 | pair01 | 0.97 | 0.97 | 0.02 |

### Weak knot found (100 replications): 2N thinning against a fixed penalty of 52.8

| weak slope | 2N thinning | fixed penalty |
|---|---|---|
| 0.008 | 0.09 | 0.02 |
| 0.012 | 0.85 | 0.46 |

Correction to sections 25 and 26 above: the 2N-level improvement in exact counts is mostly the removal of isolated false-positive knots, not only duplicates (section 28).

## 30. Variance of smooth projections under AR(1) errors: exact check of Lemma 48 (unit innovation variance; entries are v'Gamma v / omega^2 - 1)

| phi | T | hinge | hinge, second knot present | step | hinge, knot 10 dates away | T x hinge | T x step |
|---|---|---|---|---|---|---|---|
| 0.5 | 100 | -0.0482 | -0.0406 | -0.1466 | -0.0596 | -4.82 | -14.7 |
| 0.5 | 200 | -0.0229 | -0.0172 | -0.0749 | -0.0265 | -4.59 | -15.0 |
| 0.5 | 400 | -0.0112 | -0.0077 | -0.0379 | -0.0124 | -4.46 | -15.1 |
| 0.5 | 800 | -0.0055 | -0.0036 | -0.0190 | -0.0060 | -4.40 | -15.2 |
| 0.5 | 1600 | -0.0027 | -0.0018 | -0.0095 | -0.0029 | -4.37 | -15.3 |
| 0.8 | 100 | -0.1807 | -0.1913 | -0.4330 | -0.3131 | -18.07 | -43.3 |
| 0.8 | 200 | -0.0832 | -0.0760 | -0.2373 | -0.1553 | -16.63 | -47.5 |
| 0.8 | 400 | -0.0391 | -0.0311 | -0.1233 | -0.0756 | -15.63 | -49.3 |
| 0.8 | 800 | -0.0188 | -0.0135 | -0.0627 | -0.0371 | -15.06 | -50.2 |
| 0.8 | 1600 | -0.0092 | -0.0062 | -0.0316 | -0.0183 | -14.76 | -50.6 |

The relative error shrinks like 1/T and is negative: the true variance sits slightly below the long-run variance.

## 31. Size and power under AR(1) errors with unit stationary variance (T = 300, N = 10; omega^2 = (1+phi)/(1-phi)). The i.i.d. comparator is capped at 10 knots, so its mean knot count is censored.

| phi | design | method | reps | exact | extra runs | missed window | mean knots | omega2 estimate | omega2 true |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 | null | mc (iid) | 100 | 0.020 | 0.980 | 0.000 | 7.80 | - | 1.86 |
| 0.3 | null | hac_mc | 100 | 0.950 | 0.050 | 0.000 | 0.05 | 1.76 | 1.86 |
| 0.3 | 2 kinks s=0.04 | mc (iid) | 40 | 0.000 | 1.000 | 0.000 | 9.65 | - | 1.86 |
| 0.3 | 2 kinks s=0.04 | hac_mc | 40 | 0.950 | 0.050 | 0.000 | 2.05 | 1.75 | 1.86 |
| 0.3 | 2 kinks s=0.06 | mc (iid) | 40 | 0.000 | 1.000 | 0.000 | 9.80 | - | 1.86 |
| 0.3 | 2 kinks s=0.06 | hac_mc | 40 | 0.950 | 0.050 | 0.000 | 2.05 | 1.75 | 1.86 |
| 0.5 | null | mc (iid) | 100 | 0.000 | 1.000 | 0.000 | 9.97 | - | 3.00 |
| 0.5 | null | hac_mc | 100 | 0.940 | 0.060 | 0.000 | 0.06 | 2.75 | 3.00 |
| 0.5 | 2 kinks s=0.04 | mc (iid) | 40 | 0.000 | 1.000 | 0.000 | 10.00 | - | 3.00 |
| 0.5 | 2 kinks s=0.04 | hac_mc | 40 | 0.925 | 0.075 | 0.000 | 2.08 | 2.71 | 3.00 |
| 0.5 | 2 kinks s=0.06 | mc (iid) | 40 | 0.000 | 1.000 | 0.000 | 9.97 | - | 3.00 |
| 0.5 | 2 kinks s=0.06 | hac_mc | 40 | 0.950 | 0.050 | 0.000 | 2.05 | 2.72 | 3.00 |
| 0.7 | null | mc (iid) | 100 | 0.000 | 1.000 | 0.000 | 9.99 | - | 5.67 |
| 0.7 | null | hac_mc | 100 | 0.940 | 0.060 | 0.000 | 0.08 | 5.01 | 5.67 |
| 0.7 | 2 kinks s=0.04 | mc (iid) | 40 | 0.000 | 1.000 | 0.000 | 10.00 | - | 5.67 |
| 0.7 | 2 kinks s=0.04 | hac_mc | 40 | 0.625 | 0.025 | 0.350 | 1.68 | 5.89 | 5.67 |
| 0.7 | 2 kinks s=0.06 | mc (iid) | 40 | 0.000 | 1.000 | 0.000 | 9.97 | - | 5.67 |
| 0.7 | 2 kinks s=0.06 | hac_mc | 40 | 0.825 | 0.000 | 0.175 | 1.82 | 6.45 | 5.67 |

## 32. One kink at 150, slope 0.08, T = 300, N = 10, 400 replications: localization and 95% intervals

| phi | runs with one knot | empirical sd | predicted sd (long-run variance) | i.i.d. formula | mean estimated sd | coverage, long-run variance | coverage, residual variance |
|---|---|---|---|---|---|---|---|
| 0.3 | 0.953 | 1.30 | 1.25 | 0.91 | 1.21 | 0.945 | 0.769 |
| 0.5 | 0.963 | 1.59 | 1.58 | 0.91 | 1.52 | 0.938 | 0.670 |
| 0.7 | 0.975 | 2.16 | 2.18 | 0.91 | 2.04 | 0.936 | 0.515 |
| 0.8 | 0.988 | 2.64 | 2.75 | 0.91 | 2.54 | 0.932 | 0.448 |

## 33. Series with different persistence: null false-positive rate and power, pooled against series-specific scaling (T = 300, N = 10)

| design | reps | pooled | series |
|---|---|---|---|
| null, phi~U(0.0,0.4) | 200 | 0.050 | 0.060 |
| null, phi~U(0.0,0.8) | 200 | 0.110 | 0.075 |
| null, phi~U(0.5,0.9) | 200 | 0.180 | 0.240 |
| 2 kinks s=0.06, phi~U(0,0.8): exact count | 50 | 0.780 | 0.880 |
| 2 kinks s=0.06, phi~U(0,0.8): missed window | 50 | 0.180 | 0.120 |

A first run with pooled scaling only (study_kink_dep3.csv, different random draws) gave 0.060, 0.105 and 0.155 for the three persistence ranges.

## 34. Null false-positive rate at high persistence (hac_mc, pooled; comparator capped at 10 knots)

| T | phi | replications done | false positive | omega2 estimate | omega2 true | mean passes |
|---|---|---|---|---|---|---|
| 300 | 0.8 | 150 | 0.207 | 7.56 | 9.0 | 2.8 |
| 300 | 0.9 | 100 | 0.830 | 6.13 | 19.0 | 3.0 |
| 1000 | 0.8 | 60 | 0.117 | 8.18 | 9.0 | 2.7 |
| 1000 | 0.9 | 27 | 0.148 | 16.61 | 19.0 | 3.8 |

The T = 1000, phi = 0.9 row stopped at 27 replications on a time limit and is imprecise.

## 35. The same-sign trap and the long-run-variance step (AR(1), T = 300, N = 10, two kinks at 100 and 200, hac_mc, comparator uncapped runs capped at 10 knots)

Diagnosis (kink_dep_diag.R, phi = 0.7, slope 0.06, 40 runs, 7 with a missed window): in every missed run pass 0 returned a single knot near 150 (the two same-sign kinks look like one bend) with a long-run variance of 23 to 25, and later passes stayed at that knot with estimates of 11 to 18 against a true 5.7; a search with the true long-run variance found both kinks in all 7 runs. The first version of the option (extra knots in every pass, not saved) gave exact counts 0.725 and 0.575 for 2 and 4 extra knots; the table below is the two-phase version in the package.

| phi | slope | design | lrv_extra | reps | both kinks found | missed window | false positive (null) | omega2 estimate | omega2 true |
|---|---|---|---|---|---|---|---|---|---|
| 0.7 | 0.06 | two kinks | 0 | 40 | 0.825 | 0.175 | - | 6.45 | 5.67 |
| 0.7 | 0.06 | two kinks | 2 | 40 | 0.800 | 0.000 | - | 4.75 | 5.67 |
| 0.7 | 0.06 | two kinks | 4 | 40 | 0.750 | 0.000 | - | 4.68 | 5.67 |
| 0.7 | 0.04 | two kinks | 0 | 40 | 0.625 | 0.350 | - | 5.89 | 5.67 |
| 0.7 | 0.04 | two kinks | 2 | 40 | 0.775 | 0.000 | - | 4.70 | 5.67 |
| 0.5 | 0.04 | two kinks | 0 | 40 | 0.925 | 0.000 | - | 2.71 | 3.00 |
| 0.5 | 0.04 | two kinks | 2 | 40 | 0.875 | 0.000 | - | 2.69 | 3.00 |
| 0.3 | 0 | null | 0 | 150 | - | - | 0.067 | 1.76 | 1.86 |
| 0.3 | 0 | null | 2 | 150 | - | - | 0.067 | 1.75 | 1.86 |
| 0.5 | 0 | null | 0 | 150 | - | - | 0.067 | 2.76 | 3.00 |
| 0.5 | 0 | null | 2 | 150 | - | - | 0.067 | 2.75 | 3.00 |
| 0.7 | 0 | null | 0 | 150 | - | - | 0.087 | 4.99 | 5.67 |
| 0.7 | 0 | null | 2 | 150 | - | - | 0.087 | 4.95 | 5.67 |

## 36. Cost of lrv_extra = 2: share of single-kink runs (slope 0.08, T = 300, N = 10, 400 replications) that return exactly one knot

| phi | lrv_extra = 0 | lrv_extra = 2 |
|---|---|---|
| 0.3 | 0.953 | 0.927 |
| 0.5 | 0.963 | 0.907 |
| 0.7 | 0.975 | 0.855 |
| 0.8 | 0.988 | 0.988 |

Localization sd and interval coverage were unchanged by the option (study_kink_dep2_extra2.csv). The option removes the missed windows of the same-sign trap but raises the share of runs with an extra knot when a real kink is present, so the default is 0.

## 37. Long-run variance from an AR fitted to first differences (T = 300, N = 10, two kinks of slope 0.06 present, 100 replications; unit innovation variance)

First-difference AR(1) is accurate for AR(1) errors and unaffected by the kinks; the pass-0 estimate used before (residuals of the no-knot fit, bandwidth floor(T^(1/3))) is 22 to 43 for every error process. A second-difference AR fit was unstable (first run, not saved: mean 137 with sd 368 at phi = 0.7 against a true 11).

| errors | true omega2 | AR(1) mean | AR(1) median | AR(1) sd | AR(2) median | second-difference AR(1) median | old pass-0 mean |
|---|---|---|---|---|---|---|---|
| white | 1.00 | 1.01 | 1.01 | 0.12 | 1.02 | 1.03 | 22.3 |
| AR1 .3 | 2.04 | 2.08 | 2.09 | 0.28 | 2.10 | 2.14 | 23.2 |
| AR1 .5 | 4.00 | 4.10 | 4.01 | 0.60 | 4.26 | 4.24 | 24.4 |
| AR1 .7 | 11.11 | 11.79 | 11.74 | 1.98 | 11.62 | 13.46 | 27.6 |
| AR1 .8 | 25.00 | 28.42 | 27.74 | 6.46 | 27.65 | 39.94 | 31.6 |
| AR1 .9 | 100.00 | 133.65 | 122.17 | 47.22 | 127.24 | 383.13 | 42.7 |
| AR2 (.6,.2) | 25.00 | 1.70 | 1.62 | 0.46 | 29.11 | 0.60 | 29.9 |
| ARMA11 (.5,.5) | 9.00 | 44.40 | 45.15 | 4.89 | 4.95 | 1156.44 | 28.1 |
| MA1 (.6) | 2.56 | 8.33 | 8.28 | 0.69 | 1.17 | 1752.25 | 23.8 |

## 38. Estimate divided by truth for first-difference AR(p), p = 1 to 4 (80 replications; 10th percentile, median, 90th percentile)

| errors | true | p1 q10 | p1 med | p1 q90 | p2 q10 | p2 med | p2 q90 | p3 q10 | p3 med | p3 q90 | p4 q10 | p4 med | p4 q90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| white | 1.00 | 0.87 | 1.02 | 1.16 | 0.78 | 1.02 | 1.33 | 0.66 | 1.02 | 1.58 | 0.45 | 0.94 | 1.51 |
| AR1 .5 | 4.00 | 0.85 | 1.01 | 1.20 | 0.76 | 1.00 | 1.22 | 0.72 | 1.02 | 1.35 | 0.69 | 1.08 | 1.62 |
| AR1 .8 | 25.00 | 0.78 | 1.10 | 1.44 | 0.76 | 1.07 | 1.61 | 0.77 | 1.12 | 1.88 | 0.75 | 1.17 | 1.89 |
| AR2 (.6,.2) | 25.00 | 0.05 | 0.07 | 0.09 | 0.72 | 1.11 | 1.90 | 0.74 | 1.26 | 2.24 | 0.69 | 1.30 | 2.21 |
| AR2 (1.2,-.4) | 25.00 | 9.75 | 11.25 | 13.70 | 0.86 | 1.02 | 1.31 | 0.82 | 1.09 | 1.44 | 0.76 | 1.13 | 1.45 |
| ARMA11 (.5,.5) | 9.00 | 4.18 | 4.93 | 5.50 | 0.40 | 0.54 | 0.65 | 1.09 | 1.53 | 1.95 | 0.68 | 0.87 | 1.20 |
| ARMA11 (.8,-.4) | 9.00 | 0.10 | 0.13 | 0.14 | 0.18 | 0.26 | 0.35 | 0.31 | 0.49 | 0.81 | 0.32 | 0.80 | 1.60 |
| MA1 (.6) | 2.56 | 2.96 | 3.23 | 3.57 | 0.35 | 0.45 | 0.54 | 1.49 | 1.90 | 2.43 | 0.50 | 0.70 | 0.98 |
| MA1 (-.5) | 0.25 | 3.18 | 3.69 | 4.13 | 1.75 | 2.36 | 3.07 | 1.41 | 2.05 | 4.42 | 0.46 | 1.07 | 3.07 |

The package estimate is the geometric mean of the p = 2, 3, 4 estimates.

## 39. Power: three long-run variance modes (T = 300, N = 10, two kinks at 100 and 200, 40 replications; hac = residuals in every pass, ar_start = knot-free first pass, ar = knot-free in every pass)

| errors | slope | mode | reps | both kinks found | missed window | mean knots | omega2 estimate | omega2 true |
|---|---|---|---|---|---|---|---|---|
| AR1 0.7 | 0.06 | hac | 40 | 0.825 | 0.175 | 1.82 | 6.45 | 5.67 |
| AR1 0.7 | 0.06 | ar_start | 40 | 0.975 | 0.000 | 2.02 | 4.89 | 5.67 |
| AR1 0.7 | 0.06 | ar | 40 | 0.950 | 0.000 | 2.08 | 6.69 | 5.67 |
| AR1 0.7 | 0.04 | hac | 40 | 0.625 | 0.350 | 1.68 | 5.89 | 5.67 |
| AR1 0.7 | 0.04 | ar_start | 40 | 0.875 | 0.100 | 1.93 | 5.15 | 5.67 |
| AR1 0.7 | 0.04 | ar | 40 | 0.800 | 0.150 | 1.93 | 6.35 | 5.67 |
| AR1 0.5 | 0.04 | hac | 40 | 0.925 | 0.000 | 2.08 | 2.71 | 3.00 |
| AR1 0.5 | 0.04 | ar_start | 40 | 0.900 | 0.000 | 2.12 | 2.70 | 3.00 |
| AR1 0.5 | 0.04 | ar | 40 | 0.900 | 0.000 | 2.12 | 3.28 | 3.00 |
| AR1 0.3 | 0.04 | hac | 40 | 0.950 | 0.000 | 2.05 | 1.75 | 1.86 |
| AR1 0.3 | 0.04 | ar_start | 40 | 0.950 | 0.000 | 2.05 | 1.75 | 1.86 |
| AR1 0.3 | 0.04 | ar | 40 | 0.875 | 0.000 | 2.12 | 1.93 | 1.86 |
| AR1 0.8 | 0.06 | hac | 40 | 0.275 | 0.725 | 1.27 | 17.09 | 9.00 |
| AR1 0.8 | 0.06 | ar_start | 40 | 0.900 | 0.100 | 1.90 | 8.75 | 9.00 |
| AR1 0.8 | 0.06 | ar | 40 | 0.875 | 0.125 | 1.88 | 12.55 | 9.00 |
| ARMA(1,1) .5,.5 | 0.06 | hac | 40 | 0.875 | 0.125 | 1.88 | 4.82 | 3.86 |
| ARMA(1,1) .5,.5 | 0.06 | ar_start | 40 | 0.975 | 0.000 | 2.02 | 3.75 | 3.86 |
| MA(1) .6 | 0.06 | hac | 40 | 1.000 | 0.000 | 2.00 | 1.87 | 1.88 |
| MA(1) .6 | 0.06 | ar_start | 40 | 0.900 | 0.000 | 2.12 | 1.86 | 1.88 |
| AR(2) 1.2,-.4 | 0.06 | hac | 40 | 0.100 | 0.900 | 1.10 | 20.29 | 5.57 |
| AR(2) 1.2,-.4 | 0.06 | ar_start | 40 | 1.000 | 0.000 | 2.00 | 5.90 | 5.57 |
| white | 0.04 | hac | 40 | 1.000 | 0.000 | 2.00 | 1.01 | 1.00 |
| white | 0.04 | ar_start | 40 | 0.950 | 0.000 | 2.05 | 1.01 | 1.00 |

## 40. Null false-positive rate of the three modes (100 replications; hac at phi = 0.8 and 0.9 from section 34)

| errors | mode | false positive | mean knots | omega2 estimate | omega2 true |
|---|---|---|---|---|---|
| AR1 0.3 | hac | 0.050 | 0.05 | 1.76 | 1.86 |
| AR1 0.3 | ar_start | 0.050 | 0.05 | 1.76 | 1.86 |
| AR1 0.3 | ar | 0.090 | 0.12 | 1.88 | 1.86 |
| AR1 0.5 | hac | 0.060 | 0.06 | 2.75 | 3.00 |
| AR1 0.5 | ar_start | 0.040 | 0.04 | 2.76 | 3.00 |
| AR1 0.5 | ar | 0.060 | 0.06 | 3.13 | 3.00 |
| AR1 0.7 | hac | 0.060 | 0.08 | 5.01 | 5.67 |
| AR1 0.7 | ar_start | 0.030 | 0.03 | 5.03 | 5.67 |
| AR1 0.7 | ar | 0.030 | 0.03 | 5.96 | 5.67 |
| AR1 0.8 | ar_start | 0.030 | 0.04 | 7.78 | 9.00 |
| AR1 0.8 | ar | 0.040 | 0.05 | 9.27 | 9.00 |
| AR1 0.9 | ar_start | 0.000 | 0.00 | 14.75 | 19.00 |
| AR1 0.9 | ar | 0.040 | 0.04 | 19.27 | 19.00 |
| AR(2) 1.2,-.4 | hac | 0.060 | 0.08 | 5.55 | 5.57 |
| AR(2) 1.2,-.4 | ar_start | 0.050 | 0.05 | 5.58 | 5.57 |
| AR(2) 1.2,-.4 | ar | 0.070 | 0.07 | 5.70 | 5.57 |
| ARMA(1,1) .5,.5 | hac | 0.070 | 0.07 | 3.70 | 3.86 |
| ARMA(1,1) .5,.5 | ar_start | 0.070 | 0.07 | 3.70 | 3.86 |
| ARMA(1,1) .5,.5 | ar | 0.140 | 0.17 | 3.41 | 3.86 |
| MA(1) .6 | hac | 0.060 | 0.06 | 1.85 | 1.88 |
| MA(1) .6 | ar_start | 0.060 | 0.06 | 1.85 | 1.88 |
| MA(1) .6 | ar | 0.160 | 0.22 | 1.55 | 1.88 |
| ARMA(1,1) .8,-.4 | hac | 0.300 | 0.66 | 4.10 | 6.23 |
| ARMA(1,1) .8,-.4 | ar_start | 0.310 | 1.04 | 4.01 | 6.23 |
| ARMA(1,1) .8,-.4 | ar | 0.580 | 1.88 | 2.85 | 6.23 |
| AR1 0.8 (T = 300, 150 reps) | hac | 0.207 | - | 7.56 | 9.00 |
| AR1 0.9 (T = 300, 100 reps) | hac | 0.830 | - | 6.13 | 19.00 |

The tables of sections 31 to 36 were produced with lrv = "hac", the default before section 19.

## 41. Null false-positive rate with cross-sectional dependence (T = 300; Omega = equicorrelation 0.5 unless stated; none = ignore the dependence, separable = eigenvalue-weighted penalty, whiten = whitened panel, full = N x N long-run covariance from residuals; all capped at 10 knots)

| setting | N | cross | reps | false positive | mean knots |
|---|---|---|---|---|---|
| uncorrelated (control) | 10 | none | 100 | 0.060 | 0.08 |
| uncorrelated (control) | 10 | separable | 100 | 0.050 | 0.07 |
| uncorrelated (control) | 10 | whiten | 100 | 0.050 | 0.07 |
| equicorr 0.5, iid in time | 10 | none | 100 | 0.340 | 0.51 |
| equicorr 0.5, iid in time | 10 | separable | 100 | 0.050 | 0.05 |
| equicorr 0.5, iid in time | 10 | whiten | 100 | 0.050 | 0.07 |
| equicorr 0.5, AR(1) 0.5 separable | 10 | none | 100 | 0.320 | 0.48 |
| equicorr 0.5, AR(1) 0.5 separable | 10 | separable | 100 | 0.020 | 0.02 |
| equicorr 0.5, AR(1) 0.5 separable | 10 | whiten | 100 | 0.060 | 0.07 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable) | 10 | none | 100 | 0.570 | 1.63 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable) | 10 | separable | 100 | 0.310 | 0.59 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable) | 10 | whiten | 100 | 0.370 | 0.75 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable) | 10 | full | 100 | 0.130 | 0.25 |
| equicorr 0.5, N=30 | 30 | none | 80 | 0.650 | 1.38 |
| equicorr 0.5, N=30 | 30 | separable | 80 | 0.037 | 0.04 |
| equicorr 0.5, N=30 | 30 | whiten | 80 | 0.050 | 0.05 |
| equicorr 0.5, AR(1) 0.5 separable | 10 | full | 100 | 0.030 | 0.03 |

## 42. Power: two kinks at 100 and 200, equicorrelation 0.5, errors iid in time, N = 10, T = 300 (aligned jumps = same sign in all series, orthogonal = alternating signs)

| jumps | slope | cross | reps | both kinks found | missed window | mean knots |
|---|---|---|---|---|---|---|
| aligned | 0.03 | none | 40 | 0.725 | 0.025 | 2.38 |
| aligned | 0.03 | separable | 40 | 0.950 | 0.025 | 2.00 |
| aligned | 0.03 | whiten | 40 | 0.625 | 0.325 | 1.80 |
| orthogonal | 0.03 | none | 40 | 0.575 | 0.000 | 2.65 |
| orthogonal | 0.03 | separable | 40 | 0.925 | 0.025 | 2.02 |
| orthogonal | 0.03 | whiten | 40 | 0.925 | 0.000 | 2.08 |
| orthogonal | 0.015 | separable | 60 | 0.183 | 0.817 | 1.18 |
| orthogonal | 0.015 | whiten | 60 | 0.933 | 0.000 | 2.07 |
| orthogonal | 0.012 | separable | 60 | 0.017 | 0.983 | 1.03 |
| orthogonal | 0.012 | whiten | 60 | 0.917 | 0.033 | 2.02 |

## 43. One kink at 150, slope 0.08, N = 10, T = 300: localization and 95% intervals. pred_sd_sandwich uses the true long-run covariance; pred_sd_scalar is the formula that ignores the cross-sectional covariance (unit per-series variance)

| errors | jumps | cross | reps | one knot | empirical sd | predicted sd (sandwich) | predicted sd (scalar) | mean estimated sd | coverage |
|---|---|---|---|---|---|---|---|---|---|
| iid in time | aligned | none | 150 | 0.693 | 2.18 | 2.14 | 0.91 | 0.91 | 0.538 |
| iid in time | aligned | separable | 150 | 0.967 | 2.14 | 2.14 | 0.91 | 2.14 | 0.952 |
| iid in time | aligned | whiten | 150 | 0.927 | 2.09 | 2.14 | 0.91 | 2.15 | 0.950 |
| iid in time | orthogonal | none | 150 | 0.613 | 0.83 | 0.65 | 0.91 | 0.92 | 0.946 |
| iid in time | orthogonal | separable | 150 | 0.960 | 0.82 | 0.65 | 0.91 | 0.66 | 0.958 |
| iid in time | orthogonal | whiten | 150 | 0.933 | 0.83 | 0.65 | 0.91 | 0.65 | 0.943 |
| AR(1) 0.5, separable | aligned | separable | 150 | 0.960 | 3.60 | 3.71 | 1.58 | 3.52 | 0.931 |
| AR(1) 0.5, separable | aligned | whiten | 150 | 0.947 | 3.82 | 3.71 | 1.58 | 3.61 | 0.937 |
| AR(1) 0.5, separable | orthogonal | separable | 150 | 0.960 | 1.35 | 1.12 | 1.58 | 1.10 | 0.903 |
| AR(1) 0.5, separable | orthogonal | whiten | 150 | 0.947 | 1.37 | 1.12 | 1.58 | 1.10 | 0.915 |
| factor AR(1) 0.7 + white idiosyncratic | aligned | separable | 100 | 0.700 | 4.74 | 4.91 | 2.18 | 2.48 | 0.729 |
| factor AR(1) 0.7 + white idiosyncratic | aligned | full | 100 | 0.880 | 4.66 | 4.91 | 2.18 | 4.01 | 0.898 |
| factor AR(1) 0.7 + white idiosyncratic | orthogonal | separable | 100 | 0.500 | 0.74 | 0.65 | 2.18 | 1.28 | 1.000 |
| factor AR(1) 0.7 + white idiosyncratic | orthogonal | full | 100 | 0.840 | 0.85 | 0.65 | 2.18 | 0.68 | 0.929 |

For whiten the intervals are in the whitened coordinates (predicted sd in the whitened panel is 2.15 and 0.65 for aligned and orthogonal jumps at phi = 0, 3.74 and 1.12 at phi = 0.5).

## 44. Estimating the top long-run eigenvalue under a persistent common factor (N = 10, T = 300, rho = 0.5, factor AR(1) 0.7 with white idiosyncratic noise, 200 replications; ratios to the truth, top eigenvalue 28.83, idiosyncratic variance 0.5)

| estimator | median | 10th percentile | 90th percentile | mean |
|---|---|---|---|---|
| full Parzen matrix, top eigenvalue (L = 6) | 0.57 | 0.43 | 0.70 | 0.57 |
| PC score, Parzen L=12 | 0.70 | 0.50 | 0.92 | 0.72 |
| PC score, AR first differences | 0.87 | 0.39 | 1.85 | 1.04 |
| PC score, max of the two | 0.91 | 0.62 | 1.85 | 1.13 |
| bulk pooled, AR first differences | 1.04 | 0.85 | 1.30 | 1.06 |
| bulk pooled, Parzen L=8 | 1.04 | 0.96 | 1.12 | 1.04 |

## 45. Null false-positive rate of the principal-direction-wise option (cross = "pc"), 100 replications, T = 300, N = 10 (compare sections 41; fits capped at 10 knots)

| setting | variant | false positive | mean knots |
|---|---|---|---|
| uncorrelated (control) | pc (no factor found) | 0.050 | 0.07 |
| equicorr 0.5, iid in time | pc, no floor | 0.100 | 0.15 |
| equicorr 0.5, iid in time; pc with floor | pc, with floor | 0.020 | 0.02 |
| equicorr 0.5, AR(1) 0.5 separable; pc with floor | pc, with floor | 0.020 | 0.02 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable); pc_lrv=ar | pc, no floor, factor variance from AR | 0.060 | 0.08 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable); pc_lrv=max | pc, no floor, max with residual Parzen | 0.000 | 0.00 |
| nonsep pc scale 0.5 | pc, no floor, factor variances x 0.5 | 0.250 | 0.50 |
| factor AR(1) 0.7 + white idiosyncratic (non-separable); pc with floor | pc, with floor | 0.050 | 0.07 |

## 46. Power in the non-separable design (factor AR(1) 0.7 + white idiosyncratic, equicorrelation 0.5, two kinks at 100 and 200, slope 0.05, 40 replications) and the oracle that uses the true long-run covariance

| jumps | method | both kinks found | missed window | mean knots |
|---|---|---|---|---|
| aligned | full | 0.800 | 0.050 | 2.17 |
| aligned | pc | 0.500 | 0.475 | 1.65 |
| aligned | pc (pc_lrv = max) | 0.525 | 0.450 | 1.62 |
| aligned | pc with floor | 0.500 | 0.475 | 1.65 |
| orthogonal | full | 0.700 | 0.000 | 2.52 |
| orthogonal | pc | 0.400 | 0.425 | 1.82 |
| orthogonal | pc with floor | 0.400 | 0.425 | 1.82 |
| aligned | oracle (penalty 256.8) | 0.700 | 0.300 | 1.73 |
| orthogonal | oracle (penalty 256.8) | 0.525 | 0.450 | 1.57 |

Oracle null false-positive rate (100 replications): 0.000, so the oracle penalty is conservative here. The pc rows without a floor were run before the floor existed; the floor does not change them in this design.

## 47. One kink at 150, slope 0.08, factor AR(1) 0.7 + white idiosyncratic (non-separable), 100 replications: localization and 95% intervals (predicted sd with the true long-run covariance)

| jumps | cross | one knot | empirical sd | predicted sd | mean estimated sd | coverage |
|---|---|---|---|---|---|---|
| aligned | separable | 0.700 | 4.74 | 4.91 | 2.48 | 0.729 |
| aligned | full | 0.880 | 4.66 | 4.91 | 4.01 | 0.898 |
| aligned | pc | 0.980 | 4.55 | 4.91 | 5.63 | 0.898 |
| orthogonal | separable | 0.500 | 0.74 | 0.65 | 1.28 | 1.000 |
| orthogonal | full | 0.840 | 0.85 | 0.65 | 0.68 | 0.929 |
| orthogonal | pc | 0.930 | 0.83 | 0.65 | 0.76 | 0.935 |
## 48. Calibration of the weighted Monte Carlo penalty: true 95% quantile of the null maximum statistic (simulated from the dependent model) against the i.i.d.-in-time Monte Carlo quantile for the same long-run covariance (N = 10, gap 6; 2000 simulations at T = 300, 1000 at T = 1000)

| model | T | true quantile | i.i.d. MC quantile | ratio |
|---|---|---|---|---|
| iid | 300 | 25.0 | 25.5 | 0.982 |
| AR(1) 0.3 | 300 | 44.8 | 47.3 | 0.947 |
| AR(1) 0.5 | 300 | 70.9 | 76.5 | 0.927 |
| AR(1) 0.7 | 300 | 127.4 | 144.4 | 0.882 |
| AR(1) 0.8 | 300 | 192.6 | 229.4 | 0.839 |
| AR(1) 0.9 | 300 | 362.2 | 484.3 | 0.748 |
| equicorr 0.5, iid time | 300 | 47.1 | 45.7 | 1.030 |
| equicorr 0.5, AR(1) 0.5 (separable) | 300 | 128.7 | 137.2 | 0.938 |
| factor AR(1) 0.7 + white | 300 | 190.2 | 220.3 | 0.864 |
| iid | 1000 | 25.4 | 25.8 | 0.987 |
| AR(1) 0.3 | 1000 | 45.6 | 47.9 | 0.951 |
| AR(1) 0.5 | 1000 | 72.5 | 77.3 | 0.937 |
| AR(1) 0.7 | 1000 | 134.6 | 146.1 | 0.921 |
| AR(1) 0.8 | 1000 | 211.9 | 232.0 | 0.913 |
| AR(1) 0.9 | 1000 | 421.4 | 489.9 | 0.860 |
| equicorr 0.5, iid time | 1000 | 48.2 | 49.3 | 0.977 |
| equicorr 0.5, AR(1) 0.5 (separable) | 1000 | 135.3 | 148.0 | 0.914 |
| factor AR(1) 0.7 + white | 1000 | 206.6 | 238.1 | 0.868 |

The i.i.d. MC itself reads 0.98 against truth at phi = 0 (simulation noise of about 2%).

## 49. Exact finite-sample variance ratios v'Gamma v / omega^2 for AR(1) errors at T = 300 (unit residualized hinge and step at the middle) against the observed quantile ratios

| phi | hinge ratio | step ratio | observed quantile ratio |
|---|---|---|---|
| 0 | 1.000 | 1.000 | 0.982 |
| 0.3 | 0.993 | 0.980 | 0.947 |
| 0.5 | 0.986 | 0.961 | 0.927 |
| 0.7 | 0.971 | 0.920 | 0.882 |
| 0.8 | 0.951 | 0.872 | 0.839 |
| 0.9 | 0.886 | 0.741 | 0.748 |

## 50. Knot-free estimate as effective variance (omega^2 times the step ratio from the fitted autocovariances): null false-positive rate and power for ar_start and ar (T = 300, N = 10; 150 null replications, 60 power replications; 'before' = raw long-run variance, sections 39 and 40)

| errors | mode | false positive (after) | false positive (before) | omega2 estimate | true omega2 |
|---|---|---|---|---|---|
| AR(1) 0.5 | ar_start | 0.060 | 0.04 | 2.76 | 3.00 |
| AR(1) 0.5 | ar | 0.067 | 0.06 | 2.98 | 3.00 |
| AR(1) 0.7 | ar_start | 0.060 | 0.03 | 5.01 | 5.67 |
| AR(1) 0.7 | ar | 0.067 | 0.03 | 5.47 | 5.67 |
| AR(1) 0.8 | ar_start | 0.047 | 0.03 | 7.74 | 9.00 |
| AR(1) 0.8 | ar | 0.073 | 0.04 | 8.13 | 9.00 |
| AR(1) 0.9 | ar_start | 0.040 | 0.00 | 14.56 | 19.00 |
| AR(1) 0.9 | ar | 0.133 | 0.04 | 14.26 | 19.00 |

| errors | slope | mode | both kinks found | missed window |
|---|---|---|---|---|
| AR1 0.7 | 0.06 | ar_start | 0.967 | 0.000 |
| AR1 0.7 | 0.06 | ar | 0.967 | 0.000 |
| AR1 0.7 | 0.04 | ar_start | 0.850 | 0.100 |
| AR1 0.7 | 0.04 | ar | 0.800 | 0.117 |
| AR1 0.8 | 0.06 | ar_start | 0.917 | 0.067 |
| AR1 0.8 | 0.06 | ar | 0.933 | 0.067 |

Before (40 replications): phi = 0.7 slope 0.06: 0.975 / 0.95; slope 0.04: 0.875 / 0.80; phi = 0.8 slope 0.06: 0.90 / 0.875 (ar_start / ar).

## 51. Corrected pc option in the non-separable design (factor AR(1) 0.7 + white, equicorrelation 0.5)

| quantity | value | replications | before the correction |
|---|---|---|---|
| null false positive, 500 penalty draws | 0.050 | 100 | 0.05 |
| null false positive, 2000 penalty draws | 0.067 | 45 | - |
| both kinks found, aligned, slope 0.05 | 0.575 | 40 | 0.50 |
| both kinks found, orthogonal, slope 0.05 | 0.425 | 40 | 0.40 |

## 52. Monte Carlo precision of the weighted penalty (N = 10, T = 300; ratios to a 20000-draw reference; 20 seeds at 500 draws, 12 seeds at 2000)

| eigenvalues | reference | 500 draws: mean | sd | min | max | 2000 draws: mean | sd |
|---|---|---|---|---|---|---|---|
| unit | 25.2 | 0.984 | 0.026 | 0.948 | 1.035 | 0.992 | 0.012 |
| oracle_factor | 221.4 | 0.991 | 0.062 | 0.907 | 1.134 | 0.994 | 0.017 |

The default penalty for the oracle eigenvalues at 500 draws (seed 20240101) was 256.8, 16% above the 221.4 reference.

## 53. Oracle benchmark at the true 95% quantile of the null maximum (penalty 190.2; non-separable design, slope 0.05; replaces the oracle rows of section 46, which used 256.8)

| jumps | reps | both kinks found | missed window | false positive (null) |
|---|---|---|---|---|
| - | 200 | - | - | 0.030 |
| aligned | 60 | 0.833 | 0.167 | - |
| orthogonal | 60 | 0.917 | 0.050 | - |

## 54. First-pass overshoot: the trapped panel (non-separable design, orthogonal jumps at slope 0.05, seed 41002; 500 penalty draws)

| pc_lrv | lrv_extra | penalties by pass | knots found |
|---|---|---|---|
| ar | 0 | 444, 468 | 156 |
| resid | 0 | 444, 240, 236 | 99, 200 |
| resid | 2 | 444, 231, 229, 236 | 99, 200 |

## 55. Null false-positive rate of the residual-based later-pass option (cross = "pc", 500 penalty draws, T = 300, N = 10)

| errors | pc_lrv | lrv_extra | reps | false positive |
|---|---|---|---|---|
| factor AR(1) 0.7 + white (non-separable) | ar (corrected, section 51) | 0 | 100 | 0.050 |
| factor AR(1) 0.7 + white (non-separable) | resid | 0 | 100 | 0.010 |
| factor AR(1) 0.7 + white (non-separable) | resid | 2 | 100 | 0.030 |
| equicorr 0.5, iid in time | resid | 2 | 80 | 0.037 |
| uncorrelated | resid | 2 | 80 | 0.037 |
| equicorr 0.5, AR(1) 0.5 separable | resid | 2 | 80 | 0.025 |
| factor AR(1) 0.5 + white (non-separable) | resid | 2 | 80 | 0.037 |

## 56. Power with the residual-based later-pass option (two kinks at 100 and 200, slope 0.05, N = 10, T = 300; non-separable designs with equicorrelation 0.5)

| factor AR(1) | jumps | pc_lrv | lrv_extra | reps | both kinks found | missed window |
|---|---|---|---|---|---|---|
| 0.7 | aligned | ar (corrected, section 51) | 0 | 40 | 0.575 | 0.400 |
| 0.7 | orthogonal | ar (corrected, section 51) | 0 | 40 | 0.425 | 0.375 |
| 0.7 | aligned | resid | 0 | 60 | 0.567 | 0.367 |
| 0.7 | orthogonal | resid | 0 | 60 | 0.717 | 0.200 |
| 0.7 | aligned | resid | 2 | 60 | 0.750 | 0.150 |
| 0.7 | orthogonal | resid | 2 | 60 | 0.783 | 0.117 |
| 0.5 | aligned | ar | 0 | 60 | 0.800 | 0.117 |
| 0.5 | aligned | resid | 2 | 60 | 0.917 | 0.000 |
| 0.7 | aligned | oracle at the true quantile (190.2), section 53 | - | 60 | 0.833 | 0.167 |
| 0.7 | orthogonal | oracle at the true quantile (190.2), section 53 | - | 60 | 0.917 | 0.050 |

## 57. One kink at 150, slope 0.08, factor AR(1) 0.7 + white, 100 replications, resid with two extra knots: localization and 95% intervals (compare section 47)

| jumps | one knot | empirical sd | predicted sd | mean estimated sd | coverage |
|---|---|---|---|---|---|
| aligned | 0.980 | 4.52 | 4.91 | 4.27 | 0.918 |
| orthogonal | 0.930 | 0.84 | 0.65 | 0.73 | 0.935 |
