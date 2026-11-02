library(panelart)
ok <- function(cond, msg) if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)

# 1. The pooled criterion equals the reduction in the sum of squares
set.seed(3)
Y <- matrix(rnorm(60 * 4), 60, 4)
S0 <- panelart:::.cumsum0(Y)
L <- 5L; R <- 50L; ks <- 10:40
V <- panelart:::.node_V(S0, L, R, ks)
sse <- function(M) sum(sweep(M, 2, colMeans(M))^2)
bf <- sapply(ks, function(k) sse(Y[(L + 1):R, ]) - sse(Y[(L + 1):k, ]) - sse(Y[(k + 1):R, ]))
ok(max(abs(V - bf)) < 1e-10, "criterion equals SSE reduction")

# 2. Penalties
ok(abs(art_penalty_iid(10, 200, 2) - 2 * art_penalty_iid(10, 200, 1)) < 1e-9, "penalty scales with sigma2")
ok(art_penalty_mc(5, 100, 1, nsim = 100) < art_penalty_iid(5, 100, 1), "mc penalty below union bound")

# 3. Break recovery, iid
set.seed(10)
sim <- art_sim_data(T = 200, N = 10, breaks = c(0.4, 0.7), delta_bar = 0.8)
for (pen in c("iid", "mc")) {
  fit <- art_fit(sim$Y, pen, nsim = 200)
  ok(length(fit$dates) == 2 && all(abs(fit$dates - sim$dates) <= 3), paste("iid recovery with", pen))
}
cl <- art_classic(sim$Y, "bic")
ok(length(cl$dates) >= 2, "classic tree finds breaks")

# 4. No false positives on a null panel (iid penalty is conservative)
set.seed(11)
Y0 <- matrix(rnorm(150 * 8), 150, 8)
ok(length(art_fit(Y0, "iid")$dates) == 0, "no breaks in null panel")

# 5. Long-run covariance: trace near N * omega^2 for AR(1)
set.seed(12)
sa <- art_sim_data(T = 600, N = 5, breaks = numeric(0), errors = "ar1", phi = 0.5)
om <- art_lrcov(sa$Y, integer(0), "auto")
ok(om$tr / 5 > 2.5 && om$tr / 5 < 5.5, "long-run variance in plausible range")

# 6. Dependent-error fit runs and reports passes
set.seed(13)
sd2 <- art_sim_data(T = 250, N = 6, breaks = c(0.5), delta_bar = 1, errors = "ar1", phi = 0.4)
fh <- art_fit(sd2$Y, "hac_mc", nsim_cov = 50)
ok(inherits(fh, "art") && !is.null(fh$history), "hac_mc fit object")
ok(length(fh$dates) >= 1, "hac_mc finds the break")

# 7. Whitening and confidence intervals run
fw <- art_fit(sim$Y, "iid", whiten = "diff")
ok(inherits(fw, "art") && fw$whitened, "whitened fit")
ci <- art_ci(art_fit(sim$Y, "mc", nsim = 100))
ok(nrow(ci) == 2 && all(ci$upper > ci$lower), "confidence intervals")

# 8. Limit-law density integrates to one
ok(abs(integrate(art_xi_density, -Inf, Inf)$value - 1) < 1e-6, "xi density integrates to 1")

# 9. Factor projection runs
set.seed(14)
sf <- art_sim_data(T = 200, N = 12, breaks = c(0.5), delta_bar = 1, errors = "factor",
                   r = 1, jump_in = "orthogonal")
ff <- art_fit_factor(sf$Y, r = 1, nsim_cov = 50)
ok(inherits(ff, "art_factor"), "factor projection object")

# 10. Evaluation helper
ev <- art_eval(c(78, 141), c(80, 140), 200)
ok(ev[["found"]] == 2 && ev[["spurious"]] == 0, "art_eval")
# 11. Trend criterion equals the reduction in the sum of squares of per-segment line fits
set.seed(21)
Yt <- matrix(rnorm(60 * 4), 60, 4) + outer(1:60, c(.01, -.02, .03, 0))
S0t <- panelart:::.cumsum0(Yt, "trend"); Lt <- 5L; Rt <- 50L; kt <- 10:40
Vt <- panelart:::.node_V(S0t, Lt, Rt, kt)
sse_line <- function(M) { tt <- seq_len(nrow(M)); sum(sapply(seq_len(ncol(M)), function(j) sum(resid(lm(M[, j] ~ tt))^2))) }
bft <- sapply(kt, function(k) sse_line(Yt[(Lt + 1):Rt, ]) - sse_line(Yt[(Lt + 1):k, ]) - sse_line(Yt[(k + 1):Rt, ]))
ok(max(abs(Vt - bft)) < 1e-8, "trend criterion equals line-fit SSE reduction")

# 12. Trend breaks: a slope break is found as one break; a common trend alone is not a break
set.seed(22)
st <- art_sim_data(T = 200, N = 10, breaks = 0.5, delta_bar = 0, slope_bar = 0.04)
ft <- art_fit(st$Y, "mc", type = "trend", nsim = 200)
ok(ft$n_breaks == 1 && identical(ft$type, "trend"), "trend fit finds one slope break")
nb <- sapply(1:5, function(i) {
  Yn <- matrix(rnorm(200 * 10), 200, 10) + outer(1:200, rep(c(0.01, -0.01), 5))
  art_fit(Yn, "mc", type = "trend", level = 0.01, nsim = 200)$n_breaks
})
ok(mean(nb == 0) >= 0.8, "mostly no break for a pure common trend")
ok(inherits(try(art_ci(ft), silent = TRUE), "try-error"), "no intervals for trend fits")

# 13. Join-point (kink) estimator and classic ART with line fits
set.seed(23)
sk <- art_sim_data(T = 200, N = 10, breaks = 0.5, delta_bar = 0, slope_bar = 0.04)
fk <- art_kink(sk$Y, nsim = 200)
ok(fk$n_breaks == 1 && abs(fk$dates - 100) <= 8, "kink estimator finds one knot near the truth")
cik <- art_ci(fk)
ok(nrow(cik) == 1 && cik$lower < fk$dates && cik$upper > fk$dates && cik$sd > 0, "intervals for kink fits")
ct <- art_classic(sk$Y, "bic", type = "trend")
ok(is.integer(ct$dates), "classic ART with line fits runs")
# the kink criterion is invariant to a series-specific linear trend
Yl <- sk$Y + outer(1:200, seq(-0.02, 0.02, length.out = 10))
ok(identical(art_kink(Yl, nsim = 200)$dates, fk$dates) || abs(art_kink(Yl, nsim = 200)$dates - fk$dates) <= 1, "kink fit invariant to linear trend")

# 14. Kink search: start and cap options; the result does not depend on a perturbed start
set.seed(24)
s2 <- art_sim_data(T = 300, N = 10, breaks = c(1/3, 2/3), delta_bar = 0, slope_bar = 0.04)
fa <- art_kink(s2$Y, nsim = 200, start = s2$dates, max_knots = 2)
fb <- art_kink(s2$Y, nsim = 200, start = s2$dates + c(15, -15), max_knots = 2)
ok(fa$n_breaks == 2 && identical(fa$dates, fb$dates), "repartition reaches the same knots from a perturbed start")
fc <- art_kink(s2$Y, nsim = 200, max_knots = 1)
ok(fc$n_breaks <= 1, "knot cap respected")

# 15. Thinning only removes knots, and leaves one-knot fits alone
ft1 <- art_kink(s2$Y, nsim = 200); ft0 <- art_kink(s2$Y, nsim = 200, thin = FALSE)
ok(length(ft1$dates) <= length(ft0$dates) && identical(ft0$dates, ft1$dates_before_thin), "thinning only removes knots")
f1a <- art_kink(sk$Y, nsim = 200); f1b <- art_kink(sk$Y, nsim = 200, thin = FALSE)
ok(identical(f1a$dates, f1b$dates), "thinning does not act on a single knot")

# 16. Removal-cost identity for knots: cost = ||gamma_hat||^2 * ||h_tilde||^2
fid <- art_kink(s2$Y, nsim = 200, thin = FALSE)
if (length(fid$dates) >= 2L) {
  tt16 <- seq_len(nrow(s2$Y)); Bf <- function(K) cbind(1, tt16, sapply(K, function(k) pmax(tt16 - k, 0)))
  sf <- function(K) sum(qr.resid(qr(Bf(K)), s2$Y)^2)
  Kf <- fid$dates
  cst <- sf(Kf[-1]) - sf(Kf)
  hj <- qr.resid(qr(Bf(Kf[-1])), pmax(tt16 - Kf[1], 0))
  ok(abs(cst - sum(fid$gamma[1, ]^2) * sum(hj^2)) < 1e-6 * max(1, cst), "removal-cost identity for knots")
}

# 17. Pair test merges a duplicate pair around one strong kink and leaves a distant knot alone
set.seed(17); tt17 <- 1:300; g17 <- 0.1 * sample(c(-1, 1), 10, TRUE)
Y17 <- outer(pmax(tt17 - 150, 0), g17) + matrix(rnorm(3000), 300, 10)
f17 <- art_kink(Y17, nsim = 200, thin = FALSE, start = c(140, 160), pair_test = TRUE)
ok(length(f17$dates) == 1L && abs(f17$dates - 150) <= 6, "pair test merges a duplicate pair into one knot")
f17b <- art_kink(Y17, nsim = 200, thin = FALSE, start = c(140, 160))
ok(length(f17b$dates) == 2L, "without the pair test the duplicate pair is kept (control)")

# 18. Long-run variance penalty for kinks: iid panel, strong kink -> one knot, omega2 near 1; intervals run
set.seed(18); tt18 <- 1:300; g18 <- 0.08 * sample(c(-1, 1), 10, TRUE)
Y18 <- outer(pmax(tt18 - 150, 0), g18) + matrix(rnorm(3000), 300, 10)
f18 <- art_kink(Y18, "hac_mc", nsim = 200)
ok(any(abs(f18$dates - 150) <= 8) && length(f18$dates) <= 2L, "hac_mc finds the kink in iid errors")
ok(f18$omega2 > 0.6 && f18$omega2 < 1.6, "hac_mc long-run variance near 1 for iid errors")
ci18 <- art_ci(f18); j18 <- which.min(abs(f18$dates - 150)); ok(nrow(ci18) == length(f18$dates) && ci18$lower[j18] < f18$dates[j18] && ci18$upper[j18] > f18$dates[j18], "interval from long-run variance")

# 19. AR(1) null: iid-calibrated penalty over-fits, hac_mc does not
set.seed(19)
E19 <- apply(matrix(rnorm(3000), 300, 10), 2, function(z) as.numeric(stats::filter(z, 0.6, "recursive"))) * sqrt(1 - 0.36)
ok(length(art_kink(E19, "mc", nsim = 200, max_knots = 10)$dates) >= 3L, "iid-calibrated penalty over-fits AR(1) noise")
f19 <- art_kink(E19, "hac_mc", nsim = 200)
ok(length(f19$dates) == 0L, "hac_mc finds no knot in AR(1) noise")
ok(f19$omega2 > 1.5, "long-run variance estimate exceeds the marginal variance under positive autocorrelation")

# 20. Series-specific scaling returns the standardized panel and its factors
f20 <- art_kink(Y18, "hac_mc", nsim = 200, scale = "series")
ok(length(f20$scale_factors) == 10L && f20$omega2 == 1 && any(abs(f20$dates - 150) <= 8), "series scaling")

# 21. lrv_extra option runs, is recorded in the history, and keeps the final long-run variance from the selected knots
f21 <- art_kink(Y18, "hac_mc", nsim = 200, lrv_extra = 2L)
ok("extra" %in% names(f21$history) && any(abs(f21$dates - 150) <= 8), "lrv_extra option")

# 22. Knot-free autoregressive long-run variance from first differences: accurate for white and AR(1) noise
set.seed(22); pw <- matrix(rnorm(3000), 300, 10)
ok(abs(panelart:::.lrv_ar_diff(pw) - 1) < 0.4, "AR long-run variance near 1 for white noise")
pa <- apply(matrix(rnorm(3000), 300, 10), 2, function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
ok(abs(panelart:::.lrv_ar_diff(pa) / 4 - 1) < 0.4, "AR long-run variance near 4 for AR(1) 0.5 noise")
pk <- pa + outer(pmax(1:300 - 100, 0), rep(0.06, 10)) + outer(pmax(1:300 - 200, 0), rep(0.06, 10))
ok(abs(panelart:::.lrv_ar_diff(pk) / panelart:::.lrv_ar_diff(pa) - 1) < 0.1, "AR long-run variance not affected by kinks")

# 23. The same-sign trap: residual-based first pass settles at one knot, knot-free first pass finds both
set.seed(77001); tt23 <- 1:300; g23 <- 0.06 * sample(c(-1, 1), 10, TRUE)
E23 <- apply(matrix(rnorm(3000), 300, 10), 2, function(z) as.numeric(stats::filter(z, 0.7, "recursive"))) * sqrt(1 - 0.49)
Y23 <- outer(pmax(tt23 - 100, 0), g23) + outer(pmax(tt23 - 200, 0), g23) + E23
ok(length(art_kink(Y23, "hac_mc", nsim = 200, lrv = "hac")$dates) == 1L, "residual-based first pass settles at one knot (trap)")
f23 <- art_kink(Y23, "hac_mc", nsim = 200, lrv = "ar")$dates
ok(length(f23) == 2L && all(abs(f23 - c(100, 200)) <= 12), "knot-free long-run variance finds both kinks")

# 24. Weighted Monte Carlo penalty: equals the unit penalty for unit weights, larger for spiked weights
Nw <- 10; lam_spike <- c(5.5, rep(0.5, 9))
pu <- art_penalty_kink_mc(Nw, 300, 0.05, 300L, 6L); pw1 <- art_penalty_kink_mc_w(rep(1, Nw), 300, 0.05, 300L, 6L)
ok(abs(pu - pw1) < 1e-8, "weighted penalty with unit weights equals the unit penalty")
ok(art_penalty_kink_mc_w(lam_spike, 300, 0.05, 300L, 6L) > 1.5 * pu, "spiked weights raise the penalty")

# 25. Cross-sectional options: separable returns Omega, intervals use the sandwich; whiten returns W
set.seed(25); Om25 <- 0.5 * diag(10) + 0.5 * matrix(1, 10, 10); E25 <- matrix(rnorm(3000), 300, 10) %*% chol(Om25)
Y25 <- outer(pmax(1:300 - 150, 0), rep(0.08, 10)) + E25
f25 <- art_kink(Y25, "hac_mc", nsim = 200, cross = "separable")
ok(any(abs(f25$dates - 150) <= 10) && is.matrix(f25$Omega) && all(dim(f25$Omega) == 10), "separable fit returns Omega")
ok(abs(sum(diag(f25$Omega)) / (10 * f25$omega2) - 1) < 1e-8, "Omega is scaled to the pooled long-run variance")
j25 <- which.min(abs(f25$dates - 150)); ci25 <- art_ci(f25); f25n <- art_kink(Y25, "hac_mc", nsim = 200)
ci25n <- art_ci(f25n); jn <- which.min(abs(f25n$dates - 150))
ok(ci25$sd[j25] > 1.5 * ci25n$sd[jn], "sandwich interval is wider than the scalar one for jumps along the factor")
w25 <- art_kink(Y25, "hac_mc", nsim = 200, cross = "whiten")
Dw <- diff(w25$Y); Cw <- crossprod(sweep(Dw, 2, colMeans(Dw), "-")) / nrow(Dw)
ok(is.matrix(w25$W) && max(abs(Cw / mean(diag(Cw)) - diag(10))) < 1e-6, "whitened panel has identity first-difference covariance")
fu25 <- art_kink(Y25, "hac_mc", nsim = 200, cross = "full")
ok(is.matrix(fu25$Omega) && any(abs(fu25$dates - 150) <= 10), "full long-run covariance option runs")
ok(inherits(try(art_kink(Y25, "hac_mc", nsim = 100, cross = "separable", scale = "series"), silent = TRUE), "try-error"), "cross with series scaling is refused")

# 26. Principal-direction-wise option: finds one factor in a factor panel, none in an uncorrelated one, returns Omega
set.seed(26); tt26 <- 1:300
f26 <- as.numeric(stats::filter(rnorm(300), 0.7, "recursive", init = 0)) * sqrt(1 - 0.49) * sqrt(0.5)
E26 <- f26 %*% t(rep(1, 10)) + matrix(rnorm(3000), 300, 10) * sqrt(0.5)
Y26 <- outer(pmax(tt26 - 150, 0), rep(0.08, 10)) + E26
p26 <- art_kink(Y26, "hac_mc", nsim = 200, cross = "pc")
ok(p26$n_factors == 1L && is.matrix(p26$Omega) && any(abs(p26$dates - 150) <= 10), "pc finds one factor and the kink")
ev26 <- eigen(p26$Omega, symmetric = TRUE)$values
ok(ev26[1] > 5 * ev26[2], "pc long-run covariance has a dominant factor eigenvalue")
set.seed(27); u27 <- art_kink(matrix(rnorm(3000), 300, 10), "hac_mc", nsim = 200, cross = "pc")
ok(u27$n_factors == 0L, "pc finds no factor in uncorrelated series")
ok(!is.null(art_ci(p26)), "pc intervals run")
ok(art_kink(Y26, "hac_mc", nsim = 200, cross = "pc", n_factors = 0L)$n_factors == 0L, "n_factors override")
ok(sum(diag(art_kink(Y26, "hac_mc", nsim = 200, cross = "pc", pc_floor = FALSE)$Omega)) > 0, "pc without the floor runs")

# 27. Finite-sample effective variance: below the long-run variance for positive autocorrelation, equal for white noise
set.seed(28); Pw <- matrix(rnorm(3000), 300, 10)
P8 <- apply(matrix(rnorm(3000), 300, 10), 2, function(z) as.numeric(stats::filter(z, 0.8, "recursive", init = 0))) * sqrt(1 - 0.64)
rw <- panelart:::.lrv_ar_diff(Pw, effective = TRUE) / panelart:::.lrv_ar_diff(Pw, effective = FALSE)
r8 <- panelart:::.lrv_ar_diff(P8, effective = TRUE) / panelart:::.lrv_ar_diff(P8, effective = FALSE)
ok(abs(rw - 1) < 0.02, "effective variance equals the long-run variance for white noise")
ok(r8 > 0.80 && r8 < 0.94, "effective variance is about 0.87 of the long-run variance at AR(1) 0.8, T = 300")
sv <- panelart:::.step_vec(300)
ok(abs(sum(sv^2) - 1) < 1e-10 && abs(sum(sv)) < 1e-8 && abs(sum(sv * (1:300))) < 1e-6, "step vector is a unit vector orthogonal to constant and trend")
# exact ratio for AR(1) 0.8 from the Toeplitz covariance
Gm <- toeplitz(0.8^(0:299)); ex <- drop(crossprod(sv, Gm %*% sv)) / ((1 + 0.8) / (1 - 0.8))
ok(abs(ex - 0.872) < 0.01, "exact step ratio at AR(1) 0.8, T = 300 is 0.872")

# 28. pc_lrv = "resid": a first-pass overshoot of the factor variance is not locked in
N28 <- 10; T28 <- 300; tt28 <- 1:T28
ar28 <- function(M, p) apply(M, 2, function(z) as.numeric(stats::filter(z, p, "recursive", init = 0))) * sqrt(1 - p^2)
set.seed(41002); g28 <- 0.05 * rep(c(1, -1), length.out = N28)
mu28 <- outer(pmax(tt28 - 100, 0), g28) + outer(pmax(tt28 - 200, 0), g28)
f28 <- ar28(matrix(rnorm(T28), T28, 1), 0.7) * sqrt(0.5)
Y28 <- mu28 + f28 %*% t(rep(1, N28)) + matrix(rnorm(T28 * N28), T28, N28) * sqrt(0.5)
a28 <- art_kink(Y28, "hac_mc", nsim = 500, cross = "pc", pc_lrv = "ar", max_knots = 10)
r28 <- art_kink(Y28, "hac_mc", nsim = 500, cross = "pc", pc_lrv = "resid", lrv_extra = 2L, max_knots = 10)
ok(length(a28$dates) == 1L, "pc_lrv = ar settles at one knot on the trapped panel")
ok(length(r28$dates) == 2L && all(abs(r28$dates - c(100, 200)) <= 12), "pc_lrv = resid finds both kinks on the trapped panel")
ok("extra" %in% names(r28$history) && r28$history$alpha[2] < 0.75 * r28$history$alpha[1], "resid lowers the penalty after the first pass")
ok(length(panelart:::.augment_knots(Y28, c(100L, 200L), 2L, 6L)) == 4L && length(panelart:::.augment_knots(Y28, c(100L), 0L, 6L)) == 1L, "augment_knots adds the requested number of knots")

# 29. Default settings of cross = "pc": residual-based later passes with two extra knots
d29 <- art_kink(Y28, "hac_mc", nsim = 200, cross = "pc", max_knots = 10)
ok(d29$history$extra[1] == 2L, "cross = pc defaults to lrv_extra = 2")
ok(art_kink(Y28, "hac_mc", nsim = 200, cross = "pc", lrv_extra = 0L, max_knots = 10)$history$extra[1] == 0L, "explicit lrv_extra is respected")

cat("all tests passed\n")
