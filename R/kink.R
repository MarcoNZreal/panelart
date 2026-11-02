#' Monte Carlo penalty for the join-point (kink) criterion
#'
#' Simulates N(0,1) panels without breaks, regresses each series on a constant
#' and a time trend, and returns the \code{1 - level} quantile of the largest
#' gain in the sum of squares from adding one kink (a hinge regressor
#' \code{(t - k)_+}) over knots at least \code{gap} from the ends. Cached.
#'
#' @param N,T Panel dimensions.
#' @param level Nominal null rejection probability.
#' @param nsim Number of simulated panels. The null level of the search is close
#'   to \code{level} only for a large number (200 gave 9.5\% false positives at the
#'   5\% level, 2000 gave 5.3\%).
#' @param gap Minimum distance of a knot from the ends.
#' @param seed Internal seed.
#' @return The penalty for unit error variance.
#' @export
art_penalty_kink_mc <- function(N, T, level = 0.05, nsim = 2000L, gap = 4L,
                                seed = 20240101L) {
  key <- paste("kink", N, T, level, nsim, gap, sep = "_")
  if (!is.null(.pkg_cache[[key]])) return(.pkg_cache[[key]])
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  tt <- seq_len(T); ks <- gap:(T - gap)
  H <- outer(tt, ks, function(t, k) pmax(t - k, 0))
  q <- qr(cbind(1, tt))
  Hres <- qr.resid(q, H); nrm <- colSums(Hres^2)
  mx <- numeric(nsim)
  for (i in seq_len(nsim)) {
    Yres <- qr.resid(q, matrix(rnorm(T * N), T, N))
    mx[i] <- max(rowSums(crossprod(Hres, Yres)^2) / nrm)
  }
  if (had) assign(".Random.seed", old, envir = globalenv())
  else rm(".Random.seed", envir = globalenv())
  val <- unname(quantile(mx, 1 - level))
  .pkg_cache[[key]] <- val
  val
}

#' Common kinks: continuous piecewise-linear breaks in a panel
#'
#' Fits, for every series, a continuous piecewise-linear mean with common knots
#' (join points): \eqn{\mu_{it} = a_i + b_i t + \sum_j \gamma_{ij} (t - k_j)_+}.
#' The knots \code{k_j} are shared by all series; the intercept, slope and slope
#' changes are free for each series. This is the join-point (segmented
#' regression) model; compare \code{art_fit(type = "trend")}, which fits a
#' separate line in each segment and so also allows level jumps.
#'
#' The search mirrors \code{\link{art_fit}}: (1) add the knot with the largest
#' gain in the pooled sum of squares if it exceeds \code{alpha}; (2) move each
#' knot to its best position given the others; (3) drop the knot with the
#' smallest removal cost if it is below \code{alpha}; repeat. Knots stay at
#' least \code{min_gap} apart and from the ends. The penalty is the Monte
#' Carlo null quantile for one added kink (\code{penalty = "mc"}), the
#' union-bound form for independent errors (\code{"iid"}) or a fixed value.
#' Afterwards a thinning step (\code{thin = TRUE}) removes a knot if deleting it
#' and re-optimizing the other knots raises the sum of squares by less than
#' \code{alpha}. The error variance is estimated from second differences. This estimator is
#' not covered by the theory in the notes.
#'
#' @param Y Panel matrix (T x N), or a vector.
#' @param penalty One of \code{"mc"}, \code{"iid"}, \code{"fixed"}, \code{"hac_mc"}.
#'   \code{"hac_mc"} is for serially correlated errors: the Monte Carlo penalty is
#'   multiplied by an estimate of the long-run variance of the errors, obtained
#'   from the residuals of the current fit with a Parzen kernel and iterated until
#'   the knots stop changing (section 18 of the notes). The first pass uses the
#'   bandwidth \code{floor(T^(1/3))}. The long-run variance is returned as
#'   \code{omega2} and used by \code{\link{art_ci}}.
#' @param sigma2 Known error variance (default: second-difference estimate).
#' @param alpha Penalty for \code{penalty = "fixed"}.
#' @param eps Minimum spacing fraction of T between knots.
#' @param min_gap Minimum spacing in dates; default \code{max(3, ceiling(eps T))}.
#' @param level,nsim Null level and simulation size for \code{penalty = "mc"}.
#' @param C0 Constant for \code{penalty = "iid"}.
#' @param max_iter Safety cap on search cycles.
#' @param start Initial knots (default none). The search starts from these.
#' @param thin Logical; apply the thinning step described in Details.
#' @param alpha_thin Threshold for the thinning step. Default \code{NULL} uses
#'   \code{alpha}. The value \code{\"2N\"} uses
#'   \code{sigma2 * (2N + 2 sqrt(2N log T) + 2 log T)}. This acts as a stricter
#'   overall level for extra knots: in simulation it removed isolated
#'   false-positive knots as well as duplicate pairs (section 17 of the notes),
#'   and it deletes weak true knots that sit next to a strong one.
#' @param pair_test Experimental. Logical; after thinning, test every pair of
#'   neighbouring knots closer than \code{merge_dist} against the best single
#'   knot and merge the pair if the test does not reject. In simulation it
#'   merged the same duplicates as thinning at \code{alpha} and no more
#'   (section 17 of the notes). Default \code{FALSE}.
#' @param merge_dist Largest distance between two knots that are tested as a
#'   pair; default \code{ceiling(0.1 T)}.
#' @param pair_level Level of the pair test. The null distribution of the
#'   statistic (best pair minus best single in the window, in the sum of squares)
#'   is obtained by a parametric bootstrap from the fitted single-knot model.
#' @param pair_nboot Number of bootstrap draws for the pair test.
#' @param L Bandwidth of the long-run variance estimate for \code{penalty = "hac_mc"}
#'   after the first pass; default is the Andrews AR(1) plug-in.
#' @param max_pass Largest number of long-run variance iterations.
#' @param lrv_extra For \code{penalty = "hac_mc"} (default 0, or 2 when \code{cross = "pc"} with \code{pc_lrv = "resid"}): number of extra knots, chosen
#'   greedily by gain without any penalty, that are added to the selected knots
#'   before the residuals are used to estimate the long-run variance. They absorb
#'   bends that the current fit has not yet found, which would otherwise inflate
#'   the estimate and can trap the iteration at too few knots, for example when
#'   two same-sign kinks look like one bend (section 18). Setting 2 removed all
#'   missed knots in that case but raised the false-positive rate when a real kink
#'   is present, so it is off by default.
#' @param lrv For \code{penalty = "hac_mc"}: how the long-run variance is
#'   estimated (default \code{"ar_start"}). \code{"hac"} uses the residuals of the current fit in every pass
#'   (first pass: all knots unremoved, bandwidth \code{floor(T^(1/3))}).
#'   \code{"ar_start"} uses an autoregressive estimate from the first differences
#'   in the first pass, which does not depend on the knots, and the residual
#'   estimate afterwards. \code{"ar"} uses the autoregressive estimate in every
#'   pass, so the penalty does not depend on the fitted knots at all. The
#'   autoregressive estimate fits AR(p) models of order 2, 3 and 4 to the
#'   autocovariances of the first differences and takes the geometric mean of the
#'   three long-run variances (section 19 of the notes); it is accurate for
#'   autoregressive errors and within a factor of about 2 for the ARMA and MA
#'   cases tried.
#' @param cross For \code{penalty = "hac_mc"}: treatment of correlation across
#'   series. \code{"none"} assumes uncorrelated series. \code{"separable"}
#'   estimates the cross-sectional covariance shape from the first differences
#'   (which the kinks barely affect), takes its eigenvalues as weights in the
#'   Monte Carlo penalty, and \code{art_ci} then uses a sandwich variance with
#'   the long-run covariance \code{Omega}; the long-run covariance is the
#'   estimated shape scaled to the pooled long-run variance, which assumes all
#'   series share one temporal autocorrelation. \code{"whiten"} multiplies the
#'   panel by the inverse square root of that shape and then proceeds as for
#'   uncorrelated series (the returned \code{Y} is whitened and \code{W} holds
#'   the transformation). \code{"full"} starts from the separable estimate and
#'   then estimates the whole N x N long-run covariance from the residuals of the
#'   current fit with a Parzen kernel in every later pass, which does not require
#'   a common temporal autocorrelation; the penalty is recalculated from its
#'   eigenvalues each pass. \code{"pc"} estimates the long-run variance separately
#'   along each of the leading principal directions of the first-difference
#'   covariance (the common factors, whose loadings the first differences still
#'   identify although they shrink a persistent factor's variance) and pools the
#'   remainder; see \code{n_factors} and \code{pc_lrv}. Use with \code{N} well below \code{T / 10}; see
#'   section 20 of the notes.
#' @param n_factors For \code{cross = "pc"}: number of leading directions treated
#'   as factors. The default chooses it (at most 3) by the eigenvalue ratio of the
#'   first-difference covariance, accepting a factor only if the ratio is at least
#'   1.5; 0 means none.
#' @param pc_lrv For \code{cross = "pc"}: the default \code{"resid"} uses the
#'   knot-free estimate only in the first pass and the residual-based Parzen
#'   estimate afterwards (after removing the fitted knots and \code{lrv_extra}
#'   unpenalized extra knots, which default to 2 in this case), which avoids an
#'   overshoot of the first pass being locked in; \code{"ar"} uses the knot-free
#'   autoregressive estimate (section 19) for each factor direction in every
#'   pass; \code{"max"} uses the larger of that and the residual-based Parzen
#'   estimate after the first pass (section 21 and 23 of the notes).
#' @param pc_scale For \code{cross = "pc"}: multiplier applied to the factor
#'   long-run variances (default 1). Values below 1 trade size for power.
#' @param pc_floor For \code{cross = "pc"} (default \code{TRUE}): use the larger
#'   of the factor-direction estimate and the variance in that direction implied
#'   by the separable estimate, which pools over all series. This protects the
#'   size when the dependence is in fact separable and the single-series factor
#'   estimate is too low.
#' @param scale For \code{penalty = "hac_mc"}: \code{"pooled"} uses one long-run
#'   variance for all series; \code{"series"} divides each series by the square
#'   root of its own long-run variance before fitting, which keeps the null
#'   distribution exact when the series differ in persistence (the pooled version
#'   over-rejects then). The returned fit holds the standardized panel in
#'   \code{Y}, with the factors in \code{scale_factors}, and \code{omega2} is 1.
#' @param max_knots Largest number of knots; no knot is added once it is reached
#'   (moves and deletions still happen).
#' @return An object of class \code{"art"} with \code{type = "kink"}, the knot
#'   dates, the slope changes \code{gamma} (knots x series) and the fitted values.
#' @examples
#' set.seed(1)
#' sim <- art_sim_data(T = 200, N = 10, breaks = 0.5, delta_bar = 0, slope_bar = 0.03)
#' art_kink(sim$Y)
#' @export
art_kink <- function(Y, penalty = c("mc", "iid", "fixed", "hac_mc"), sigma2 = NULL,
                     alpha = NULL, eps = 0.02, min_gap = NULL, level = 0.05,
                     nsim = 2000L, C0 = 10, max_iter = 100L,
                     start = integer(0), max_knots = Inf, thin = TRUE,
                     alpha_thin = NULL, pair_test = FALSE, merge_dist = NULL,
                     pair_level = 0.01, pair_nboot = 200L, L = NULL, max_pass = 8L,
                     scale = c("pooled", "series"), lrv_extra = 0L,
                     lrv = c("ar_start", "hac", "ar"),
                     cross = c("none", "separable", "whiten", "full", "pc"),
                     n_factors = NULL, pc_lrv = c("resid", "ar", "max"), pc_scale = 1, pc_floor = TRUE) {
  cl <- match.call(); penalty <- match.arg(penalty); scale <- match.arg(scale); lrv <- match.arg(lrv); cross <- match.arg(cross); pc_lrv <- match.arg(pc_lrv)
  if (cross == "pc" && pc_lrv == "resid" && missing(lrv_extra)) lrv_extra <- 2L
  Y <- .as_panel(Y); Tn <- nrow(Y); N <- ncol(Y)
  if (Tn < 20L) stop("need at least 20 dates")
  if (penalty == "hac_mc")
    return(.art_kink_hac(Y, cl, sigma2, eps, min_gap, level, nsim, max_iter, start, max_knots,
                         thin, alpha_thin, pair_test, merge_dist, pair_level, pair_nboot, L, max_pass, scale, lrv_extra, lrv, cross, n_factors, pc_lrv, pc_scale, pc_floor))
  gap <- if (is.null(min_gap)) max(3L, as.integer(ceiling(eps * Tn))) else as.integer(min_gap)
  if (is.null(sigma2) && penalty != "fixed") sigma2 <- art_sigma2_diff(Y, 2L)
  alpha <- switch(penalty,
    mc = sigma2 * art_penalty_kink_mc(N, Tn, level, nsim, gap),
    iid = art_penalty_iid(N, Tn, sigma2, C0),
    fixed = { if (is.null(alpha)) stop("alpha required"); alpha })
  tt <- seq_len(Tn)
  Hall <- outer(tt, seq_len(Tn), function(t, k) pmax(t - k, 0))
  basis <- function(K) cbind(1, tt, Hall[, K, drop = FALSE])
  allowed <- function(O) {
    ks <- gap:(Tn - gap)
    for (o in O) ks <- ks[abs(ks - o) >= gap]
    ks
  }
  # gain in sum of squares from adding each candidate knot to the knots O
  gains <- function(O, ks) {
    if (!length(ks)) return(numeric(0))
    q <- qr(basis(O))
    Hres <- qr.resid(q, Hall[, ks, drop = FALSE])
    nrm <- colSums(Hres^2)
    Yres <- qr.resid(q, Y)
    g <- rowSums(crossprod(Hres, Yres)^2) / pmax(nrm, 1e-12)
    g[nrm < 1e-8] <- -Inf
    g
  }
  K <- sort(as.integer(start))
  for (it in seq_len(max_iter)) {
    changed <- FALSE
    ks <- allowed(K); g <- gains(K, ks)
    if (length(g) && max(g) > alpha && length(K) < max_knots) {
      K <- sort(c(K, ks[which.max(g)])); changed <- TRUE
    }
    if (length(K)) {
      for (j in seq_along(K)) {
        O <- K[-j]; ks <- allowed(O); g <- gains(O, ks)
        cur <- g[match(K[j], ks)]
        if (!is.na(cur) && max(g) > cur + 1e-8) { K[j] <- ks[which.max(g)]; changed <- TRUE }
      }
      K <- sort(K)
      repeat {
        if (!length(K)) break
        cost <- vapply(seq_along(K), function(j) gains(K[-j], K[j]), numeric(1))
        j <- which.min(cost)
        if (cost[j] < alpha) { K <- K[-j]; changed <- TRUE } else break
      }
    }
    if (!changed) break
  }
  K0 <- K
  # Thinning: removal cost with re-optimization. For each knot, delete it, let
  # the remaining knots move to their best positions (coordinate ascent), and
  # record the rise in the sum of squares. A knot is dropped if the smallest
  # such cost is below alpha. The neighbour of the deleted knot moves first, so
  # the cost never exceeds the cost of merging the two knots into the best one. This removes pairs of knots that jointly mimic
  # one kink from the same side, a local optimum that plain removal costs miss.
  sse <- function(Kset) sum(qr.resid(qr(basis(Kset)), Y)^2)
  reopt <- function(Kset, first = NULL) {
    for (sw in 1:20) {
      moved <- FALSE
      ord <- if (is.null(first)) seq_along(Kset) else c(first, setdiff(seq_along(Kset), first))
      for (j in ord) {
        O <- Kset[-j]; ks <- allowed(O); g <- gains(O, ks)
        cur <- g[match(Kset[j], ks)]
        if (length(g) && !is.na(cur) && max(g) > cur + 1e-8) { Kset[j] <- ks[which.max(g)]; moved <- TRUE }
      }
      if (!moved) break
    }
    sort(Kset)
  }
  athin <- if (is.null(alpha_thin)) alpha else if (identical(alpha_thin, "2N")) {
    s2t <- if (!is.null(sigma2)) sigma2 else art_sigma2_diff(Y, 2L)
    s2t * (2 * N + 2 * sqrt(2 * N * log(Tn)) + 2 * log(Tn))
  } else as.numeric(alpha_thin)
  if (thin && length(K) >= 2L) {
    repeat {
      if (length(K) < 2L) break
      cur <- sse(K)
      cand <- lapply(seq_along(K), function(j) {
        Kp <- K[-j]
        reopt(Kp, first = which.min(abs(Kp - K[j])))
      })
      cost <- vapply(cand, function(Kp) sse(Kp) - cur, numeric(1))
      j <- which.min(cost)
      if (cost[j] < athin) K <- cand[[j]] else break
    }
  }
  pair_info <- NULL
  if (pair_test && length(K) >= 2L) {
    md <- if (is.null(merge_dist)) as.integer(ceiling(0.1 * Tn)) else as.integer(merge_dist)
    s2p <- if (!is.null(sigma2)) sigma2 else art_sigma2_diff(Y, 2L)
    # statistic for the window around the pair (k1, k2): best two-knot gain
    # minus best one-knot gain, other knots O held fixed
    pair_stat <- function(O, k1, k2, nboot) {
      ks <- allowed(O); ks <- ks[ks >= k1 - md & ks <= k2 + md]
      w <- length(ks); if (w < 3L) return(NULL)
      qO <- qr(basis(O))
      Hres <- qr.resid(qO, Hall[, ks, drop = FALSE])
      G <- crossprod(Hres); nrm <- diag(G)
      dk <- abs(outer(ks, ks, "-"))
      ok <- (dk >= gap) & (dk <= md) & upper.tri(dk)
      den <- outer(nrm, nrm) - G^2
      ok <- ok & (den > 1e-8 * outer(nrm, nrm))
      if (!any(ok)) return(NULL)
      stat <- function(Yres) {
        R <- crossprod(Hres, Yres); S <- tcrossprod(R); Sd <- diag(S)
        g1 <- Sd / pmax(nrm, 1e-12); g1[nrm < 1e-8] <- -Inf
        num <- outer(Sd, nrm) - 2 * G * S + outer(nrm, Sd)
        g2 <- max(num[ok] / den[ok])
        list(T = max(0, g2 - max(g1)), khat = ks[which.max(g1)])
      }
      ob <- stat(qr.resid(qO, Y))
      Fit <- Y - qr.resid(qr(basis(sort(c(O, ob$khat)))), Y)
      had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
      set.seed(20240202L)
      Tb <- vapply(seq_len(nboot), function(b)
        stat(qr.resid(qO, Fit + sqrt(s2p) * matrix(rnorm(Tn * N), Tn, N)))$T, numeric(1))
      if (had) assign(".Random.seed", old, envir = globalenv())
      else rm(".Random.seed", envir = globalenv())
      list(Tobs = ob$T, q = unname(quantile(Tb, 1 - pair_level)), khat = ob$khat)
    }
    pair_info <- list()
    for (iter in 1:10) {
      merged <- FALSE
      for (c in which(diff(K) <= md)) {
        O <- K[-c(c, c + 1L)]
        r <- pair_stat(O, K[c], K[c + 1L], as.integer(pair_nboot))
        if (is.null(r)) next
        pair_info[[length(pair_info) + 1L]] <- c(k1 = K[c], k2 = K[c + 1L], Tobs = r$Tobs,
                                                  q = r$q, merged = as.numeric(r$Tobs < r$q))
        if (r$Tobs < r$q) { K <- sort(c(O, as.integer(r$khat))); merged <- TRUE; break }
      }
      if (!merged || length(K) < 2L) break
    }
    if (length(pair_info)) pair_info <- do.call(rbind, pair_info)
  }
  B <- basis(K); qrB <- qr(B)
  fitted <- Y - qr.resid(qrB, Y)
  coefs <- qr.coef(qrB, Y)
  gamma <- if (length(K)) coefs[-(1:2), , drop = FALSE] else matrix(0, 0, N)
  b <- c(0L, K, Tn)
  structure(list(
    dates = as.integer(K), n_breaks = length(K),
    segments = data.frame(start = b[-length(b)] + 1L, end = b[-1L]),
    means = .segment_means(Y, K), alpha = alpha, penalty = penalty,
    sigma2 = sigma2, omega = NULL, history = NULL, whitened = FALSE, W = NULL,
    Y = Y, Yw = Y, eps = eps, min_len = gap, C0 = C0, L = NULL, call = cl,
    type = "kink", gamma = gamma, fitted = fitted, dates_before_thin = as.integer(K0),
    alpha_thin = athin, pair_info = pair_info),
    class = "art")
}


# Long-run variance of the residuals of the kink fit with knots K, pooled over
# series (Parzen kernel, degrees-of-freedom correction for the 2 + |K| fitted
# regressors).
.kink_lrv <- function(Y, K, L, extra = 0L, gap = 6L) {
  Tn <- nrow(Y); N <- ncol(Y); tt <- seq_len(Tn)
  hing <- function(KK) if (length(KK)) outer(tt, KK, function(t, k) pmax(t - k, 0)) else NULL
  K2 <- K
  if (extra > 0L) {
    Hall <- outer(tt, seq_len(Tn), function(t, k) pmax(t - k, 0))
    for (e in seq_len(extra)) {
      q <- qr(cbind(1, tt, hing(K2)))
      ks <- gap:(Tn - gap); for (o in K2) ks <- ks[abs(ks - o) >= gap]
      if (!length(ks)) break
      Hres <- qr.resid(q, Hall[, ks, drop = FALSE]); nrm <- colSums(Hres^2)
      g <- rowSums(crossprod(Hres, qr.resid(q, Y))^2) / pmax(nrm, 1e-12); g[nrm < 1e-8] <- -Inf
      K2 <- c(K2, ks[which.max(g)])
    }
  }
  B <- cbind(1, tt, hing(K2))
  E <- qr.resid(qr(B), Y)
  d <- colSums(E^2) / Tn
  wsum <- 1
  for (l in seq_len(L)) {
    w <- .kern(l / (L + 1), "parzen")
    g <- colSums(E[1:(Tn - l), , drop = FALSE] * E[(l + 1):Tn, , drop = FALSE])
    d <- d + 2 * w * g / Tn
    wsum <- wsum + 2 * w
  }
  fac <- 1 / max(1 - ncol(B) * wsum / Tn, 0.5)
  list(omega2 = max(mean(d) * fac, 1e-12), omega2_series = pmax(d * fac, 1e-12), E = E)
}

.art_kink_hac <- function(Y, cl, sigma2, eps, min_gap, level, nsim, max_iter, start, max_knots,
                          thin, alpha_thin, pair_test, merge_dist, pair_level, pair_nboot, L, max_pass, scale = "pooled", lrv_extra = 0L, lrv = "hac", cross = "none",
                          n_factors = NULL, pc_lrv = "ar", pc_scale = 1, pc_floor = TRUE) {
  Tn <- nrow(Y); N <- ncol(Y)
  gap <- if (is.null(min_gap)) max(3L, as.integer(ceiling(eps * Tn))) else as.integer(min_gap)
  if (cross != "none" && scale == "series") stop("cross != 'none' cannot be combined with scale = 'series'")
  Cs <- NULL; Wm <- NULL
  if (cross != "none") {
    Cs <- .cross_cov(Y)
    if (cross == "whiten") { Wm <- .inv_sqrt(Cs); Y <- Y %*% Wm }
  }
  pcinfo <- NULL
  if (cross == "pc") {
    eg <- eigen((Cs + t(Cs)) / 2, symmetric = TRUE); rmax <- min(3L, N - 2L)
    r <- if (!is.null(n_factors)) as.integer(n_factors) else {
      ratio <- eg$values[1:rmax] / eg$values[2:(rmax + 1L)]
      kbest <- which.max(ratio); if (ratio[kbest] >= 1.5) kbest else 0L }
    pcinfo <- list(eg = eg, r = r)
  }
  unit <- if (cross %in% c("separable", "full", "pc")) {
    ev <- pmax(eigen(Cs, symmetric = TRUE, only.values = TRUE)$values, 0)
    art_penalty_kink_mc_w(ev * N / sum(ev), Tn, level, nsim, gap)
  } else art_penalty_kink_mc(N, Tn, level, nsim, gap)
  K <- sort(as.integer(start)); Lcur <- max(1L, floor(Tn^(1 / 3)))
  hist <- list(); fit <- NULL; first <- TRUE
  # Phase 1 estimates the long-run variance from residuals that also remove the
  # best lrv_extra unpenalized knots, which gets the iteration out of a too-few-knots
  # fixed point; phase 2 re-estimates it from the selected knots alone.
  for (ex in unique(c(as.integer(lrv_extra), 0L))) {
    seen <- character(0)
    for (pass in 0:max_pass) {
      lr <- NULL
      if (lrv == "ar" || (lrv == "ar_start" && first)) {
        oar <- .lrv_ar_diff(Y)
        if (is.finite(oar)) lr <- list(omega2 = oar, omega2_series = rep(oar, N))
      }
      if (is.null(lr)) lr <- .kink_lrv(Y, K, Lcur, ex, gap)
      om2 <- lr$omega2
      alpha_pass <- om2 * unit; Om_pass <- NULL
      if (cross == "pc") {
        Om_pass <- .pc_omega(Y, K, Lcur, ex, gap, first, Cs, pcinfo, pc_lrv, om2, pc_scale, pc_floor)
        alpha_pass <- art_penalty_kink_mc_w(pmax(eigen(Om_pass, symmetric = TRUE, only.values = TRUE)$values, 0), Tn, level, nsim, gap)
        om2 <- sum(diag(Om_pass)) / N
      }
      if (cross == "full") {
        Om_pass <- if (first) Cs * (N * om2 / sum(diag(Cs))) else .kink_lrcov(Y, K, Lcur, ex, gap)
        alpha_pass <- art_penalty_kink_mc_w(pmax(eigen(Om_pass, symmetric = TRUE, only.values = TRUE)$values, 0), Tn, level, nsim, gap)
        om2 <- sum(diag(Om_pass)) / N
      }
      Yuse <- if (scale == "series") sweep(Y, 2L, sqrt(lr$omega2_series), "/") else Y
      if (scale == "series") om2 <- 1
      fit <- art_kink(Yuse, "fixed", alpha = alpha_pass, sigma2 = om2, eps = eps, min_gap = min_gap,
                      max_iter = max_iter, start = K, max_knots = max_knots, thin = thin,
                      alpha_thin = alpha_thin, pair_test = pair_test, merge_dist = merge_dist,
                      pair_level = pair_level, pair_nboot = pair_nboot)
      Knew <- as.integer(fit$dates)
      hist[[length(hist) + 1L]] <- data.frame(pass = length(hist), L = Lcur, extra = ex, omega2 = om2, alpha = alpha_pass,
                                              n_knots = length(Knew), dates = paste(Knew, collapse = ","))
      conv <- !first && identical(Knew, as.integer(K))
      key <- paste(Knew, collapse = ",")
      K <- Knew; first <- FALSE
      if (conv || key %in% seen) break
      seen <- c(seen, key)
      E <- .kink_lrv(Y, K, 1L, ex, gap)$E
      Lcur <- if (is.null(L)) .auto_bandwidth(E, list(seq_len(Tn)), "parzen") else as.integer(L)
    }
  }
  fin <- NULL
  if (lrv == "ar") { oar <- .lrv_ar_diff(Y); if (is.finite(oar)) fin <- list(omega2 = oar, omega2_series = rep(oar, N)) }
  if (is.null(fin)) fin <- .kink_lrv(Y, K, Lcur, 0L, gap)
  if (scale == "series") {
    fit <- art_kink(sweep(Y, 2L, sqrt(fin$omega2_series), "/"), "fixed", alpha = unit, sigma2 = 1, eps = eps,
                    min_gap = min_gap, max_iter = max_iter, start = K, max_knots = max_knots, thin = FALSE)
    fit$scale_factors <- sqrt(fin$omega2_series); fit$omega2 <- 1; fit$sigma2 <- 1
  } else { fit$omega2 <- fin$omega2; fit$sigma2 <- fin$omega2 }
  if (cross == "separable") fit$Omega <- Cs * (N * fit$omega2 / sum(diag(Cs)))
  if (cross == "pc") { fit$Omega <- .pc_omega(Y, K, Lcur, 0L, gap, FALSE, Cs, pcinfo, pc_lrv, fin$omega2, pc_scale, pc_floor); fit$omega2 <- fit$sigma2 <- sum(diag(fit$Omega)) / N; fit$n_factors <- pcinfo$r }
  if (cross == "full") { fit$Omega <- .kink_lrcov(Y, K, Lcur, 0L, gap); fit$omega2 <- fit$sigma2 <- sum(diag(fit$Omega)) / N }
  if (cross == "whiten") fit$W <- Wm
  fit$penalty <- "hac_mc"; fit$history <- do.call(rbind, hist); fit$call <- cl
  fit$L <- Lcur
  fit
}


# ---- autoregressive long-run variance from first differences -------------------
# AR coefficients and autocorrelations (lags 0..lm) from partial autocorrelations r
.pacf2acf <- function(r, lm) {
  p <- length(r); phi <- r[1]; L <- max(lm, p); rho <- numeric(L + 1); rho[1] <- 1; rho[2] <- r[1]
  if (p >= 2) for (k in 2:p) {
    pk <- phi
    rho[k + 1] <- r[k] * (1 - sum(pk * rho[2:k])) + sum(pk * rho[k:2])
    phi <- c(pk - r[k] * rev(pk), r[k])
  }
  if (L > p) for (l in (p + 1):L) rho[l + 1] <- sum(phi * rho[l:(l - p + 1)])
  list(phi = phi, rho = rho)
}
# autocovariances (lags 0..H) of the first difference of an AR process
.acov_diff1 <- function(r, s2u, H) {
  p <- length(r); a <- .pacf2acf(r, H + 2); phi <- a$phi; rho <- a$rho
  ge <- s2u / (1 - sum(phi * rho[2:(p + 1)])) * rho
  vapply(0:H, function(h) 2 * ge[h + 1] - ge[abs(h - 1) + 1] - ge[h + 2], numeric(1))
}
# Minimum-distance fit of AR(p) to the pooled autocovariances gh (lags 0..H) of the
# first differences; returns the long-run variance of the AR process
.lrv_ar_fit <- function(gh, p, Tn = NULL) {
  H <- length(gh) - 1L
  obj <- function(th) {
    v <- sum((gh - .acov_diff1(tanh(th[1:p]), exp(th[p + 1]), H))^2)
    if (is.finite(v)) v else 1e12
  }
  best <- NULL
  for (st in list(c(atanh(0.3), rep(0, p - 1)), c(atanh(0.8), rep(0, p - 1)))) {
    o <- try(stats::optim(c(st, log(gh[1] / 2)), obj, method = "BFGS", control = list(maxit = 200)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(c(NA_real_, NA_real_))
  r <- pmax(pmin(tanh(best$par[1:p]), 0.97), -0.97)
  phi <- .pacf2acf(r, 1L)$phi; s2u <- exp(best$par[p + 1])
  om <- s2u / (1 - sum(phi))^2
  ratio <- 1
  if (!is.null(Tn)) {
    a <- .pacf2acf(r, Tn); ge <- s2u / (1 - sum(phi * a$rho[2:(p + 1)])) * a$rho[1:Tn]
    v <- .step_vec(Tn)
    ratio <- drop(crossprod(v, stats::toeplitz(ge) %*% v)) / om
    if (!is.finite(ratio) || ratio <= 0) ratio <- 1
  }
  c(om, ratio)
}
# unit step at the middle, residualized on a constant, a trend and the hinge there: the
# increment direction that drives the maximum of the kink criterion
.step_vec <- function(Tn) {
  key <- paste("stepvec", Tn)
  if (!is.null(.pkg_cache[[key]])) return(.pkg_cache[[key]])
  tt <- seq_len(Tn); b <- round(Tn / 2)
  v <- qr.resid(qr(cbind(1, tt, pmax(tt - b, 0))), as.numeric(tt > b)); v <- v / sqrt(sum(v^2))
  .pkg_cache[[key]] <- v
  v
}
# Long-run variance of the series in Y from AR fits to the first differences (orders 2 to 4,
# geometric mean). With effective = TRUE the estimate is multiplied, for each fit, by the
# finite-sample ratio v'Gamma v / omega^2 of the variance of a unit step-like increment to the
# long-run variance (section 22): the variance that governs the maximum of the criterion and
# the knot's localization error at sample size T, which is below the long-run variance for
# positively autocorrelated errors.
.lrv_ar_diff <- function(Y, orders = 2:4, H = 10L, effective = TRUE) {
  D <- diff(Y); D <- sweep(D, 2L, colMeans(D), "-"); n <- nrow(D)
  if (n < 5 * H) return(NA_real_)
  gh <- vapply(0:H, function(h) mean(colSums(D[1:(n - h), , drop = FALSE] * D[(1 + h):n, , drop = FALSE]) / n), numeric(1))
  est <- vapply(orders, function(p) { f <- .lrv_ar_fit(gh, p, if (effective) nrow(Y) else NULL); f[1] * f[2] }, numeric(1))
  est <- est[is.finite(est) & est > 0]
  if (!length(est)) NA_real_ else exp(mean(log(est)))
}

# ---- cross-sectional dependence --------------------------------------------------
# Covariance of the demeaned first differences: shape of the cross-sectional covariance,
# barely affected by piecewise-linear means (contamination of rank at most the number of
# knots and norm at most R^2 / 4, section 19)
.cross_cov <- function(Y) {
  D <- diff(Y); D <- sweep(D, 2L, colMeans(D), "-")
  crossprod(D) / nrow(D)
}
.inv_sqrt <- function(C) {
  e <- eigen((C + t(C)) / 2, symmetric = TRUE)
  v <- pmax(e$values, 1e-8 * max(e$values))
  e$vectors %*% (t(e$vectors) / sqrt(v))
}

#' Monte Carlo penalty for the kink criterion with correlated series
#'
#' Like \code{\link{art_penalty_kink_mc}} but the series have covariance
#' eigenvalues \code{lambda} (errors independent over time): the statistic is
#' the largest over knots of the weighted sum of squared projections
#' \eqn{\sum_j \lambda_j (v_k' \xi_j)^2}. With \code{lambda = rep(1, N)} it
#' reduces to the unit-variance penalty.
#'
#' @param lambda Eigenvalues of the cross-sectional covariance (mean 1 for the
#'   penalty per unit long-run variance).
#' @param T,level,nsim,gap,seed As in \code{art_penalty_kink_mc}.
#' @return The penalty.
#' @export
art_penalty_kink_mc_w <- function(lambda, T, level = 0.05, nsim = 2000L, gap = 4L,
                                  seed = 20240101L) {
  key <- paste("kinkw", paste(round(lambda, 3), collapse = ","), T, level, nsim, gap, sep = "_")
  if (!is.null(.pkg_cache[[key]])) return(.pkg_cache[[key]])
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  N <- length(lambda); tt <- seq_len(T); ks <- gap:(T - gap)
  H <- outer(tt, ks, function(t, k) pmax(t - k, 0))
  q <- qr(cbind(1, tt)); Hres <- qr.resid(q, H); nrm <- colSums(Hres^2)
  mx <- numeric(nsim)
  for (i in seq_len(nsim)) {
    Yres <- qr.resid(q, matrix(rnorm(T * N), T, N))
    P <- crossprod(Hres, Yres)^2
    mx[i] <- max((P %*% lambda) / nrm)
  }
  if (had) assign(".Random.seed", old, envir = globalenv())
  else rm(".Random.seed", envir = globalenv())
  val <- unname(quantile(mx, 1 - level))
  .pkg_cache[[key]] <- val
  val
}


# Parzen long-run covariance matrix (N x N) of the residuals of the kink fit with knots K,
# with the degrees-of-freedom factor of .kink_lrv
.kink_lrcov <- function(Y, K, L, extra = 0L, gap = 6L) {
  Tn <- nrow(Y); N <- ncol(Y)
  E <- .kink_lrv(Y, K, L, extra, gap)$E
  Om <- crossprod(E) / Tn; wsum <- 1
  for (l in seq_len(L)) {
    w <- .kern(l / (L + 1), "parzen")
    G <- crossprod(E[1:(Tn - l), , drop = FALSE], E[(l + 1):Tn, , drop = FALSE]) / Tn
    Om <- Om + w * (G + t(G)); wsum <- wsum + 2 * w
  }
  p <- 2L + length(K) + extra
  Om <- Om / max(1 - p * wsum / Tn, 0.5)
  (Om + t(Om)) / 2
}


# Residual-based Parzen long-run variance of one score series (automatic bandwidth, df factor)
.score_lrv_resid <- function(sc, K, gap = 6L) {
  Tn <- length(sc); tt <- seq_len(Tn)
  B <- cbind(1, tt, if (length(K)) outer(tt, K, function(t, k) pmax(t - k, 0)) else NULL)
  e <- qr.resid(qr(B), sc)
  L <- .auto_bandwidth(matrix(e, ncol = 1), list(seq_len(Tn)), "parzen")
  v <- sum(e^2) / Tn; wsum <- 1
  for (l in seq_len(L)) { w <- .kern(l / (L + 1), "parzen"); v <- v + 2 * w * sum(e[1:(Tn - l)] * e[(l + 1):Tn]) / Tn; wsum <- wsum + 2 * w }
  max(v / max(1 - ncol(B) * wsum / Tn, 0.5), 1e-12)
}

# Long-run covariance from principal-direction-wise long-run variances: the leading r
# eigenvectors of the first-difference covariance each get their own variance (knot-free
# autoregressive estimate, optionally maxed with the residual estimate after the first
# pass); the remaining directions share one pooled variance and the shape of the
# first-difference covariance.
.pc_omega <- function(Y, K, L, ex, gap, first, Cs, pcinfo, pc_lrv, om2_scalar, pc_scale = 1, pc_floor = TRUE) {
  N <- ncol(Y); r <- pcinfo$r
  if (r == 0L) return(Cs * (N * om2_scalar / sum(diag(Cs))))
  Uf <- pcinfo$eg$vectors[, 1:r, drop = FALSE]
  S <- Y %*% Uf; Yp <- Y - S %*% t(Uf)
  omf <- vapply(seq_len(r), function(j) {
    a <- .lrv_ar_diff(S[, j, drop = FALSE])
    if (!is.finite(a)) a <- .score_lrv_resid(S[, j], K, gap)
    if (!first && pc_lrv == "max") a <- max(a, .score_lrv_resid(S[, j], K, gap))
    if (!first && pc_lrv == "resid") a <- .score_lrv_resid(S[, j], .augment_knots(Y, K, ex, gap), gap)
    a * pc_scale }, numeric(1))
  if (pc_floor) {
    sepf <- colSums(Uf * ((Cs * (N * om2_scalar / sum(diag(Cs)))) %*% Uf))
    omf <- pmax(omf, sepf)
  }
  ob <- if (first) .lrv_ar_diff(Yp) else .kink_lrv(Yp, K, L, ex, gap)$omega2
  if (!is.finite(ob)) ob <- .kink_lrv(Yp, K, L, ex, gap)$omega2
  ob <- ob * N / (N - r)
  Pp <- diag(N) - Uf %*% t(Uf)
  Cb <- Pp %*% Cs %*% Pp; Cb <- Cb * (N - r) / sum(diag(Cb))
  Om <- Uf %*% (omf * t(Uf)) + ob * Cb
  (Om + t(Om)) / 2
}


# Selected knots plus the `extra` best unpenalized knots, chosen greedily by gain on the whole panel
.augment_knots <- function(Y, K, extra = 0L, gap = 6L) {
  if (extra <= 0L) return(K)
  Tn <- nrow(Y); tt <- seq_len(Tn)
  Hall <- outer(tt, seq_len(Tn), function(t, k) pmax(t - k, 0))
  K2 <- K
  for (e in seq_len(extra)) {
    q <- qr(cbind(1, tt, if (length(K2)) outer(tt, K2, function(t, k) pmax(t - k, 0)) else NULL))
    ks <- gap:(Tn - gap); for (o in K2) ks <- ks[abs(ks - o) >= gap]
    if (!length(ks)) break
    Hres <- qr.resid(q, Hall[, ks, drop = FALSE]); nrm <- colSums(Hres^2)
    g <- rowSums(crossprod(Hres, qr.resid(q, Y))^2) / pmax(nrm, 1e-12); g[nrm < 1e-8] <- -Inf
    K2 <- c(K2, ks[which.max(g)])
  }
  K2
}
