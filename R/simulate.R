#' Simulate a panel with common mean breaks
#'
#' Generates \code{Y = mu + E}. The means \code{mu} are piecewise constant with
#' jumps \code{delta_j} (N-vectors) at the break dates. The per-series average
#' squared jump is \code{delta_bar^2}, so the pooled squared jump is
#' \code{N * delta_bar^2}.
#'
#' Error models (all with variance parameter \code{sigma}, the innovation
#' standard deviation):
#' \describe{
#'   \item{\code{"iid"}}{Independent over series and time.}
#'   \item{\code{"ar1"}}{Independent series, each AR(1) with coefficient \code{phi}.}
#'   \item{\code{"equicorr"}}{AR(1) with coefficient \code{phi} whose innovations
#'     have equicorrelation \code{rho} across series (a strong common factor).}
#'   \item{\code{"factor"}}{\code{e = Lambda f + u} with \code{r} factors, loadings
#'     with iid N(0,1) entries (so the cross product of the loadings divided by N is about the identity),
#'     factor variance \code{factor_var}, AR(1) dynamics with coefficient
#'     \code{phi} for factors and idiosyncratic terms.}
#' }
#'
#' Jump placement:
#' \describe{
#'   \item{\code{"random"}}{A random subset of \code{ceiling(density * N)} series
#'     jumps by \code{+- delta_bar}.}
#'   \item{\code{"factor"}}{The jump vector lies in the column space of the loadings.}
#'   \item{\code{"orthogonal"}}{The jump vector is orthogonal to the loadings.}
#' }
#' The last two use the loadings drawn for the error model (drawn anyway if the
#' errors are not \code{"factor"}).
#'
#' @param T,N Dimensions.
#' @param breaks Break dates (integers) or fractions of T in (0, 1).
#' @param delta_bar Per-series root-mean-square jump size. May be a vector with
#'   one value per break (heterogeneous jump sizes); it is recycled.
#' @param sigma Innovation standard deviation.
#' @param density Fraction of series that jump (for \code{jump_in = "random"}).
#' @param errors Error model.
#' @param phi AR(1) coefficient.
#' @param rho Cross-sectional innovation correlation for \code{"equicorr"}.
#' @param r Number of factors.
#' @param factor_var Factor variance.
#' @param jump_in Placement of the jumps.
#' @param same_direction If \code{TRUE}, every break has the same jump vector
#'   (a staircase). With one factor and \code{jump_in} equal to \code{"factor"},
#'   independent directions give offsetting pulses half the time, which are
#'   hard to detect by construction.
#' @param slope_bar Per-series root-mean-square change in slope at each break
#'   (vector recycled over breaks); the mean path then has a continuous kink
#'   at the break, in addition to any level jump from \code{delta_bar}. Use
#'   \code{delta_bar = 0} for a pure slope break.
#' @param seed Optional seed.
#' @return A list with \code{Y}, \code{mu}, \code{E}, \code{dates},
#'   \code{delta} (m x N), \code{Lambda} and the call arguments.
#' @examples
#' sim <- art_sim_data(T = 150, N = 8, breaks = c(0.5), delta_bar = 0.8, seed = 1)
#' sim$dates
#' @export
art_sim_data <- function(T, N, breaks, delta_bar = 1, sigma = 1, density = 1,
                         errors = c("iid", "ar1", "equicorr", "factor"),
                         phi = 0, rho = 0.5, r = 1L, factor_var = 1,
                         jump_in = c("random", "factor", "orthogonal"),
                         same_direction = FALSE, slope_bar = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  errors <- match.arg(errors); jump_in <- match.arg(jump_in)
  dates <- if (length(breaks) && all(breaks < 1)) as.integer(round(breaks * T)) else as.integer(breaks)
  m <- length(dates)
  Lam <- NULL
  if (errors == "factor" || jump_in != "random") Lam <- matrix(rnorm(N * r), N, r)
  delta <- matrix(0, m, N)
  db <- rep_len(delta_bar, max(m, 1L))
  for (j in seq_len(m)) {
    if (jump_in == "random") {
      idx <- sample.int(N, max(1L, ceiling(density * N)))
      delta[j, idx] <- db[j] * sample(c(-1, 1), length(idx), replace = TRUE)
    } else {
      d <- if (jump_in == "factor") as.numeric(Lam %*% rnorm(r)) else {
        z <- rnorm(N); z - as.numeric(Lam %*% solve(crossprod(Lam), crossprod(Lam, z)))
      }
      delta[j, ] <- d * db[j] * sqrt(N) / sqrt(sum(d^2))
    }
  }
  if (same_direction && m > 1L)
    for (j in 2:m) delta[j, ] <- delta[1L, ] * (if (db[1L] > 0) db[j] / db[1L] else 0)
  mu <- matrix(0, T, N)
  for (j in seq_len(m)) {
    rows <- (dates[j] + 1L):T
    mu[rows, ] <- mu[rows, ] + matrix(delta[j, ], length(rows), N, byrow = TRUE)
  }
  gamma <- matrix(0, m, N)
  if (any(slope_bar > 0) && m > 0L) {
    sb <- rep_len(slope_bar, m)
    for (j in seq_len(m)) {
      gamma[j, ] <- sb[j] * sample(c(-1, 1), N, replace = TRUE)
      rows <- (dates[j] + 1L):T
      mu[rows, ] <- mu[rows, ] + outer(seq_along(rows), gamma[j, ])
    }
  }
  burn <- 100L; Tt <- T + burn
  ar <- function(Z, phi) {
    if (phi == 0) return(Z)
    for (t in 2:nrow(Z)) Z[t, ] <- phi * Z[t - 1L, ] + Z[t, ]
    Z
  }
  E <- switch(errors,
    iid = matrix(rnorm(T * N, sd = sigma), T, N),
    ar1 = ar(matrix(rnorm(Tt * N, sd = sigma), Tt, N), phi)[-seq_len(burn), , drop = FALSE],
    equicorr = {
      Z <- sqrt(1 - rho) * matrix(rnorm(Tt * N, sd = sigma), Tt, N) +
        sqrt(rho) * matrix(rnorm(Tt, sd = sigma), Tt, N)
      ar(Z, phi)[-seq_len(burn), , drop = FALSE]
    },
    factor = {
      sdf <- sqrt(factor_var * (1 - phi^2))
      Fm <- ar(matrix(rnorm(Tt * r, sd = if (phi == 0) sqrt(factor_var) else sdf), Tt, r), phi)
      U <- ar(matrix(rnorm(Tt * N, sd = sigma), Tt, N), phi)
      (Fm %*% t(Lam) + U)[-seq_len(burn), , drop = FALSE]
    })
  list(Y = mu + E, mu = mu, E = E, dates = dates, delta = delta, gamma = gamma, Lambda = Lam,
       args = list(T = T, N = N, delta_bar = delta_bar, sigma = sigma,
                   errors = errors, phi = phi, jump_in = jump_in))
}

#' Score estimated break dates against the truth
#'
#' A true break is "found" if some estimated date lies within \code{tol} dates
#' of it. An estimated date is "spurious" if no true break lies within
#' \code{tol} dates of it.
#'
#' @param est,truth Integer vectors of dates.
#' @param T Number of dates.
#' @param tol Tolerance in dates; default \code{max(5, ceiling(0.02 * T))}.
#' @return A named numeric vector: \code{n_est}, \code{n_true},
#'   \code{exact_count}, \code{found}, \code{missed}, \code{spurious},
#'   \code{hausdorff} (as a fraction of T), \code{mae} (mean absolute error of
#'   the nearest estimate for each true break, in dates, \code{NA} if none).
#' @export
art_eval <- function(est, truth, T, tol = NULL) {
  if (is.null(tol)) tol <- max(5, ceiling(0.02 * T))
  est <- as.numeric(est); truth <- as.numeric(truth)
  m <- length(truth)
  dist_to <- function(a, b) if (!length(b)) rep(Inf, length(a)) else
    vapply(a, function(x) min(abs(x - b)), numeric(1))
  dt <- dist_to(truth, est); de <- dist_to(est, truth)
  haus <- if (!m && !length(est)) 0 else max(c(dt, de)[is.finite(c(dt, de))], 0)
  if ((m && !length(est)) || (!m && length(est))) haus <- 1 * T
  c(n_est = length(est), n_true = m, exact_count = as.numeric(length(est) == m),
    found = sum(dt <= tol), missed = sum(dt > tol), spurious = sum(de > tol),
    hausdorff = haus / T,
    mae = if (m && length(est)) mean(dt) else NA_real_)
}

#' Default set of methods for a simulation study
#'
#' Each method is a function taking the panel \code{Y} and returning a vector
#' of dates.
#' \itemize{
#'   \item \code{classic_bic}: classic ART on the pooled panel with BIC pruning.
#'   \item \code{panel_iid}: \code{art_fit} with the iid (union-bound) penalty.
#'   \item \code{panel_mc}: \code{art_fit} with the Monte Carlo penalty.
#'   \item \code{trend_mc}: as \code{panel_mc} but with trend breaks
#'     (\code{type = "trend"}).
#'   \item \code{classic_trend}: classic ART with a line fitted in each node.
#'   \item \code{kink_mc}: continuous join-point estimator \code{\link{art_kink}}.
#'   \item \code{panel_hac}: \code{art_fit} with the iterated long-run
#'     covariance and Monte Carlo penalty (\code{"hac_mc"}).
#'   \item \code{panel_hac_bound}: the same with the union-bound penalty
#'     (\code{"hac"}).
#'   \item \code{panel_whiten}: Monte Carlo penalty after whitening with the
#'     first-difference covariance.
#'   \item \code{panel_factor}: factor projection.
#' }
#' @return A named list of functions.
#' @export
art_methods <- function() {
  list(
    classic_bic = function(Y) art_classic(Y, "bic")$dates,
    panel_iid = function(Y) art_fit(Y, "iid")$dates,
    panel_mc = function(Y) art_fit(Y, "mc")$dates,
    trend_mc = function(Y) art_fit(Y, "mc", type = "trend")$dates,
    classic_trend = function(Y) art_classic(Y, "bic", type = "trend")$dates,
    kink_mc = function(Y) art_kink(Y, "mc")$dates,
    panel_hac = function(Y) art_fit(Y, "hac_mc")$dates,
    panel_hac_bound = function(Y) art_fit(Y, "hac")$dates,
    panel_whiten = function(Y) art_fit(Y, "mc", whiten = "diff")$dates,
    panel_factor = function(Y) art_fit_factor(Y)$dates
  )
}

#' Monte Carlo study over a grid of scenarios and methods
#'
#' \code{scenarios} is a data frame; each row is a scenario. Recognized columns
#' are the arguments of \code{\link{art_sim_data}} (\code{T}, \code{N},
#' \code{delta_bar}, \code{sigma}, \code{density}, \code{errors}, \code{phi},
#' \code{rho}, \code{r}, \code{factor_var}, \code{jump_in}) and \code{m}, the
#' number of breaks, placed at the fractions \code{j / (m + 1)}. Missing columns
#' take the defaults of \code{art_sim_data}; \code{T} and \code{N} are required.
#' Set \code{m = 0} for a null scenario (to measure false positives).
#'
#' @param scenarios Data frame of scenarios.
#' @param methods Named list of functions \code{function(Y)} returning dates.
#' @param nrep Replications per scenario.
#' @param seed Base seed; replication \code{i} of scenario \code{s} uses a
#'   seed derived from it, so every method sees the same data.
#' @param cores Number of cores for \code{parallel::mclapply}.
#' @return A data frame with one row per scenario, replication and method
#'   (columns of \code{\link{art_eval}} plus \code{time} and \code{failed}).
#' @examples
#' \donttest{
#' sc <- data.frame(T = 150, N = c(5, 10), m = 2, delta_bar = 0.5)
#' res <- art_sim_study(sc, art_methods()[c("classic_bic", "panel_iid")], nrep = 5)
#' art_sim_summary(res)
#' }
#' @export
art_sim_study <- function(scenarios, methods = art_methods(), nrep = 50L,
                          seed = 1L, cores = 1L) {
  stopifnot(is.data.frame(scenarios), all(c("T", "N") %in% names(scenarios)))
  jobs <- expand.grid(scenario = seq_len(nrow(scenarios)), rep = seq_len(nrep))
  run_one <- function(j) {
    s <- jobs$scenario[j]; rp <- jobs$rep[j]
    sc <- as.list(scenarios[s, , drop = FALSE])
    m <- if (!is.null(sc$m)) sc$m else 2
    sc$m <- NULL
    breaks <- if (m > 0) (seq_len(m)) / (m + 1) else numeric(0)
    set.seed(seed * 100003L + s * 1009L + rp)
    dat <- do.call(art_sim_data, c(list(breaks = breaks), sc))
    rows <- lapply(names(methods), function(nm) {
      t0 <- proc.time()[["elapsed"]]
      est <- tryCatch(methods[[nm]](dat$Y), error = function(e) NULL)
      el <- proc.time()[["elapsed"]] - t0
      if (is.null(est)) {
        data.frame(scenario = s, rep = rp, method = nm, t(setNames(rep(NA_real_, 8),
                   c("n_est", "n_true", "exact_count", "found", "missed",
                     "spurious", "hausdorff", "mae"))), time = el, failed = TRUE)
      } else {
        ev <- art_eval(est, dat$dates, sc$T)
        data.frame(scenario = s, rep = rp, method = nm, t(ev), time = el, failed = FALSE)
      }
    })
    do.call(rbind, rows)
  }
  res <- if (cores > 1L) mclapply(seq_len(nrow(jobs)), run_one, mc.cores = cores)
         else lapply(seq_len(nrow(jobs)), run_one)
  do.call(rbind, res)
}

#' Summarize a simulation study
#'
#' @param res Output of \code{\link{art_sim_study}}.
#' @return A data frame by scenario and method: mean number of estimated breaks,
#'   probability of the correct count, probability of at least one spurious
#'   break, mean fraction of true breaks found, mean Hausdorff distance (as a
#'   fraction of T), mean absolute error, mean time and the failure rate.
#' @export
art_sim_summary <- function(res) {
  key <- paste(res$scenario, res$method, sep = "|")
  parts <- split(res, key)
  out <- do.call(rbind, lapply(parts, function(d) {
    ok <- !d$failed
    data.frame(scenario = d$scenario[1], method = d$method[1],
               n_est = mean(d$n_est, na.rm = TRUE),
               n_true = d$n_true[ok][1],
               p_exact_count = mean(d$exact_count, na.rm = TRUE),
               p_any_spurious = mean(d$spurious > 0, na.rm = TRUE),
               frac_found = if (isTRUE(d$n_true[ok][1] > 0))
                 mean(d$found / d$n_true, na.rm = TRUE) else NA_real_,
               hausdorff = mean(d$hausdorff, na.rm = TRUE),
               mae = mean(d$mae, na.rm = TRUE),
               time = mean(d$time), fail_rate = mean(d$failed))
  }))
  rownames(out) <- NULL
  out[order(out$scenario, out$method), ]
}
