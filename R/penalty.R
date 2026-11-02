#' Penalty (stopping threshold) for the pooled split criterion
#'
#' Under Gaussian errors the pooled criterion at a split point in a node without
#' breaks has mean \code{tr(Omega)}. The penalty adds fluctuation allowances.
#' With the true covariance (\code{hat = FALSE}) the penalty is
#' \code{tr + 4 * fro * sqrt(log T) + C0 * op * log T}; with an estimated
#' covariance (\code{hat = TRUE}) the constants 4 and C0 are replaced by 5 and
#' 1.1 * C0, as in the iterated scheme of the accompanying notes.
#'
#' The proofs require \code{C0 >= 100}; this is conservative. The default
#' \code{C0 = 10} is a practical choice that should be checked by simulation
#' (see \code{\link{art_sim_study}}).
#'
#' @param tr,fro,op Trace, Frobenius norm and operator norm of the (long-run)
#'   error covariance matrix.
#' @param T Number of dates.
#' @param C0 Constant multiplying the log T term.
#' @param hat Logical; use the constants for an estimated covariance.
#' @return A single number, the penalty on the scale of the sum of squares.
#' @export
art_penalty <- function(tr, fro, op, T, C0 = 10, hat = FALSE) {
  lt <- log(T)
  if (hat) tr + 5 * fro * sqrt(lt) + 1.1 * C0 * op * lt
  else tr + 4 * fro * sqrt(lt) + C0 * op * lt
}

#' Penalty for independent errors with common variance
#'
#' @param N Number of series.
#' @param T Number of dates.
#' @param sigma2 Error variance.
#' @param C0 See \code{\link{art_penalty}}.
#' @param df Parameters per series and segment: 1 for mean breaks, 2 for trend
#'   breaks (intercept and slope).
#' @return The penalty \code{N sigma2 + 4 sigma2 sqrt(N log T) + C0 sigma2 log T}.
#' @export
art_penalty_iid <- function(N, T, sigma2 = 1, C0 = 10, df = 1) {
  art_penalty(df * N * sigma2, sigma2 * sqrt(df * N), sigma2, T, C0, hat = FALSE)
}

#' Monte Carlo calibrated penalty
#'
#' Simulates independent N(0,1) panels of size T x N under the null of no
#' breaks and returns \code{sigma2} times the \code{1 - level} quantile of the
#' largest value of the pooled criterion over all split points (at least
#' \code{min_len} from the ends). This controls the probability of a false
#' split of a break-free segment at about \code{level}. It is much less
#' conservative than the union-bound penalty of \code{\link{art_penalty}}, which
#' is what the proofs need, and is the practical default for independent errors.
#' Results are cached per (N, T, level, nsim, min_len). The user's random
#' number stream is left untouched.
#'
#' @param N,T Panel dimensions.
#' @param sigma2 Error variance.
#' @param level Nominal null rejection probability for one break-free segment.
#' @param nsim Number of simulated null panels.
#' @param min_len Minimum distance of a split from the ends.
#' @param type \code{"mean"} or \code{"trend"} criterion.
#' @param seed Internal seed for the calibration simulation.
#' @return The penalty.
#' @export
art_penalty_mc <- function(N, T, sigma2 = 1, level = 0.05, nsim = 1000L,
                           min_len = 2L, type = c("mean", "trend"),
                           seed = 20240101L) {
  type <- match.arg(type)
  key <- paste(N, T, level, nsim, min_len, type, sep = "_")
  if (!is.null(.pkg_cache[[key]])) return(sigma2 * .pkg_cache[[key]])
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  mx <- numeric(nsim)
  for (i in seq_len(nsim)) {
    S0 <- .cumsum0(matrix(rnorm(T * N), T, N), type)
    mx[i] <- .best_split(S0, 0L, as.integer(T), as.integer(min_len),
                         as.integer(T - min_len))$V
  }
  if (had) assign(".Random.seed", old, envir = globalenv())
  else rm(".Random.seed", envir = globalenv())
  q <- unname(quantile(mx, 1 - level))
  .pkg_cache[[key]] <- q
  sigma2 * q
}

.pkg_cache <- new.env(parent = emptyenv())

#' Monte Carlo penalty for a given long-run covariance
#'
#' Like \code{\link{art_penalty_mc}}, but the null panels are drawn from
#' N(0, Omega) independently over time. For nodes that are long compared with
#' the serial dependence, the pooled criterion depends on the errors only
#' through the long-run covariance, so this approximates the null law under
#' serial and cross-sectional dependence. It is not cached, because Omega
#' changes between passes; a fixed internal seed is used so that the penalty
#' varies smoothly with Omega.
#'
#' @param Omega N x N long-run covariance matrix, or a vector of per-series
#'   long-run variances (diagonal covariance).
#' @param T Number of dates.
#' @param level Nominal null rejection probability for one break-free segment.
#' @param nsim Number of simulated panels.
#' @param min_len Minimum distance of a split from the ends.
#' @param type \code{"mean"} or \code{"trend"} criterion.
#' @param seed Internal seed.
#' @return The penalty.
#' @export
art_penalty_mc_cov <- function(Omega, T, level = 0.05, nsim = 100L,
                               min_len = 2L, type = c("mean", "trend"),
                               seed = 20240101L) {
  type <- match.arg(type)
  if (is.matrix(Omega)) {
    N <- nrow(Omega)
    ev <- eigen((Omega + t(Omega)) / 2, symmetric = TRUE)
    A <- ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0)))
  } else {
    N <- length(Omega); A <- NULL; sd_d <- sqrt(pmax(Omega, 0))
  }
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  mx <- numeric(nsim)
  for (i in seq_len(nsim)) {
    Z <- matrix(rnorm(T * N), T, N)
    Z <- if (!is.null(A)) Z %*% A else sweep(Z, 2L, sd_d, "*")
    mx[i] <- .best_split(.cumsum0(Z, type), 0L, as.integer(T), as.integer(min_len),
                         as.integer(T - min_len))$V
  }
  if (had) assign(".Random.seed", old, envir = globalenv())
  else rm(".Random.seed", envir = globalenv())
  unname(quantile(mx, 1 - level))
}

#' First-difference estimate of the error variance
#'
#' Valid for serially independent errors. Differencing removes the means except
#' at the break dates, so no knowledge of the breaks is needed. Under serial
#' correlation it estimates a short-run variance, not the long-run variance.
#'
#' @param Y Panel matrix (T x N).
#' @param order Order of differencing: 1 removes level shifts (mean breaks),
#'   2 also removes a linear trend (trend breaks).
#' @return Average over series of the variance estimate.
#' @export
art_sigma2_diff <- function(Y, order = 1L) {
  Y <- .as_panel(Y)
  D <- diff(Y, differences = order)
  sum(D^2) / (choose(2 * order, order) * ncol(Y) * nrow(D))
}

#' First-difference estimate of the cross-sectional error covariance
#' @param Y Panel matrix (T x N).
#' @param order Order of differencing, as in \code{\link{art_sigma2_diff}}.
#' @return N x N matrix.
#' @export
art_sigma_diff <- function(Y, order = 1L) {
  Y <- .as_panel(Y)
  D <- diff(Y, differences = order)
  crossprod(D) / (choose(2 * order, order) * nrow(D))
}

#' Symmetric inverse square root of a covariance matrix
#'
#' @param Omega Positive definite matrix.
#' @param tol Eigenvalues below \code{tol * max(eigenvalue)} are raised to that
#'   level.
#' @param bias_n If not \code{NULL}, the inverse is multiplied by
#'   \code{(bias_n - N - 1) / bias_n} (Wishart bias correction for a covariance
#'   estimated from \code{bias_n} observations) when \code{bias_n > N + 1}.
#' @return Matrix \code{Omega^(-1/2)}; use \code{Y \%*\% W} to whiten data.
#' @export
art_whiten_matrix <- function(Omega, tol = 1e-8, bias_n = NULL) {
  e <- eigen((Omega + t(Omega)) / 2, symmetric = TRUE)
  vals <- pmax(e$values, tol * max(e$values))
  scale <- 1
  N <- nrow(Omega)
  if (!is.null(bias_n) && bias_n > N + 1) scale <- (bias_n - N - 1) / bias_n
  e$vectors %*% (t(e$vectors) * (sqrt(scale) / sqrt(vals)))
}

#' Check the bandwidth window for long-run covariance estimation
#'
#' Reports the rough constraints (B1) and (B2) of the accompanying notes with
#' unit constants: the bandwidth should exceed about \code{sqrt(N / log T)} and
#' stay below \code{T / (N (mbar + 1) log T)}. The notes do not specify the
#' constants, so treat this as a diagnostic, not a rule.
#'
#' @param T,N Dimensions of the panel.
#' @param L Bandwidth.
#' @param mbar Upper bound on the number of breaks.
#' @return A list with the lower bound, upper bound, the ratio
#'   \code{N L log T / T}, and a logical \code{ok}.
#' @export
art_bandwidth_check <- function(T, N, L, mbar = 3) {
  lt <- log(T)
  lower <- sqrt(N / lt)
  upper <- T / (N * (mbar + 1) * lt)
  list(lower = lower, upper = upper, ratio = N * L * lt / T,
       ok = (L >= lower) && (L <= upper))
}
