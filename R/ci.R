# Distribution of xi = argmax_v { W(v) - |v|/2 } -------------------------------

#' Density of the limit of the scaled break-date error
#'
#' Density of \eqn{\xi = \arg\max_v \{W(v) - |v|/2\}} with \eqn{W} a two-sided
#' standard Brownian motion, in the closed form
#' \eqn{g(x) = \tfrac32 e^{|x|}\Phi(-\tfrac32\sqrt{|x|}) -
#' \tfrac12\Phi(-\tfrac12\sqrt{|x|})} (Yao, 1987; Bai, 1997). The formula
#' was checked by simulation in the package tests.
#'
#' @param x Numeric vector.
#' @return Density values.
#' @export
art_xi_density <- function(x) {
  ax <- abs(x)
  1.5 * exp(ax + pnorm(-1.5 * sqrt(ax), log.p = TRUE)) - 0.5 * pnorm(-0.5 * sqrt(ax))
}

#' Quantile of |xi|
#' @param level Probability in (0, 1).
#' @return The value q with P(|xi| <= q) = level.
#' @export
art_xi_quantile <- function(level = 0.95) {
  cdf <- function(q) 2 * integrate(art_xi_density, 0, q)$value
  uniroot(function(q) cdf(q) - level, c(1e-6, 500))$root
}

#' Confidence intervals for break dates from the limit law
#'
#' For an isolated break, \eqn{\Delta (\hat k - b) / \omega^2 \Rightarrow \xi}
#' where \eqn{\Delta = \|\delta\|^2} is the squared jump and
#' \eqn{\omega^2 = \delta' \Omega \delta / \|\delta\|^2} the long-run variance
#' in the direction of the jump. The interval is
#' \code{date +- q * omega2 / Delta}, with the jump estimated by the difference
#' of adjacent segment means. For whitened fits the covariance is the identity.
#' Without a covariance estimate the first-difference variance estimate is used.
#'
#' The intervals are asymptotic for a shrinking jump and a well separated
#' break. The raw estimate of \code{Delta} is biased upward by the noise in the
#' two segment means, by about \code{(1/n_left + 1/n_right) * tr(Omega)}. With
#' \code{bias_correct = TRUE} (default) this is subtracted, with a floor at a
#' quarter of the raw value.
#'
#' @param fit An object from \code{\link{art_fit}} (type \code{"mean"}) or
#'   \code{\link{art_kink}}. Not available for \code{type = "trend"}.
#' @param level Confidence level.
#' @param bias_correct Subtract the noise bias from the squared jump (mean fits)
#'   or from the squared slope change (kink fits).
#' @return A data frame with the dates, lower and upper limits, the estimated
#'   squared jump and the long-run variance used.
#' @export
art_ci <- function(fit, level = 0.95, bias_correct = TRUE) {
  stopifnot(inherits(fit, "art"))
  if (identical(fit$type, "kink")) return(.ci_kink(fit, level, bias_correct))
  if (!identical(fit$type, "mean"))
    stop("art_ci is not available for type = 'trend': the limit law differs")
  d <- fit$dates
  if (!length(d)) return(data.frame(date = integer(0), lower = numeric(0),
                                    upper = numeric(0), Delta = numeric(0),
                                    omega2 = numeric(0)))
  Yw <- fit$Yw; Tn <- nrow(Yw)
  q <- art_xi_quantile(level)
  b <- c(0L, d, Tn)
  out <- data.frame(date = d, lower = NA_real_, upper = NA_real_,
                    Delta = NA_real_, omega2 = NA_real_)
  for (i in seq_along(d)) {
    mL <- colMeans(Yw[(b[i] + 1L):b[i + 1L], , drop = FALSE])
    mR <- colMeans(Yw[(b[i + 1L] + 1L):b[i + 2L], , drop = FALSE])
    delta <- mR - mL
    D <- sum(delta^2)
    if (bias_correct) {
      trv <- if (fit$whitened) ncol(Yw)
             else if (!is.null(fit$omega)) fit$omega$tr
             else ncol(Yw) * (if (!is.null(fit$sigma2)) fit$sigma2 else art_sigma2_diff(Yw))
      nL <- b[i + 1L] - b[i]; nR <- b[i + 2L] - b[i + 1L]
      D <- max(D - (1 / nL + 1 / nR) * trv, 0.25 * D)
    }
    om2 <- if (fit$whitened) 1
           else if (!is.null(fit$omega) && !is.null(fit$omega$Omega))
             sum(delta * (fit$omega$Omega %*% delta)) / D
           else if (!is.null(fit$sigma2)) fit$sigma2
           else art_sigma2_diff(Yw)
    hw <- q * om2 / D
    out[i, c("lower", "upper", "Delta", "omega2")] <- c(d[i] - hw, d[i] + hw, D, om2)
  }
  out
}

# Intervals for kinks. For a knot k_j the Gauss-Newton information is
#   ||gamma_j||^2 * ||(I - P_B) step_j||^2 / sigma^2,
# where B is the full design (1, t, all hinges) and step_j = 1(t > k_j) is minus
# the derivative of the hinge with respect to the knot. The date error is
# approximately normal with variance 1 / information, times the factor
#   1 + N sigma^2 / (||gamma_j||^2 ||h_j||^2)
# from the noise in the quadratic part of the criterion, where h_j is the
# hinge of knot j residualized on the other regressors (inflate = TRUE).
.ci_kink <- function(fit, level, bias_correct, inflate = TRUE) {
  K <- fit$dates; Y <- fit$Y; Tn <- nrow(Y); N <- ncol(Y)
  if (!length(K)) return(data.frame(date = integer(0), lower = numeric(0),
                                    upper = numeric(0), gamma2 = numeric(0), sd = numeric(0)))
  tt <- seq_len(Tn); Kn <- length(K)
  B <- cbind(1, tt, outer(tt, K, function(t, k) pmax(t - k, 0)))
  qB <- qr(B)
  s2 <- if (!is.null(fit$Omega)) sum(diag(fit$Omega)) / N else if (!is.null(fit$omega2)) fit$omega2 else sum(qr.resid(qB, Y)^2) / (N * (Tn - ncol(B)))
  z <- stats::qnorm(1 - (1 - level) / 2)
  R <- matrix(vapply(K, function(k) qr.resid(qB, as.numeric(tt > k)), numeric(Tn)), Tn, Kn)
  G <- matrix(fit$gamma, Kn, N)
  M <- (G %*% t(G)) * crossprod(R)           # joint information, up to 1 / sigma^2
  g2 <- diag(G %*% t(G)); infl <- rep(1, Kn)
  for (j in seq_len(Kn)) {
    hj <- qr.resid(qr(B[, -(2 + j), drop = FALSE]), B[, 2 + j])
    nh <- sum(hj^2)
    if (bias_correct) { g2[j] <- max(g2[j] - N * s2 / nh, 0.25 * g2[j]); M[j, j] <- g2[j] * sum(R[, j]^2) }
    if (inflate) infl[j] <- 1 + N * s2 / (g2[j] * nh)
  }
  V <- tryCatch(solve(M), error = function(e) diag(1 / diag(M), Kn))
  sdj <- sqrt(s2 * pmax(diag(V), 0) * infl)
  if (!is.null(fit$Omega)) {                      # sandwich with the long-run covariance of the errors
    Om <- fit$Omega
    Sc <- (G %*% Om %*% t(G)) * crossprod(R)
    Vs <- V %*% Sc %*% V
    sdj <- sqrt(pmax(diag(Vs), 0) * infl)
  }
  data.frame(date = K, lower = K - z * sdj, upper = K + z * sdj, gamma2 = g2, sd = sdj)
}
