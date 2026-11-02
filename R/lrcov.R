#' Long-run covariance estimate from segment residuals
#'
#' Kernel estimate computed from residuals after removing the mean of
#' each segment of the partition given by \code{dates}. With \code{dates =
#' integer(0)} the estimate is the raw one that leaves the breaks in the data.
#' Only within-segment pairs of residuals are used.
#'
#' For large panels (\code{N > full_max}) the full matrix is not formed; the
#' function then uses per-series long-run variances (a diagonal approximation,
#' valid under cross-sectional independence) and reports the trace, the
#' Frobenius norm and the operator norm of that diagonal matrix.
#'
#' @param Y Panel matrix (T x N).
#' @param dates Break dates defining the segments.
#' @param L Bandwidth (number of lags). Default \code{floor(T^(1/3))}. The
#'   value \code{"auto"} uses Andrews' (1991) AR(1) plug-in rule for the Bartlett
#'   kernel, pooled over series and computed from the segment residuals.
#' @param kernel Lag window: \code{"parzen"} (default; bias of order 1/L^2) or
#'   \code{"bartlett"} (bias of order 1/L, as in the theory notes). Both give
#'   positive semidefinite estimates.
#' @param bias_correct Logical; divide by \code{1 - (number of segments) *
#'   (sum of kernel weights) / T}, which removes the first-order downward bias
#'   from demeaning each segment.
#' @param type \code{"mean"}: residuals are deviations from segment means;
#'   \code{"trend"}: residuals from a line fitted in each segment.
#' @param full Logical; form the full N x N matrix. Default: \code{N <= full_max}.
#' @param full_max Largest N for which the full matrix is formed by default.
#' @return A list with \code{Omega} (matrix or \code{NULL}), \code{tr},
#'   \code{fro}, \code{op}, \code{L} and \code{full}.
#' @export
art_lrcov <- function(Y, dates = integer(0), L = NULL, full = NULL, full_max = 200L,
                      kernel = c("parzen", "bartlett"), bias_correct = TRUE,
                      type = c("mean", "trend")) {
  kernel <- match.arg(kernel); type <- match.arg(type)
  Y <- .as_panel(Y)
  Tn <- nrow(Y); N <- ncol(Y)
  b <- c(0L, sort(as.integer(dates)), Tn)
  E <- Y
  segs <- vector("list", length(b) - 1L)
  for (i in seq_len(length(b) - 1L)) {
    idx <- (b[i] + 1L):b[i + 1L]
    segs[[i]] <- idx
    E[idx, ] <- .resid_seg(Y[idx, , drop = FALSE], type)
  }
  if (is.null(L)) L <- max(1L, floor(Tn^(1 / 3)))
  else if (identical(L, "auto")) L <- .auto_bandwidth(E, segs, kernel)
  L <- as.integer(L)
  fac <- 1
  if (bias_correct) {
    wsum <- 1 + 2 * sum(vapply(seq_len(L), function(l) .kern(l / (L + 1), kernel), numeric(1)))
    fac <- 1 / max(1 - (if (type == "trend") 2 else 1) * length(segs) * wsum / Tn, 0.5)
  }
  if (is.null(full)) full <- N <= full_max
  if (full) {
    Om <- crossprod(E) / Tn
    for (l in seq_len(L)) {
      w <- .kern(l / (L + 1), kernel)
      G <- matrix(0, N, N)
      for (idx in segs) {
        n <- length(idx)
        if (n > l) G <- G + crossprod(E[idx[1:(n - l)], , drop = FALSE],
                                      E[idx[(l + 1):n], , drop = FALSE])
      }
      G <- G / Tn
      Om <- Om + w * (G + t(G))
    }
    Om <- Om * fac
    ev <- eigen((Om + t(Om)) / 2, symmetric = TRUE, only.values = TRUE)$values
    list(Omega = Om, tr = sum(diag(Om)), fro = sqrt(sum(Om^2)),
         op = max(ev), L = L, full = TRUE)
  } else {
    d <- colSums(E^2) / Tn
    for (l in seq_len(L)) {
      w <- .kern(l / (L + 1), kernel)
      g <- numeric(N)
      for (idx in segs) {
        n <- length(idx)
        if (n > l) g <- g + colSums(E[idx[1:(n - l)], , drop = FALSE] *
                                    E[idx[(l + 1):n], , drop = FALSE])
      }
      d <- d + 2 * w * g / Tn
    }
    d <- pmax(d * fac, 0)
    list(Omega = NULL, diag = d, tr = sum(d), fro = sqrt(sum(d^2)), op = max(d),
         L = L, full = FALSE)
  }
}

.kern <- function(x, kernel) {
  if (kernel == "bartlett") return(1 - x)
  if (x <= 0.5) 1 - 6 * x^2 + 6 * x^3 else 2 * (1 - x)^3
}

# Andrews (1991) AR(1) plug-in bandwidth, pooled over series, from
# segment-demeaned residuals E. Segments are index vectors.
.auto_bandwidth <- function(E, segs, kernel = "parzen") {
  Tn <- nrow(E); N <- ncol(E)
  num <- 0; den <- 0
  for (j in seq_len(N)) {
    e <- E[, j]
    v <- sum(e^2) / Tn
    if (v <= 0) next
    c1 <- 0
    for (idx in segs) if (length(idx) > 1L)
      c1 <- c1 + sum(e[idx[-1L]] * e[idx[-length(idx)]])
    rho <- max(min(c1 / (v * Tn), 0.95), -0.95)
    s2 <- v * (1 - rho^2)
    if (kernel == "bartlett")
      num <- num + 4 * rho^2 * s2^2 / ((1 - rho)^6 * (1 + rho)^2)
    else
      num <- num + 4 * rho^2 * s2^2 / (1 - rho)^8
    den <- den + s2^2 / (1 - rho)^4
  }
  a <- if (den > 0) num / den else 0
  L <- if (kernel == "bartlett") ceiling(1.1447 * (a * Tn)^(1 / 3))
       else ceiling(2.6614 * (a * Tn)^(1 / 5))
  as.integer(max(1L, min(L, floor(Tn / 5))))
}

# Residuals of the columns of Ysub after removing a constant (or a line)
.resid_seg <- function(Ysub, type = "mean") {
  if (type == "mean" || nrow(Ysub) < 3L)
    return(sweep(Ysub, 2L, colMeans(Ysub), "-"))
  u <- seq_len(nrow(Ysub)) - (nrow(Ysub) + 1) / 2
  Yc <- sweep(Ysub, 2L, colMeans(Ysub), "-")
  Yc - outer(u, as.numeric(crossprod(u, Yc)) / sum(u^2))
}
