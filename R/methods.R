#' @export
print.art <- function(x, ...) {
  cat("panelart fit (", x$penalty, " penalty", if (x$whitened) ", whitened", if (identical(x$type, "trend")) ", trend breaks" else if (identical(x$type, "kink")) ", kinks", ")\n", sep = "")
  cat("  T =", nrow(x$Y), " N =", ncol(x$Y), "\n")
  cat("  penalty alpha =", format(x$alpha, digits = 5), "\n")
  cat("  breaks:", if (x$n_breaks) paste(x$dates, collapse = ", ") else "none", "\n")
  invisible(x)
}

#' @export
summary.art <- function(object, ...) {
  b <- object$segments
  b$length <- b$end - b$start + 1L
  b$mean_over_series <- rowMeans(object$means)
  out <- list(call = object$call, segments = b, alpha = object$alpha,
              history = object$history, ci = art_ci(object))
  class(out) <- "summary.art"
  out
}

#' @export
print.summary.art <- function(x, ...) {
  print(x$segments, row.names = FALSE)
  cat("\nPenalty:", format(x$alpha, digits = 5), "\n")
  if (nrow(x$ci)) { cat("\n95% intervals from the limit law:\n"); print(x$ci, row.names = FALSE) }
  if (!is.null(x$history)) { cat("\nCovariance passes:\n"); print(x$history, row.names = FALSE) }
  invisible(x)
}

#' Plot a fit
#'
#' Draws the series (thin grey lines), the cross-sectional mean (black) and
#' the segment means of that mean (red), with dashed vertical lines at the
#' estimated break dates.
#'
#' @param x An object from \code{\link{art_fit}}.
#' @param max_series Maximum number of individual series drawn.
#' @param ... Unused.
#' @export
plot.art <- function(x, max_series = 20L, ...) {
  Y <- x$Y; Tn <- nrow(Y)
  k <- seq_len(min(ncol(Y), max_series))
  matplot(Y[, k, drop = FALSE], type = "l", lty = 1,
          col = adjustcolor("grey50", 0.5), xlab = "time", ylab = "value")
  lines(rowMeans(Y), lwd = 2)
  b <- c(0L, x$dates, Tn)
  for (i in seq_len(length(b) - 1L)) {
    idx <- (b[i] + 1L):b[i + 1L]
    lines(idx, rep(mean(rowMeans(Y)[idx]), length(idx)), col = "red", lwd = 2)
  }
  abline(v = x$dates + 0.5, lty = 2)
  invisible(x)
}

#' @export
print.art_factor <- function(x, ...) {
  cat("panelart factor-projection fit: r =", x$r, ", rbar =", x$rbar, "\n")
  cat("  factor-pipeline breaks:",
      if (length(x$factor_fit$dates)) paste(x$factor_fit$dates, collapse = ", ") else "none", "\n")
  cat("  orthogonal-pipeline breaks:",
      if (length(x$orth_fit$dates)) paste(x$orth_fit$dates, collapse = ", ") else "none", "\n")
  cat("  merged breaks:", if (length(x$dates)) paste(x$dates, collapse = ", ") else "none", "\n")
  invisible(x)
}
