# Search steps: repartition, elimination, forward cycle, thinning, phase C ----

# Move each date to the maximizer of the criterion in the node spanned by its
# neighbours, keeping at least m dates from the node ends. Repeats until stable.
.repartition <- function(S0, dates, m, max_sweep = 50L) {
  Tn <- nrow(S0) - 1L
  if (!length(dates)) return(dates)
  dates <- sort(as.integer(dates))
  for (sw in seq_len(max_sweep)) {
    moved <- FALSE
    for (i in seq_along(dates)) {
      p <- if (i == 1L) 0L else dates[i - 1L]
      q <- if (i == length(dates)) Tn else dates[i + 1L]
      lo <- p + m; hi <- q - m
      if (hi < lo) next
      ks <- lo:hi
      V <- .node_V(S0, p, q, ks)
      cur <- dates[i]
      Vcur <- if (cur >= lo && cur <= hi) V[cur - lo + 1L] else -Inf
      j <- which.max(V)
      if (V[j] > Vcur + 1e-9) { dates[i] <- as.integer(ks[j]); moved <- TRUE }
    }
    if (!moved) break
  }
  dates
}

# Delete the date with the smallest removal cost while that cost is below
# (strict) or at most (non-strict) the penalty alpha; repartition after each
# deletion.
.eliminate <- function(S0, dates, alpha, m, strict = TRUE) {
  repeat {
    if (!length(dates)) break
    cst <- .costs(S0, dates)
    j <- which.min(cst)
    drop <- if (strict) cst[j] < alpha else cst[j] <= alpha
    if (!drop) break
    dates <- dates[-j]
    dates <- .repartition(S0, dates, m)
  }
  dates
}

#' Forward-repartition-eliminate cycle
#'
#' Starting from any set of dates, repeatedly (1) inserts a date at the best
#' split of the segment whose best split has the largest criterion value, if it
#' exceeds \code{alpha}; (2) moves every date to the maximizer of the criterion
#' in the node between its neighbours; (3) deletes the date with the smallest
#' removal cost if that cost is below \code{alpha}. Each effective step lowers
#' \code{SSE + alpha * (number of dates)}, so the cycle terminates.
#'
#' @param S0 Cumulative-sum matrix, with a leading zero row, of the (possibly
#'   whitened) panel.
#' @param alpha Penalty.
#' @param dates Starting dates.
#' @param min_len Minimum distance of a split from the ends of its node.
#' @param max_iter Safety cap on the number of cycles.
#' @return Sorted integer vector of dates.
#' @keywords internal
art_forward <- function(S0, alpha, dates = integer(0), min_len = 2L, max_iter = 500L) {
  Tn <- nrow(S0) - 1L
  dates <- sort(unique(as.integer(dates)))
  for (it in seq_len(max_iter)) {
    changed <- FALSE
    b <- c(0L, dates, Tn)
    bestV <- -Inf; bestk <- NA_integer_
    for (i in seq_len(length(b) - 1L)) {
      L <- b[i]; R <- b[i + 1L]
      bs <- .best_split(S0, L, R, L + min_len, R - min_len)
      if (bs$V > bestV) { bestV <- bs$V; bestk <- bs$k }
    }
    if (is.finite(bestV) && bestV > alpha) {
      dates <- sort(c(dates, bestk)); changed <- TRUE
    }
    d2 <- .repartition(S0, dates, min_len)
    if (!identical(d2, dates)) { dates <- d2; changed <- TRUE }
    d3 <- .eliminate(S0, dates, alpha, min_len)
    if (!identical(d3, dates)) { dates <- d3; changed <- TRUE }
    dates <- sort(unique(dates))
    if (!changed) break
  }
  as.integer(dates)
}

#' Hull-merge thinning
#'
#' Links dates closer than \code{3 * eps * T}. Every linked group with at least
#' two dates is replaced, in one parallel round, by the maximizer of the
#' criterion over the node between its two outside neighbours, restricted to
#' within \code{eps * T} of the group's span. A cleanup loop then merges any
#' remaining pair closer than \code{eps * T}.
#'
#' @param S0 Cumulative-sum matrix with a leading zero row.
#' @param dates Dates to thin.
#' @param eps Spacing fraction.
#' @return Sorted integer vector with spacing at least \code{ceiling(eps * T)}.
#' @keywords internal
art_thin <- function(S0, dates, eps = 0.02) {
  Tn <- nrow(S0) - 1L
  gap <- ceiling(3 * eps * Tn); pad <- ceiling(eps * Tn)
  dates <- sort(unique(as.integer(dates)))
  if (length(dates) < 2L) return(dates)
  comp <- cumsum(c(TRUE, diff(dates) >= gap))
  new <- integer(0)
  for (cc in unique(comp)) {
    idx <- which(comp == cc)
    if (length(idx) == 1L) { new <- c(new, dates[idx]); next }
    p_out <- if (min(idx) == 1L) 0L else dates[min(idx) - 1L]
    q_out <- if (max(idx) == length(dates)) Tn else dates[max(idx) + 1L]
    lo <- max(p_out + 1L, dates[min(idx)] - pad)
    hi <- min(q_out - 1L, dates[max(idx)] + pad)
    ks <- lo:hi
    V <- .node_V(S0, p_out, q_out, ks)
    new <- c(new, ks[which.max(V)])
  }
  dates <- sort(unique(as.integer(new)))
  # cleanup: enforce spacing by merging the closest pair
  while (length(dates) >= 2L && min(diff(dates)) < pad) {
    i <- which.min(diff(dates))
    p_out <- if (i == 1L) 0L else dates[i - 1L]
    q_out <- if (i + 1L == length(dates)) Tn else dates[i + 2L]
    lo <- max(p_out + 1L, dates[i] - pad)
    hi <- min(q_out - 1L, dates[i + 1L] + pad)
    ks <- lo:hi
    V <- .node_V(S0, p_out, q_out, ks)
    dates <- sort(unique(c(dates[-c(i, i + 1L)], ks[which.max(V)])))
  }
  as.integer(dates)
}

# Phase C: repartition and eliminate with minimum spacing eps * T
.phaseC <- function(S0, dates, alpha, eps, min_len = 2L) {
  Tn <- nrow(S0) - 1L
  m <- max(as.integer(min_len), as.integer(ceiling(eps * Tn)))
  dates <- .repartition(S0, dates, m)
  .eliminate(S0, dates, alpha, m, strict = FALSE)
}
