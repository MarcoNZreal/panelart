# Internal helpers ---------------------------------------------------------

.as_panel <- function(Y) {
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1L)
  if (!is.numeric(Y)) stop("Y must be numeric")
  if (anyNA(Y)) stop("Y contains missing values")
  storage.mode(Y) <- "double"
  Y
}

# (T+1) x N matrix of cumulative sums with a leading row of zeros. For
# type = "trend" the attribute "S1" holds cumulative sums of t * y.
.cumsum0 <- function(Y, type = "mean") {
  S <- Y
  for (j in seq_len(ncol(Y))) S[, j] <- cumsum(Y[, j])
  S <- rbind(0, S)
  attr(S, "type") <- type
  if (type == "trend") {
    S1 <- Y * seq_len(nrow(Y))
    for (j in seq_len(ncol(Y))) S1[, j] <- cumsum(S1[, j])
    attr(S, "S1") <- rbind(0, S1)
  }
  S
}

# Explained sum of squares (pooled over series) of a linear fit on the segments
# (a, b]; a and b are vectors of equal length (or scalars).
.expl_trend <- function(S0, S1, a, b) {
  n <- b - a
  Sy <- S0[b + 1L, , drop = FALSE] - S0[a + 1L, , drop = FALSE]
  Sty <- S1[b + 1L, , drop = FALSE] - S1[a + 1L, , drop = FALSE]
  tbar <- (a + 1 + b) / 2
  Suy <- Sty - tbar * Sy
  Suu <- n * (n^2 - 1) / 12
  rowSums(Sy^2) / n + rowSums(Suy^2) / Suu
}

# Reduction in the sum of squares from splitting (L, R] at k when each side
# gets its own line (intercept and slope) in every series.
.node_V_trend <- function(S0, L, R, ks) {
  S1 <- attr(S0, "S1")
  left <- .expl_trend(S0, S1, rep(L, length(ks)), ks)
  right <- .expl_trend(S0, S1, ks, rep(R, length(ks)))
  whole <- .expl_trend(S0, S1, L, R)
  left + right - whole
}

# Pooled split criterion V_I(k) for split points ks inside the node (L, R].
# S0 is the cumulative-sum matrix from .cumsum0. Implements
#   V = ||W||^2 / nu,  W = partial sum - (k-L)/n * node sum,
#   nu = (k-L)(R-k)/n
# which equals the reduction in the total sum of squares from splitting at k.
.node_V <- function(S0, L, R, ks) {
  if (identical(attr(S0, "type"), "trend")) return(.node_V_trend(S0, L, R, ks))
  n <- R - L
  B <- S0[R + 1L, ] - S0[L + 1L, ]
  A <- S0[ks + 1L, , drop = FALSE]
  A <- sweep(A, 2L, S0[L + 1L, ], "-")
  W <- A - outer((ks - L) / n, B)
  nu <- (ks - L) * (R - ks) / n
  rowSums(W^2) / nu
}

.best_split <- function(S0, L, R, lo, hi) {
  if (hi < lo) return(list(k = NA_integer_, V = -Inf))
  ks <- lo:hi
  V <- .node_V(S0, L, R, ks)
  j <- which.max(V)
  list(k = as.integer(ks[j]), V = V[j])
}

# Removal cost of each date: criterion value at the date inside the node
# spanned by its two neighbours.
.costs <- function(S0, dates) {
  Tn <- nrow(S0) - 1L
  b <- c(0L, dates, Tn)
  vapply(seq_along(dates),
         function(i) .node_V(S0, b[i], b[i + 2L], dates[i]), numeric(1))
}

.segment_slopes <- function(Y, dates) {
  b <- c(0L, dates, nrow(Y))
  M <- matrix(0, length(b) - 1L, ncol(Y))
  for (i in seq_len(length(b) - 1L)) {
    idx <- (b[i] + 1L):b[i + 1L]
    u <- idx - mean(idx)
    M[i, ] <- as.numeric(crossprod(u, Y[idx, , drop = FALSE])) / sum(u^2)
  }
  M
}

.segment_means <- function(Y, dates) {
  b <- c(0L, dates, nrow(Y))
  M <- matrix(0, length(b) - 1L, ncol(Y))
  for (i in seq_len(length(b) - 1L))
    M[i, ] <- colMeans(Y[(b[i] + 1L):b[i + 1L], , drop = FALSE])
  M
}
