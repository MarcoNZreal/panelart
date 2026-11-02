#' Classic ART: greedy regression tree on time with pruning
#'
#' Replicates the procedure of Rea, Reale, Cappelli and Brown (2010) in spirit,
#' for a single series or for a pooled panel: a least-squares regression tree
#' with time as the only predictor. Nodes are split at the point of largest
#' reduction in the (pooled) sum of squares until a node is smaller than
#' \code{minsize}, would create a child smaller than \code{mincut}, or has
#' deviance below \code{mindev} times the root deviance (the defaults mimic the
#' R package \pkg{tree}). The full tree is then pruned along the weakest-link
#' sequence.
#'
#' \code{pruning = "bic"} picks the subtree minimizing a BIC-type criterion
#' \code{N T log(SSE/(N T)) + (N (q + 1) + q) log(N T)} with \code{q} splits. This
#' is an approximation of the BIC rule used with ART, not a copy of any specific
#' implementation. \code{pruning = "alpha"} picks the subtree minimizing
#' \code{SSE + alpha * q}.
#'
#' @param Y Panel matrix (T x N) or a vector.
#' @param pruning One of \code{"bic"}, \code{"alpha"}, \code{"none"}.
#' @param alpha Penalty for \code{pruning = "alpha"}.
#' @param minsize Minimum node size to attempt a split.
#' @param mincut Minimum child size.
#' @param mindev Minimum node deviance as a fraction of the root deviance.
#' @param type \code{"mean"}: piecewise-constant fits; \code{"trend"}: a line is
#'   fitted in each node (two parameters per series and segment).
#' @return A list with \code{dates}, \code{sequence} (data frame with \code{q},
#'   \code{sse}, \code{bic}) and \code{tree} (node table).
#' @examples
#' set.seed(2)
#' y <- c(rnorm(100), rnorm(100, mean = 1.5))
#' art_classic(y)$dates
#' @export
art_classic <- function(Y, pruning = c("bic", "alpha", "none"), alpha = NULL,
                        minsize = 10L, mincut = 5L, mindev = 0.01,
                        type = c("mean", "trend")) {
  pruning <- match.arg(pruning); type <- match.arg(type)
  df <- if (type == "trend") 2 else 1
  if (type == "trend") mincut <- max(mincut, 3L)
  Y <- .as_panel(Y); Tn <- nrow(Y); N <- ncol(Y)
  S0 <- .cumsum0(Y, type); S2 <- .cumsum0(Y^2)
  dev <- function(L, R) {
    tot <- sum(S2[R + 1L, ] - S2[L + 1L, ])
    if (type == "trend") tot - .expl_trend(S0, attr(S0, "S1"), L, R)
    else tot - sum((S0[R + 1L, ] - S0[L + 1L, ])^2) / (R - L)
  }
  root_dev <- dev(0L, Tn)
  nodes <- data.frame(L = 0L, R = Tn, k = NA_integer_, gain = NA_real_,
                      left = NA_integer_, right = NA_integer_)
  i <- 1L
  while (i <= nrow(nodes)) {
    L <- nodes$L[i]; R <- nodes$R[i]
    if ((R - L) >= minsize && dev(L, R) >= mindev * root_dev) {
      bs <- .best_split(S0, L, R, L + mincut, R - mincut)
      if (is.finite(bs$V) && bs$V > 1e-12) {
        nodes$k[i] <- bs$k; nodes$gain[i] <- bs$V
        nodes <- rbind(nodes,
                       data.frame(L = L, R = bs$k, k = NA_integer_, gain = NA_real_,
                                  left = NA_integer_, right = NA_integer_),
                       data.frame(L = bs$k, R = R, k = NA_integer_, gain = NA_real_,
                                  left = NA_integer_, right = NA_integer_))
        nodes$left[i] <- nrow(nodes) - 1L; nodes$right[i] <- nrow(nodes)
      }
    }
    i <- i + 1L
  }
  active <- !is.na(nodes$k)
  branch <- function(h) {
    if (!active[h]) return(c(0, 0))
    l <- branch(nodes$left[h]); r <- branch(nodes$right[h])
    c(nodes$gain[h] + l[1] + r[1], 1 + l[2] + r[2])
  }
  reach <- function(h) {
    if (!active[h]) return(integer(0))
    c(h, reach(nodes$left[h]), reach(nodes$right[h]))
  }
  seqs <- list()
  repeat {
    r <- reach(1L)
    br <- branch(1L)
    seqs[[length(seqs) + 1L]] <- list(dates = sort(nodes$k[r]), q = length(r),
                                      sse = root_dev - br[1])
    if (!length(r)) break
    g <- vapply(r, function(h) { b <- branch(h); b[1] / b[2] }, numeric(1))
    active[r[which.min(g)]] <- FALSE
  }
  q <- vapply(seqs, function(s) s$q, numeric(1))
  sse <- pmax(vapply(seqs, function(s) s$sse, numeric(1)), 1e-12)
  bic <- N * Tn * log(sse / (N * Tn)) + (df * N * (q + 1) + q) * log(N * Tn)
  pick <- switch(pruning,
                 none = 1L,
                 bic = which.min(bic),
                 alpha = {
                   if (is.null(alpha)) stop("alpha required for pruning = 'alpha'")
                   which.min(sse + alpha * q)
                 })
  list(dates = as.integer(seqs[[pick]]$dates),
       sequence = data.frame(q = q, sse = sse, bic = bic), tree = nodes)
}
