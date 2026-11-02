library(panelart)
Tn <- 300; tt <- 1:Tn; b <- 100
B <- function(K) cbind(1, tt, if (length(K)) sapply(K, function(k) pmax(tt-k,0)) else NULL)
hb <- pmax(tt-b,0)
d2 <- function(K) sum(qr.resid(qr(B(K)), hb)^2)   # misfit per unit ||gamma||^2 (one true knot, others absent)
cat("(1) misfit per unit ||gamma||^2 of approximating a kink at b by knots placed elsewhere\n")
for (pr in list(c(109,124), c(105,115), c(110,130), c(88,94), c(80,88))) {
  cat(sprintf("  pair at b%+d, b%+d: %.0f   | best single knot in [b-30,b+30] other than b: %.0f at b+1 | single at the nearer pair knot: %.0f\n",
      pr[1]-b, pr[2]-b, d2(pr), min(sapply(c((b-30):(b-1),(b+1):(b+30)), function(k) d2(k))), d2(pr[which.min(abs(pr-b))])))
}
cat("  with ||gamma||^2 = 10 * 0.04^2 = 0.016, multiply by 0.016 to get the misfit in sigma^2 units\n")
## (2) correct re-optimized costs for duplicates that survive the search (Monte Carlo penalty, thin = FALSE)
N <- 10; s <- 0.04; rows <- list()
sse <- function(K, Y) sum(qr.resid(qr(B(K)), Y)^2)
for (i in 1:200) {
  set.seed(23000 + i)
  g1 <- s * sample(c(-1,1), N, TRUE)
  Y <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g1) + matrix(rnorm(Tn*N), Tn, N)
  f <- art_kink(Y, "mc", nsim=500, thin=FALSE)
  K <- f$dates; if (length(K) < 3) next
  dup <- any(sapply(c(100,200), function(bb) sum(abs(K-bb) <= 30) >= 2))
  if (!dup) next
  base <- sse(K, Y)
  reo <- sapply(seq_along(K), function(j) {
    g <- art_kink(Y, "fixed", alpha=0, start=K[-j], max_knots=length(K)-1, thin=FALSE)
    sse(g$dates, Y) - base })
  rows[[length(rows)+1]] <- data.frame(rep=i, knots=paste(K, collapse=","),
     dist=paste(sapply(K, function(k) min(abs(k-c(100,200)))), collapse=","),
     reopt_cost=paste(round(reo,1), collapse=","), min_reopt=round(min(reo),1), alpha=round(f$alpha,1))
}
out <- do.call(rbind, rows); print(out, row.names=FALSE)
nu2 <- 2*N + 2*sqrt(2*N*log(Tn)) + 2*log(Tn)
cat(sprintf("duplicate runs: %d of 200; min re-optimized cost below alpha in %d of them; 2N-dof noise level 2N + 2 sqrt(2N log T) + 2 log T = %.1f\n",
            nrow(out), sum(out$min_reopt < out$alpha), nu2))
write.csv(out, "study_kink_traps.csv", row.names=FALSE)
