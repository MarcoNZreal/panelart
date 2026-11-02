library(panelart)
Tn <- 300; N <- 10; s <- 0.04; tt <- 1:Tn; nrep <- 200
tested <- list(); kinds <- character(0)
for (i in 1:nrep) {
  set.seed(23000 + i)
  g1 <- s * sample(c(-1,1), N, TRUE)
  Y <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g1) + matrix(rnorm(Tn*N), Tn, N)
  f0 <- art_kink(Y, "mc", nsim=500, thin=FALSE)          # search only
  K <- f0$dates
  near <- sapply(c(100,200), function(b) sum(abs(K-b) <= 30))
  far  <- sum(sapply(K, function(k) all(abs(k-c(100,200)) > 30)))
  kinds <- c(kinds, if (length(K) == 2 && all(near == 1)) "exact" else if (far > 0 && any(near >= 2)) "far+dup" else if (far > 0) "far" else if (any(near >= 2)) "dup" else "missing")
  if (length(K) >= 2 && any(diff(K) <= 30)) {
    f <- art_kink(Y, "mc", nsim=500, thin=FALSE, pair_test=TRUE, pair_level=0.01)
    if (is.matrix(f$pair_info)) tested[[length(tested)+1]] <- data.frame(rep=i, search_knots=paste(K, collapse=","),
        final=paste(f$dates, collapse=","), f$pair_info)
  }
}
cat("Search output (no thinning), 200 runs:\n"); print(table(kinds))
t <- do.call(rbind, tested); print(t, row.names=FALSE, digits=3)
write.csv(t, "study_kink_pair_diag.csv", row.names=FALSE)
cat(sprintf("runs with a close pair: %d; pair tests run: %d; merged: %d\n", length(unique(t$rep)), nrow(t), sum(t$merged)))
