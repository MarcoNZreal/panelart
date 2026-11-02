library(panelart)
nrep <- as.integer(Sys.getenv("NREP","150")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn
h1 <- pmax(tt-100,0); h2 <- pmax(tt-200,0)
meth <- list(alpha=list(), twoN=list(alpha_thin="2N"), pair01=list(pair_test=TRUE, pair_level=0.01))
rows <- list()
for (s2 in c(0.006, 0.008, 0.012)) {
  gain2 <- N * s2^2 * sum(qr.resid(qr(cbind(1,tt,h1)), h2)^2)
  det <- cnt <- matrix(NA, nrep, length(meth))
  for (i in seq_len(nrep)) {
    set.seed(63000 + i)
    sg <- sample(c(-1,1), N, TRUE)
    Y <- outer(h1, 0.04*sg) + outer(h2, s2*sg) + matrix(rnorm(Tn*N), Tn, N)
    for (v in seq_along(meth)) {
      f <- do.call(art_kink, c(list(Y=Y, penalty="mc", nsim=500), meth[[v]]))
      det[i, v] <- any(abs(f$dates - 200) <= 30); cnt[i, v] <- length(f$dates)
    }
  }
  for (v in seq_along(meth)) rows[[length(rows)+1]] <- data.frame(weak_slope=s2, weak_gain=round(gain2,1),
      method=names(meth)[v], weak_detected=mean(det[,v]), exact_two=mean(cnt[,v]==2), extra=mean(cnt[,v]>2))
  cat("weak slope", s2, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console()
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_pair2.csv", row.names=FALSE)
options(width=150); print(r, row.names=FALSE, digits=3)
