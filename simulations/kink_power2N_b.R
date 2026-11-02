library(panelart)
nrep <- as.integer(Sys.getenv("NREP","150")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn
h1 <- pmax(tt-100,0); h2 <- pmax(tt-200,0)
rows <- list()
for (s2 in c(0.004, 0.005, 0.006, 0.008, 0.012)) {
  gain2 <- N * s2^2 * sum(qr.resid(qr(cbind(1,tt,h1)), h2)^2)
  det <- cnt <- matrix(NA, nrep, 2); strong <- matrix(NA, nrep, 2)
  for (i in seq_len(nrep)) {
    set.seed(63000 + i)
    sg <- sample(c(-1,1), N, TRUE)
    Y <- outer(h1, 0.04*sg) + outer(h2, s2*sg) + matrix(rnorm(Tn*N), Tn, N)
    for (v in 1:2) {
      f <- art_kink(Y, "mc", nsim=500, alpha_thin = if (v == 1) NULL else "2N")
      d <- f$dates
      det[i, v] <- any(abs(d - 200) <= 30); strong[i, v] <- any(abs(d - 100) <= 30); cnt[i, v] <- length(d)
    }
  }
  for (v in 1:2) rows[[length(rows)+1]] <- data.frame(weak_slope=s2, weak_gain=round(gain2,1),
      thinning=if (v==1) "alpha" else "2N level", weak_detected=mean(det[,v]), strong_detected=mean(strong[,v]),
      exact_two=mean(cnt[,v]==2), extra=mean(cnt[,v]>2))
  cat("weak slope", s2, "done", round(proc.time()[["elapsed"]]-t0), "s\n")
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_power2N_b.csv", row.names=FALSE)
options(width=160); print(r, row.names=FALSE, digits=3)
