library(panelart)
nrep <- as.integer(Sys.getenv("NREP","100")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn; b <- 150
h <- pmax(tt-b,0); hres <- qr.resid(qr(cbind(1,tt)), h)
rows <- list()
for (s in c(0.004, 0.005, 0.006, 0.008)) {
  gain <- N * s^2 * sum(hres^2)
  det <- matrix(NA, nrep, 2); cnt <- matrix(NA, nrep, 2)
  for (i in seq_len(nrep)) {
    set.seed(61000 + i)
    g <- s * sample(c(-1,1), N, TRUE)
    Y <- outer(h, g) + matrix(rnorm(Tn*N), Tn, N)
    for (v in 1:2) {
      f <- art_kink(Y, "mc", nsim=500, alpha_thin = if (v == 1) NULL else "2N")
      d <- f$dates
      det[i, v] <- any(abs(d - b) <= 30); cnt[i, v] <- length(d)
    }
  }
  for (v in 1:2) rows[[length(rows)+1]] <- data.frame(slope=s, population_gain=round(gain,1),
      thinning=if (v==1) "alpha" else "2N level", detect=mean(det[,v]), exact_count=mean(cnt[,v]==1), extra=mean(cnt[,v]>1))
  cat("slope", s, "done", round(proc.time()[["elapsed"]]-t0), "s\n")
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_power2N.csv", row.names=FALSE)
options(width=150); print(r, row.names=FALSE, digits=3)
