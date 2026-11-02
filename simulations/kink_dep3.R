library(panelart)
nrep <- as.integer(Sys.getenv("NREP","100")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn
rows <- list()
for (rng in list(c(0,0.4), c(0,0.8), c(0.5,0.9))) {
  fp <- numeric(nrep); om <- numeric(nrep)
  for (i in seq_len(nrep)) {
    set.seed(99000+i); phis <- runif(N, rng[1], rng[2])
    E <- sapply(seq_len(N), function(j) as.numeric(stats::filter(rnorm(Tn), phis[j], "recursive")) * sqrt(1-phis[j]^2))
    f <- art_kink(E, "hac_mc", nsim=500); fp[i] <- length(f$dates) > 0; om[i] <- f$omega2
  }
  trueom <- mean((1+seq(rng[1],rng[2],length.out=1000))/(1-seq(rng[1],rng[2],length.out=1000)))
  rows[[length(rows)+1]] <- data.frame(phi_range=sprintf("U(%.1f,%.1f)", rng[1], rng[2]), reps=nrep, false_positive=mean(fp),
        omega2_hat=mean(om), omega2_true_mean=trueom)
  cat(rows[[length(rows)]]$phi_range, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console()
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_dep3.csv", row.names=FALSE); options(width=140); print(r, row.names=FALSE, digits=3)
