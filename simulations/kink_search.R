library(panelart)
nrep <- as.integer(Sys.getenv("NREP","30")); Tn <- 300; N <- 10; t0 <- proc.time()[["elapsed"]]
## E2: repartition from perturbed starts (no insertions: max_knots = 2), slope change 0.04 at 100 and 200
rows <- list()
for (pat in c("toward", "same")) for (off in c(0, 10, 20, 40)) {
  err <- matrix(NA_real_, nrep, 2); nk <- rep(NA_real_, nrep)
  for (i in seq_len(nrep)) {
    set.seed(21000 + i)
    d <- art_sim_data(Tn, N, breaks=c(1/3,2/3), delta_bar=0, slope_bar=0.04)
    st <- if (pat == "toward") d$dates + c(off, -off) else d$dates + c(off, off)
    f <- art_kink(d$Y, "mc", nsim=200, start=st, max_knots=2)
    nk[i] <- f$n_breaks
    if (f$n_breaks == 2) err[i,] <- f$dates - d$dates
  }
  rows[[length(rows)+1]] <- data.frame(pattern=pat, offset=off, two_knots=mean(nk==2), mean_abs_err=mean(abs(err),na.rm=TRUE),
                                       sd_err=sd(as.vector(err),na.rm=TRUE))
}
e2 <- do.call(rbind, rows); write.csv(e2, "study_kink_repartition.csv", row.names=FALSE)
cat("E2 repartition from perturbed starts (T=300, N=10, slope change 0.04 at 100 and 200; predicted sd about 2 dates)\n")
print(e2, digits=3, row.names=FALSE); cat("elapsed", round(proc.time()[["elapsed"]]-t0), "s\n")
## E3: offsetting kinks (tent) vs same-direction kinks, from scratch, no cap
rows <- list(); tt <- 1:Tn
for (sg in c(1, -1)) for (s in c(0.02, 0.03, 0.04)) {
  nk <- numeric(nrep); f1 <- numeric(nrep); f2 <- numeric(nrep)
  for (i in seq_len(nrep)) {
    set.seed(23000 + i)
    g1 <- s * sample(c(-1,1), N, TRUE); g2 <- sg * g1
    mu <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g2)
    Y <- mu + matrix(rnorm(Tn*N), Tn, N)
    f <- art_kink(Y, "mc", nsim=200)
    nk[i] <- f$n_breaks; f1[i] <- any(abs(f$dates-100) <= 5); f2[i] <- any(abs(f$dates-200) <= 5)
  }
  rows[[length(rows)+1]] <- data.frame(pattern=if (sg==1) "same direction (staircase of slopes)" else "opposite (tent)", slope=s,
     pooled_each=N*s^2*Tn^3, mean_knots=mean(nk), exact_two=mean(nk==2), found_100=mean(f1), found_200=mean(f2), max_knots=max(nk))
}
e3 <- do.call(rbind, rows); write.csv(e3, "study_kink_tent.csv", row.names=FALSE)
cat("\nE3 same-direction versus offsetting kinks (from scratch)\n"); print(e3, digits=3, row.names=FALSE)
cat("elapsed", round(proc.time()[["elapsed"]]-t0), "s\n")
