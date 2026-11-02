library(panelart)
nrep <- as.integer(Sys.getenv("NREP","300")); Tn <- 200; pooled <- 32000
out <- list(); t0 <- proc.time()[["elapsed"]]
for (N in c(10, 50, 200, 500)) {
  sl <- sqrt(pooled/(N*Tn^3))
  e <- rep(NA_real_, nrep); hit <- rep(NA_real_, nrep); sdj <- rep(NA_real_, nrep)
  for (i in seq_len(nrep)) {
    set.seed(17000 + i + N)
    d <- art_sim_data(Tn, N, breaks=0.5, delta_bar=0, slope_bar=sl)
    f <- art_kink(d$Y, "mc", nsim=200)
    if (f$n_breaks == 1) { ci <- art_ci(f); e[i] <- f$dates - d$dates; sdj[i] <- ci$sd; hit[i] <- as.numeric(abs(d$dates - f$dates) <= qnorm(0.975)*ci$sd) }
  }
  th <- sqrt(4/(N*sl^2*Tn*0.25))
  out[[length(out)+1]] <- data.frame(N=N, slope=round(sl,4), per_series_g2T3=N*0+sl^2*Tn^3, found_one=mean(!is.na(e)),
     theory_sd=th, empirical_sd=sd(e,na.rm=TRUE), robust_sd=IQR(e,na.rm=TRUE)/1.349, ratio=sd(e,na.rm=TRUE)/th, predicted_ratio=sqrt(1+192*N/pooled), mean_ci_sd=mean(sdj,na.rm=TRUE), coverage95=mean(hit,na.rm=TRUE))
  cat("N", N, "done", round(proc.time()[["elapsed"]]-t0), "s\n")
}
r <- do.call(rbind, out); write.csv(r, "study_kink_perseries.csv", row.names=FALSE)
options(width=150); print(round(r,3), row.names=FALSE)
