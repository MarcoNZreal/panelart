library(panelart)
nrep <- as.integer(Sys.getenv("NREP","120"))
cfg <- list(c(200,10,0.04), c(200,10,0.02), c(200,50,0.015), c(400,10,0.015))
out <- list(); t0 <- proc.time()[["elapsed"]]
for (cf in cfg) {
  Tn <- cf[1]; N <- cf[2]; sl <- cf[3]
  e <- rep(NA_real_, nrep); hit <- rep(NA_real_, nrep); sdp <- rep(NA_real_, nrep); w <- rep(NA_real_, nrep)
  for (i in seq_len(nrep)) {
    set.seed(15000 + i + 31*N + round(1000*sl) + Tn)
    d <- art_sim_data(Tn, N, breaks=0.5, delta_bar=0, slope_bar=sl)
    f <- art_kink(d$Y, "mc", nsim=300)
    if (f$n_breaks == 1) { ci <- art_ci(f); e[i] <- f$dates - d$dates; sdp[i] <- ci$sd
      hit[i] <- as.numeric(d$dates >= ci$lower && d$dates <= ci$upper); w[i] <- ci$upper-ci$lower }
  }
  gam2 <- N*sl^2; tau <- 0.5
  theory_sd <- sqrt(4/(gam2*Tn*tau*(1-tau)))   # sigma = 1
  out[[length(out)+1]] <- data.frame(T=Tn, N=N, slope=sl, gamma_ratio=gam2*Tn^3, per_series=sl^2*Tn^3,
    found_one=mean(!is.na(e)), theory_sd=theory_sd, empirical_sd=sd(e,na.rm=TRUE), mean_estimated_sd=mean(sdp,na.rm=TRUE),
    coverage95=mean(hit,na.rm=TRUE), mean_width=mean(w,na.rm=TRUE))
  cat("T",Tn,"N",N,"slope",sl,"done",round(proc.time()[["elapsed"]]-t0),"s\n")
}
r <- do.call(rbind, out); write.csv(r, "study_kink_ci.csv", row.names=FALSE)
options(width=170); print(round(r,3), row.names=FALSE)
