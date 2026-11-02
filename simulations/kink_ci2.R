library(panelart)
nrep <- as.integer(Sys.getenv("NREP","300")); Tn <- 300; N <- 10; sl <- 0.03
e <- matrix(NA_real_, nrep, 2); hit <- matrix(NA_real_, nrep, 2); sdm <- matrix(NA_real_, nrep, 2); t0 <- proc.time()[["elapsed"]]
for (i in seq_len(nrep)) {
  set.seed(19000 + i)
  d <- art_sim_data(Tn, N, breaks=c(1/3,2/3), delta_bar=0, slope_bar=sl)
  f <- art_kink(d$Y, "mc", nsim=200)
  if (f$n_breaks == 2) { ci <- art_ci(f); e[i,] <- f$dates - d$dates; sdm[i,] <- ci$sd
    hit[i,] <- as.numeric(d$dates >= ci$lower & d$dates <= ci$upper) }
  if (proc.time()[["elapsed"]]-t0 > 240) { cat("stopped early at rep", i, "\n"); break }
}
ok <- !is.na(e[,1]); cat("reps run:", sum(!is.na(sdm[,1])|TRUE), " exactly two knots:", round(mean(ok[seq_len(i)]),3), "\n")
cat("empirical sd per knot:", round(apply(e,2,sd,na.rm=TRUE),2), " robust:", round(apply(e,2,function(x) IQR(x,na.rm=TRUE)/1.349),2), " mean CI sd:", round(colMeans(sdm,na.rm=TRUE),2), "\n")
cat("coverage per knot:", round(colMeans(hit,na.rm=TRUE),3), "\n")
write.csv(data.frame(knot=1:2, empirical_sd=apply(e,2,sd,na.rm=TRUE), robust_sd=apply(e,2,function(x) IQR(x,na.rm=TRUE)/1.349), mean_ci_sd=colMeans(sdm,na.rm=TRUE), coverage95=colMeans(hit,na.rm=TRUE), found_two=mean(ok[seq_len(i)]), reps=i), "study_kink_ci2.csv", row.names=FALSE)
