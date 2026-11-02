library(panelart)
Tn <- as.integer(Sys.getenv("TN","300")); phi <- as.numeric(Sys.getenv("PHI","0.8")); nrep <- as.integer(Sys.getenv("NREP","100"))
N <- 10; t0 <- proc.time()[["elapsed"]]; fp <- numeric(0); om <- numeric(0); npass <- numeric(0)
for (i in seq_len(nrep)) {
  if (proc.time()[["elapsed"]] - t0 > 230) break
  set.seed(55000+i)
  E <- apply(matrix(rnorm(Tn*N),Tn,N), 2, function(z) as.numeric(stats::filter(z, phi, "recursive"))) * sqrt(1-phi^2)
  f <- art_kink(E, "hac_mc", nsim=500, max_knots=10, max_pass=6L); fp <- c(fp, length(f$dates) > 0); om <- c(om, f$omega2); npass <- c(npass, nrow(f$history)) }
r <- data.frame(T=Tn, phi=phi, reps_done=length(fp), false_positive=mean(fp), omega2_hat=mean(om), omega2_true=(1+phi)/(1-phi), mean_passes=mean(npass))
write.csv(r, sprintf("study_kink_dep5_T%d_phi%s.csv", Tn, phi), row.names=FALSE); options(width=140); print(r, row.names=FALSE, digits=3)
cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
