library(panelart)
phi <- as.numeric(Sys.getenv("PHI","0.5")); nrep <- as.integer(Sys.getenv("NREP","100")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn; b <- 150; s <- 0.08
om2 <- (1+phi)/(1-phi)
hinge <- pmax(tt-b,0)
hres <- qr.resid(qr(cbind(1,tt)), hinge); Rstep <- qr.resid(qr(cbind(1,tt,hinge)), as.numeric(tt > b))
g2 <- N*s^2
sd_iid <- sqrt(1/(g2*sum(Rstep^2)) * (1 + N/(g2*sum(hres^2))))            # what the i.i.d. formula would give with sigma^2 = 1
sd_dep <- sqrt(om2/(g2*sum(Rstep^2)) * (1 + N*om2/(g2*sum(hres^2))))      # with the long-run variance
est <- cov_dep <- cov_naive <- sdhat <- numeric(0); nk <- integer(nrep)
for (i in seq_len(nrep)) {
  set.seed(88000+i); g <- s*sample(c(-1,1),N,TRUE)
  E <- apply(matrix(rnorm(Tn*N),Tn,N), 2, function(z) as.numeric(stats::filter(z, phi, "recursive"))) * sqrt(1-phi^2)
  Y <- outer(hinge, g) + E
  f <- art_kink(Y, "hac_mc", nsim=500); nk[i] <- length(f$dates)
  if (length(f$dates) == 1L) {
    ci <- art_ci(f, 0.95); est <- c(est, f$dates)
    cov_dep <- c(cov_dep, ci$lower <= b & b <= ci$upper); sdhat <- c(sdhat, ci$sd)
    f2 <- f; f2$omega2 <- NULL; ci2 <- art_ci(f2, 0.95); cov_naive <- c(cov_naive, ci2$lower <= b & b <= ci2$upper)
  }
}
r <- data.frame(phi=phi, reps=nrep, one_knot=mean(nk==1), emp_sd=sd(est-b), emp_rmse=sqrt(mean((est-b)^2)), mean_error=mean(est-b),
   pred_sd_longrun=sd_dep, pred_sd_iid_formula=sd_iid, mean_est_sd=mean(sdhat),
   coverage_longrun_variance=mean(cov_dep), coverage_naive_residual_variance=mean(cov_naive))
write.csv(r, sprintf("study_kink_dep2_phi%s.csv", phi), row.names=FALSE)
options(width=190); print(r, row.names=FALSE, digits=3); cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
