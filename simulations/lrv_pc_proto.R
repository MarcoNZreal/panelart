library(panelart)
N <- 10; Tn <- 300; rho <- 0.5; phi <- 0.7; nrep <- 200; t0 <- proc.time()[["elapsed"]]
truth_top <- N*rho*(1+phi)/(1-phi) + (1-rho); truth_bulk <- (1-rho)
ar1 <- function(z, p) as.numeric(stats::filter(z, p, "recursive", init=0)) * sqrt(1-p^2)
parzen <- function(x, L) { x <- x - mean(x); n <- length(x); v <- sum(x^2)/n; w <- 1
  for (l in 1:L) { k <- panelart:::.kern(l/(L+1), "parzen"); v <- v + 2*k*sum(x[1:(n-l)]*x[(l+1):n])/n; w <- w + 2*k }
  v }
res <- matrix(NA, nrep, 6); colnames(res) <- c("full_top","pc_parzen_top","pc_ar_top","pc_max_top","bulk_pool_ar","bulk_pool_parzen")
for (i in 1:nrep) {
  set.seed(61000+i)
  f <- ar1(rnorm(Tn), phi) * sqrt(rho); E <- f %*% t(rep(1,N)) + matrix(rnorm(Tn*N),Tn,N) * sqrt(1-rho)
  Cs <- panelart:::.cross_cov(E); u <- eigen(Cs, symmetric=TRUE)$vectors[,1]; if (sum(u) < 0) u <- -u
  s <- drop(E %*% u); Eperp <- E - outer(s, u)
  Om <- panelart:::.kink_lrcov(E, integer(0), floor(Tn^(1/3)), 0L, 6L)
  res[i,1] <- eigen(Om, symmetric=TRUE, only.values=TRUE)$values[1]
  res[i,2] <- parzen(s, 12)
  res[i,3] <- panelart:::.lrv_ar_diff(matrix(s, ncol=1))
  res[i,4] <- max(res[i,2], res[i,3])
  res[i,5] <- panelart:::.lrv_ar_diff(Eperp[, 1:N]) * N/(N-1)   # idiosyncratic pooled (residual panel has N-1 dims)
  res[i,6] <- mean(apply(Eperp, 2, parzen, L = 8)) * N/(N-1)
}
q <- function(x, tr) round(c(median = median(x)/tr, q10 = quantile(x, .1)/tr, q90 = quantile(x, .9)/tr, mean = mean(x)/tr), 2)
out <- rbind(`full Parzen matrix, top eigenvalue (L = 6)` = q(res[,1], truth_top), `PC score, Parzen L=12` = q(res[,2], truth_top),
             `PC score, AR first differences` = q(res[,3], truth_top), `PC score, max of the two` = q(res[,4], truth_top),
             `bulk pooled, AR first differences` = q(res[,5], truth_bulk), `bulk pooled, Parzen L=8` = q(res[,6], truth_bulk))
print(out); write.csv(cbind(estimator = rownames(out), out), "study_lrv_pc_accuracy.csv", row.names = FALSE)
cat("true top LR eigenvalue", round(truth_top,2), "bulk", truth_bulk, "; sec", round(proc.time()[["elapsed"]]-t0), "\n")
