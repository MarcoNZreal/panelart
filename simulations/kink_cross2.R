library(panelart)
N <- 10; Tn <- 300; tt <- 1:Tn; b <- 150; s <- as.numeric(Sys.getenv("S","0.08")); rho <- 0.5
phi <- as.numeric(Sys.getenv("PHI","0")); nonsep <- Sys.getenv("NONSEP","0") == "1"; orth <- Sys.getenv("ORTH","0") == "1"; nrep <- as.integer(Sys.getenv("NREP","150"))
modes <- strsplit(Sys.getenv("MODES","none,separable,whiten"), ",")[[1]]; tag <- Sys.getenv("TAG","x"); t0 <- proc.time()[["elapsed"]]
Om <- (1-rho)*diag(N) + rho*matrix(1,N,N); Lc <- chol(Om); omt <- (1+phi)/(1-phi); OmLR <- if (nonsep) rho*omt*matrix(1,N,N) + (1-rho)*diag(N) else omt * Om
g <- s * (if (orth) rep(c(1,-1), N/2) else rep(1, N)); hinge <- pmax(tt-b,0)
hres <- qr.resid(qr(cbind(1,tt)), hinge); Rstep <- qr.resid(qr(cbind(1,tt,hinge)), as.numeric(tt > b))
sd_sand <- sqrt(drop(t(g) %*% OmLR %*% g) / (sum(g^2)^2 * sum(Rstep^2)) * (1 + sum(diag(OmLR)) / (sum(g^2) * sum(hres^2))))
gt <- drop(crossprod(g, solve(Om, g))); sd_white <- if (nonsep) NA else sqrt(omt / (gt * sum(Rstep^2)) * (1 + N * omt / (gt * sum(hres^2))))
sd_scalar <- sqrt(omt / (sum(g^2) * sum(Rstep^2)) * (1 + N * omt / (sum(g^2) * sum(hres^2))))   # what the scalar formula gives with the true per-series variance 1
ar1 <- function(M, p) if (p == 0) M else apply(M, 2, function(z) as.numeric(stats::filter(z, p, "recursive", init=0))) * sqrt(1-p^2)
rows <- list()
for (m in modes) {
  est <- cov <- sdh <- numeric(0); nk <- integer(0)
  for (i in seq_len(nrep)) {
    if (proc.time()[["elapsed"]] - t0 > 255) break
    set.seed(52000+i)
    E <- if (!nonsep) ar1(matrix(rnorm(Tn*N),Tn,N), phi) %*% Lc else { f <- ar1(matrix(rnorm(Tn),Tn,1), phi) * sqrt(rho); f %*% t(rep(1,N)) + matrix(rnorm(Tn*N),Tn,N) * sqrt(1-rho) }
    Y <- outer(hinge, g) + E
    f <- art_kink(Y, "hac_mc", nsim=500, cross=m, max_knots=10, pc_lrv=Sys.getenv("PCLRV","ar"), lrv_extra=as.integer(Sys.getenv("LRVEXTRA","0"))); nk <- c(nk, length(f$dates))
    if (length(f$dates) == 1L) { ci <- art_ci(f, 0.95); est <- c(est, f$dates); cov <- c(cov, ci$lower <= b & b <= ci$upper); sdh <- c(sdh, ci$sd) } }
  rows[[length(rows)+1]] <- data.frame(jumps=if (orth) "orthogonal" else "aligned", phi=phi, cross=m, reps=length(nk), one_knot=mean(nk==1),
     emp_sd=sd(est-b), pred_sd_sandwich=sd_sand, pred_sd_whitened=if (m=="whiten") sd_white else NA, pred_sd_scalar=sd_scalar, mean_est_sd=mean(sdh), coverage=mean(cov)) }
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_cross2_%s.csv", tag), row.names=FALSE)
options(width=190); print(r, row.names=FALSE, digits=3); cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
