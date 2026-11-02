library(panelart)
N <- 10; gap <- 6L; level <- 0.05; t0 <- proc.time()[["elapsed"]]
ar1 <- function(M, p) if (p == 0) M else apply(M, 2, function(z) as.numeric(stats::filter(z, p, "recursive", init=0)))
# stationary-variance-1 AR(1) with burn-in: simulate T+200 and drop
gen <- list(
  indep = function(Tn, p, rho = 0, ...) { Z <- matrix(rnorm((Tn+200)*N), Tn+200, N); ar1(Z, p)[201:(Tn+200), ] * sqrt(1-p^2) },
  equi  = function(Tn, p, rho, ...) { Z <- matrix(rnorm((Tn+200)*N), Tn+200, N); (ar1(Z, p)[201:(Tn+200), ] * sqrt(1-p^2)) %*% chol((1-rho)*diag(N) + rho*matrix(1,N,N)) },
  nonsep = function(Tn, p, rho, ...) { f <- ar1(matrix(rnorm(Tn+200), ncol=1), p)[201:(Tn+200), 1] * sqrt(1-p^2) * sqrt(rho)
              f %*% t(rep(1,N)) + matrix(rnorm(Tn*N), Tn, N) * sqrt(1-rho) })
true_q <- function(Tn, type, p, rho, nsim, seed) {
  set.seed(seed); tt <- seq_len(Tn); ks <- gap:(Tn-gap)
  H <- outer(tt, ks, function(t,k) pmax(t-k,0)); q <- qr(cbind(1,tt)); Hres <- qr.resid(q, H); nrm <- colSums(Hres^2)
  mx <- numeric(nsim)
  for (i in seq_len(nsim)) { Y <- gen[[type]](Tn, p, rho); P <- crossprod(Hres, qr.resid(q, Y))^2; mx[i] <- max(rowSums(P)/nrm) }
  unname(quantile(mx, 1 - level)) }
lam <- function(type, p, rho) { om <- (1+p)/(1-p)
  if (type == "indep") rep(om, N) else if (type == "equi") om * eigen((1-rho)*diag(N) + rho*matrix(1,N,N), symmetric=TRUE, only.values=TRUE)$values
  else eigen(rho*om*matrix(1,N,N) + (1-rho)*diag(N), symmetric=TRUE, only.values=TRUE)$values }
cases <- list(list("indep",0,0,"iid"), list("indep",.3,0,"AR(1) 0.3"), list("indep",.5,0,"AR(1) 0.5"), list("indep",.7,0,"AR(1) 0.7"),
              list("indep",.8,0,"AR(1) 0.8"), list("indep",.9,0,"AR(1) 0.9"), list("equi",0,.5,"equicorr 0.5, iid time"),
              list("equi",.5,.5,"equicorr 0.5, AR(1) 0.5 (separable)"), list("nonsep",.7,.5,"factor AR(1) 0.7 + white"))
Tn <- as.integer(Sys.getenv("TN","300")); nsim <- as.integer(Sys.getenv("NSIM","2000")); rows <- list()
for (cs in cases) {
  qt <- true_q(Tn, cs[[1]], cs[[2]], cs[[3]], nsim, 7001)
  qm <- art_penalty_kink_mc_w(lam(cs[[1]], cs[[2]], cs[[3]]), Tn, level, nsim, gap, seed = 8001)
  rows[[length(rows)+1]] <- data.frame(model = cs[[4]], T = Tn, true_quantile = qt, iid_mc_quantile = qm, ratio = qt/qm)
  cat(cs[[4]], round(qt,1), round(qm,1), round(qt/qm,3), "| sec", round(proc.time()[["elapsed"]]-t0), "\n"); flush.console() }
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_calib1_T%d.csv", Tn), row.names=FALSE)
