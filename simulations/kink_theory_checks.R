library(panelart)
Tn <- 300; N <- 10; tt <- 1:Tn
## (1) identity: removal cost of a knot = ||gamma_hat||^2 * ||h~||^2 (h~ = hinge residualized on 1, t, other hinges),
##     and the local bound ||h~||^2 <= L^2 R^2 / (3 (L + R)) with L, R the distances to the neighbouring knots (or ends)
set.seed(5); rel <- numeric(0); slack <- numeric(0)
for (i in 1:40) {
  d <- art_sim_data(Tn, N, breaks=c(1/3,2/3), delta_bar=0, slope_bar=0.04)
  f <- art_kink(d$Y, "mc", nsim=200, thin=FALSE)
  K <- f$dates; if (length(K) < 2) next
  B <- function(K) cbind(1, tt, sapply(K, function(k) pmax(tt-k,0)))
  sse <- function(K) sum(qr.resid(qr(B(K)), d$Y)^2)
  for (j in seq_along(K)) {
    cost <- sse(K[-j]) - sse(K)
    hj <- qr.resid(qr(B(K[-j])), pmax(tt-K[j],0))
    g2 <- sum(f$gamma[j,]^2)
    rel <- c(rel, abs(cost - g2*sum(hj^2))/cost)
    Lj <- K[j] - c(0,K)[j]; Rj <- c(K,Tn)[j+1] - K[j]
    slack <- c(slack, (Lj^2*Rj^2/(3*(Lj+Rj)))/sum(hj^2))
  }
}
cat("identity: max relative difference between removal cost and ||gamma||^2 ||h||^2 =", signif(max(rel),3), " over", length(rel), "knots\n")
cat("local bound L^2R^2/(3(L+R)) divided by ||h||^2: min", round(min(slack),2), " median", round(median(slack),2), " (must be >= 1)\n")
## (2) the theory's sufficient condition for covering versus the simulated design
x <- (3+3)*log(Tn); sprime <- 6*sqrt(N*x) + 64*x; lambda <- 1/3
need <- 6 * sprime * 2 * (2/lambda)^6   # required ||gamma||^2 T^3 for p = 2 comparable knots, w = lambda T / 2
have <- N*0.04^2*Tn^3
cat(sprintf("covering theorem (m-bar = 3, C0' = 64): needs ||gamma||^2 T^3 >= %.2e ; simulated design has %.2e (ratio %.0f)\n", need, have, need/have))
## (3) conservative (union-bound) penalty versus the Monte Carlo penalty, two kinks, 200 replications
nrep <- 200
far  <- function(d, truth=c(100,200), r=30) any(sapply(d, function(k) all(abs(k-truth)>r)))
res <- list()
for (s in c(0, 0.03, 0.04)) for (pen in c("mc","iid")) {
  ex <- far_ <- both <- numeric(nrep)
  for (i in 1:nrep) {
    set.seed(41000 + i)
    g1 <- s * sample(c(-1,1), N, TRUE)
    Y <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g1) + matrix(rnorm(Tn*N), Tn, N)
    f <- art_kink(Y, pen, nsim=500)
    truth_n <- if (s == 0) 0 else 2
    ex[i] <- length(f$dates)==truth_n; far_[i] <- if (length(f$dates)) far(f$dates, if (s==0) c(-1e6) else c(100,200)) else FALSE
    both[i] <- if (s==0) NA else (any(abs(f$dates-100)<=5) && any(abs(f$dates-200)<=5))
  }
  res[[length(res)+1]] <- data.frame(slope=s, penalty=pen, alpha=round(f$alpha,1), exact_count=mean(ex), far_knot_runs=mean(far_), both_found=mean(both))
}
r <- do.call(rbind, res); write.csv(r, "study_kink_union_vs_mc.csv", row.names=FALSE); print(r, row.names=FALSE)
