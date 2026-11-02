library(panelart)
N <- 10; Tn <- 300; tt <- 1:Tn; rho <- 0.5; phi <- 0.7; nrep <- as.integer(Sys.getenv("NREP","40")); t0 <- proc.time()[["elapsed"]]
OmLR <- rho*(1+phi)/(1-phi)*matrix(1,N,N) + (1-rho)*diag(N)
alpha_or <- as.numeric(Sys.getenv("ALPHA","256.8"))
cat("oracle penalty used:", round(alpha_or,1), "\n")
ar1 <- function(z, p) as.numeric(stats::filter(z, p, "recursive", init=0)) * sqrt(1-p^2)
gen <- function(seed, s, orth) { set.seed(seed); g <- s*(if (orth) rep(c(1,-1), length.out=N) else rep(1,N))
  mu <- if (s > 0) outer(pmax(tt-100,0), g) + outer(pmax(tt-200,0), g) else matrix(0,Tn,N)
  f <- ar1(rnorm(Tn), phi) * sqrt(rho); mu + f %*% t(rep(1,N)) + matrix(rnorm(Tn*N),Tn,N) * sqrt(1-rho) }
rows <- list()
for (cfg in list(list(s=0, orth=FALSE, n=200), list(s=0.05, orth=FALSE, n=nrep), list(s=0.05, orth=TRUE, n=nrep))) {
  ex <- miss <- fp <- nk <- numeric(0)
  for (i in seq_len(cfg$n)) {
    if (proc.time()[["elapsed"]] - t0 > 250) break
    f <- art_kink(gen(41000+i, cfg$s, cfg$orth), "fixed", alpha=alpha_or, sigma2=1, max_knots=10); K <- f$dates; nk <- c(nk, length(K))
    if (cfg$s > 0) { near <- sapply(c(100,200), function(b) sum(abs(K-b) <= 30)); ex <- c(ex, length(K)==2 && all(near==1)); miss <- c(miss, any(near==0)) } else fp <- c(fp, length(K) > 0) }
  rows[[length(rows)+1]] <- data.frame(slope=cfg$s, jumps=if (cfg$s==0) "-" else if (cfg$orth) "orthogonal" else "aligned", cross="oracle", reps=length(nk),
     both_found=if (cfg$s>0) mean(ex) else NA, missed_window=if (cfg$s>0) mean(miss) else NA, false_positive=if (cfg$s==0) mean(fp) else NA, mean_knots=mean(nk)) }
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_cross1_oracle_nonsep_alpha%s.csv", Sys.getenv("ALPHA","256.8")), row.names=FALSE); options(width=150); print(r, row.names=FALSE, digits=3); cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
