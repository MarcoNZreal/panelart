library(panelart)
nrep <- as.integer(Sys.getenv("NREP","300")); t0 <- proc.time()[["elapsed"]]
cfgs <- list(list(name="2 kinks, T=300, s=0.04", T=300, kn=c(100,200), s=0.04),
             list(name="2 kinks, T=300, s=0.06", T=300, kn=c(100,200), s=0.06),
             list(name="3 kinks, T=400, s=0.05", T=400, kn=c(100,200,300), s=0.05))
rows <- list()
for (cf in cfgs) {
  Tn <- cf$T; N <- 10; tt <- 1:Tn; m <- length(cf$kn)
  extras <- miss <- numeric(nrep); Rb <- Ra <- numeric(nrep)
  for (i in seq_len(nrep)) {
    set.seed(47000 + i)
    g <- cf$s * sample(c(-1,1), N, TRUE)
    mu <- Reduce(`+`, lapply(cf$kn, function(b) outer(pmax(tt-b,0), g)))
    Y <- mu + matrix(rnorm(Tn*N), Tn, N)
    f <- art_kink(Y, "iid")
    Rfun <- function(K) { B <- cbind(1, tt, if (length(K)) sapply(K, function(k) pmax(tt-k,0)) else NULL); sum(qr.resid(qr(B), mu)^2) }
    Rb[i] <- Rfun(f$dates_before_thin); Ra[i] <- Rfun(f$dates)
    extras[i] <- length(f$dates) > m
    miss[i] <- any(sapply(cf$kn, function(b) all(abs(f$dates - b) > 30)))
  }
  x <- (floor(1/0.02) + 3) * log(Tn); sprime <- 6*sqrt(N*x) + 64*x
  Gam <- m * N * cf$s^2
  rho0 <- sqrt(2 * sprime * m * Tn^3 * Gam / 192)
  rows[[length(rows)+1]] <- data.frame(design=cf$name, reps=nrep, extras=mean(extras), window_missed=mean(miss),
    misfit_median=median(Ra), misfit_max=max(Ra), misfit_max_before_thin=max(Rb), bound_rho0=rho0)
  cat(cf$name, "done", round(proc.time()[["elapsed"]]-t0), "s\n")
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_extras_iid300.csv", row.names=FALSE)
options(width=170); print(r, row.names=FALSE, digits=4)
