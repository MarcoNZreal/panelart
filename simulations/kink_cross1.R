library(panelart)
N <- as.integer(Sys.getenv("N","10")); Tn <- as.integer(Sys.getenv("TN","300")); rho <- as.numeric(Sys.getenv("RHO","0.5")); phi <- as.numeric(Sys.getenv("PHI","0"))
nonsep <- Sys.getenv("NONSEP","0") == "1"; s <- as.numeric(Sys.getenv("S","0")); orth <- Sys.getenv("ORTH","0") == "1"
nrep <- as.integer(Sys.getenv("NREP","80")); modes <- strsplit(Sys.getenv("MODES","none,separable,whiten"), ",")[[1]]; tag <- Sys.getenv("TAG","x"); label <- Sys.getenv("LABEL","")
tt <- 1:Tn; kn <- if (s > 0) c(100,200) else numeric(0); t0 <- proc.time()[["elapsed"]]
Om <- (1-rho)*diag(N) + rho*matrix(1,N,N); Lc <- chol(Om)
ar1 <- function(M, p) if (p == 0) M else apply(M, 2, function(z) as.numeric(stats::filter(z, p, "recursive", init=0))) * sqrt(1-p^2)
gen <- function(seed) {
  set.seed(seed); sg <- if (orth) rep(c(1,-1), length.out=N) else rep(1,N); g <- s*sg
  mu <- if (length(kn)) Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g))) else matrix(0,Tn,N)
  E <- if (!nonsep) ar1(matrix(rnorm(Tn*N),Tn,N), phi) %*% Lc
       else { f <- ar1(matrix(rnorm(Tn),Tn,1), phi) * sqrt(rho); u <- matrix(rnorm(Tn*N),Tn,N) * sqrt(1-rho); f %*% t(rep(1,N)) + u }
  mu + E }
rows <- list()
for (m in modes) {
  ex <- miss <- fp <- nk <- numeric(0)
  for (i in seq_len(nrep)) {
    if (proc.time()[["elapsed"]] - t0 > 255) break
    f <- art_kink(gen(41000+i), "hac_mc", nsim=as.integer(Sys.getenv("NSIM","500")), cross=m, max_knots=10, pc_lrv=Sys.getenv("PCLRV","ar"), pc_scale=as.numeric(Sys.getenv("PCSCALE","1")), pc_floor=Sys.getenv("PCFLOOR","1")=="1", lrv_extra=as.integer(Sys.getenv("LRVEXTRA","0"))); K <- f$dates; nk <- c(nk, length(K))
    if (length(kn)) { near <- sapply(kn, function(b) sum(abs(K-b) <= 30)); ex <- c(ex, length(K)==2 && all(near==1)); miss <- c(miss, any(near==0)) } else fp <- c(fp, length(K) > 0) }
  rows[[length(rows)+1]] <- data.frame(setting=label, N=N, T=Tn, slope=s, jumps=if (s>0) (if (orth) "orthogonal" else "aligned") else "-", cross=m, reps=length(nk),
     both_found=if (length(kn)) mean(ex) else NA, missed_window=if (length(kn)) mean(miss) else NA, false_positive=if (!length(kn)) mean(fp) else NA, mean_knots=mean(nk)) }
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_cross1_%s.csv", tag), row.names=FALSE)
options(width=170); print(r, row.names=FALSE, digits=3); cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
