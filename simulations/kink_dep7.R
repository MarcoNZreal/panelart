library(panelart)
ar <- as.numeric(strsplit(Sys.getenv("AR",""), ",")[[1]]); ma <- as.numeric(strsplit(Sys.getenv("MA",""), ",")[[1]])
s <- as.numeric(Sys.getenv("S","0.06")); nrep <- as.integer(Sys.getenv("NREP","40")); modes <- strsplit(Sys.getenv("MODES","hac,ar_start,ar"), ",")[[1]]
label <- Sys.getenv("LABEL","AR1 0.7"); tag <- Sys.getenv("TAG","x")
Tn <- 300; N <- 10; tt <- 1:Tn; kn <- if (s > 0) c(100,200) else numeric(0); t0 <- proc.time()[["elapsed"]]
psi <- if (length(ar) || length(ma)) ARMAtoMA(ar, ma, 600) else numeric(0); vr <- 1 + sum(psi^2)
tr <- ((1 + sum(ma)) / (1 - sum(ar)))^2 / vr            # true long-run variance at unit marginal variance
gen <- function(seed) { set.seed(seed); g <- s*sample(c(-1,1),N,TRUE)
  mu <- if (length(kn)) Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g))) else matrix(0,Tn,N)
  E <- if (length(ar) == 1 && !length(ma)) apply(matrix(rnorm(Tn*N),Tn,N), 2, function(z) as.numeric(stats::filter(z, ar, "recursive", init=rep(0,1)))) * sqrt(1-ar^2)
       else sapply(1:N, function(j) as.numeric(arima.sim(list(ar=ar, ma=ma), n=Tn, n.start=200)) / sqrt(vr))
  mu + E }
rows <- list()
for (m in modes) {
  ex <- miss <- fp <- om <- nk <- numeric(0)
  for (i in seq_len(nrep)) {
    if (proc.time()[["elapsed"]] - t0 > 255) break
    f <- art_kink(gen(77000+i), "hac_mc", nsim=500, lrv=m, max_knots=10); K <- f$dates; om <- c(om, f$omega2); nk <- c(nk, length(K))
    if (length(kn)) { near <- sapply(kn, function(b) sum(abs(K-b) <= 30)); ex <- c(ex, length(K)==2 && all(near==1)); miss <- c(miss, any(near==0)) } else fp <- c(fp, length(K) > 0) }
  rows[[length(rows)+1]] <- data.frame(errors=label, slope=s, lrv=m, reps=length(om), both_found=if (length(kn)) mean(ex) else NA,
     missed_window=if (length(kn)) mean(miss) else NA, false_positive=if (!length(kn)) mean(fp) else NA, mean_knots=mean(nk), omega2_hat=mean(om), omega2_true=tr) }
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_dep7_%s.csv", tag), row.names=FALSE)
options(width=170); print(r, row.names=FALSE, digits=3); cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
