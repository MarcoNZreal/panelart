library(panelart)
phi <- as.numeric(Sys.getenv("PHI","0.7")); s <- as.numeric(Sys.getenv("S","0.06")); nrep <- as.integer(Sys.getenv("NREP","40"))
ex <- as.integer(strsplit(Sys.getenv("EXTRA","0,2,4"), ",")[[1]]); mode <- Sys.getenv("MODE","alt")
Tn <- 300; N <- 10; tt <- 1:Tn; kn <- if (mode == "alt") c(100,200) else numeric(0); t0 <- proc.time()[["elapsed"]]
gen <- function(seed) { set.seed(seed); g <- s*sample(c(-1,1),N,TRUE)
  mu <- if (length(kn)) Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g))) else matrix(0,Tn,N)
  E <- apply(matrix(rnorm(Tn*N),Tn,N), 2, function(z) as.numeric(stats::filter(z, phi, "recursive", init=rep(0,1))))
  mu + E*sqrt(1-phi^2) }
rows <- list()
for (e in ex) {
  exact <- miss <- fp <- om <- numeric(0)
  for (i in seq_len(nrep)) {
    if (proc.time()[["elapsed"]] - t0 > 250) break
    f <- art_kink(gen(77000+i), "hac_mc", nsim=500, lrv_extra=e, max_knots=10); K <- f$dates; om <- c(om, f$omega2)
    if (mode == "alt") { near <- sapply(kn, function(b) sum(abs(K-b) <= 30)); exact <- c(exact, length(K)==2 && all(near==1)); miss <- c(miss, any(near==0)) }
    else fp <- c(fp, length(K) > 0)
  }
  rows[[length(rows)+1]] <- data.frame(phi=phi, slope=if (mode=="alt") s else 0, lrv_extra=e, reps=length(om),
     exact=if (mode=="alt") mean(exact) else NA, missed_window=if (mode=="alt") mean(miss) else NA,
     false_positive=if (mode=="null") mean(fp) else NA, omega2_hat=mean(om), omega2_true=(1+phi)/(1-phi))
}
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_dep6_%s_phi%s_s%s.csv", mode, phi, s), row.names=FALSE)
options(width=150); print(r, row.names=FALSE, digits=3); cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
