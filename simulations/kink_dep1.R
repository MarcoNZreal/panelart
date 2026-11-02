library(panelart)
phi <- as.numeric(Sys.getenv("PHI","0.5")); nrep <- as.integer(Sys.getenv("NREP","40")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn
gen <- function(seed, s, kn) { set.seed(seed); g <- s*sample(c(-1,1),N,TRUE)
  mu <- if (length(kn)) Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g))) else matrix(0,Tn,N)
  E <- apply(matrix(rnorm(Tn*N),Tn,N), 2, function(z) as.numeric(stats::filter(z, phi, "recursive", init=rep(0,1)))) 
  mu + E*sqrt(1-phi^2) }  # stationary-variance 1 scale: innovation sd sqrt(1-phi^2) so that Var(e)=1 (approximately after burn-in)
## NOTE: simulated with innovation sd sqrt(1-phi^2), so long-run variance omega^2 = (1-phi^2)/(1-phi)^2 = (1+phi)/(1-phi)
om2_true <- (1+phi)/(1-phi)
kinds <- function(K, kn) { m <- length(kn)
  near <- if (m) sapply(kn, function(b) sum(abs(K-b) <= 30)) else integer(0)
  c(exact = length(K) == m && (m == 0 || all(near == 1)), extra = length(K) > m, missed = m > 0 && any(near == 0), nk = length(K)) }
designs <- list(list(name="null", s=0, kn=numeric(0)), list(name="2 kinks s=0.04", s=0.04, kn=c(100,200)),
                list(name="2 kinks s=0.06", s=0.06, kn=c(100,200)))
sel <- strsplit(Sys.getenv("DES","null,alt04,alt06"), ",")[[1]]; tag <- Sys.getenv("TAG","a")
designs <- designs[c(null=1, alt04=2, alt06=3)[sel]]
rows <- list()
for (d in designs) {
  R <- list(mc=matrix(NA,nrep,4), hac=matrix(NA,nrep,4)); om <- numeric(nrep)
  for (i in seq_len(nrep)) {
    Y <- gen(77000+i, d$s, d$kn)
    R$mc[i,]  <- kinds(art_kink(Y, "mc", nsim=500, max_knots=10)$dates, d$kn)
    f <- art_kink(Y, "hac_mc", nsim=500); R$hac[i,] <- kinds(f$dates, d$kn); om[i] <- f$omega2
  }
  for (m in names(R)) rows[[length(rows)+1]] <- data.frame(phi=phi, design=d$name, method=if (m=="mc") "mc (iid)" else "hac_mc",
     reps=nrep, exact=mean(R[[m]][,1]), extra_runs=mean(R[[m]][,2]), missed_window=mean(R[[m]][,3]), mean_knots=mean(R[[m]][,4]),
     omega2_hat=if (m=="hac") mean(om) else NA, omega2_true=om2_true)
  cat(d$name, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console()
}
r <- do.call(rbind, rows); write.csv(r, sprintf("study_kink_dep1_phi%s_%s.csv", phi, tag), row.names=FALSE)
options(width=170); print(r, row.names=FALSE, digits=3)
