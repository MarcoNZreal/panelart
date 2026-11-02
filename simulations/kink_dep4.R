library(panelart)
t0 <- proc.time()[["elapsed"]]; Tn <- 300; N <- 10; tt <- 1:Tn
gen <- function(seed, rng, s, kn) { set.seed(seed); phis <- runif(N, rng[1], rng[2]); g <- s*sample(c(-1,1),N,TRUE)
  E <- sapply(seq_len(N), function(j) as.numeric(stats::filter(rnorm(Tn), phis[j], "recursive")) * sqrt(1-phis[j]^2))
  mu <- if (length(kn)) Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g))) else matrix(0,Tn,N); mu + E }
rows <- list()
for (rng in list(c(0,0.4), c(0,0.8), c(0.5,0.9))) {
  nrep <- 200; fp <- matrix(NA, nrep, 2)
  for (i in seq_len(nrep)) { Y <- gen(99000+i, rng, 0, numeric(0))
    fp[i,1] <- length(art_kink(Y, "hac_mc", nsim=500)$dates) > 0
    fp[i,2] <- length(art_kink(Y, "hac_mc", nsim=500, scale="series")$dates) > 0 }
  rows[[length(rows)+1]] <- data.frame(design=sprintf("null, phi~U(%.1f,%.1f)", rng[1], rng[2]), reps=nrep, pooled=mean(fp[,1]), series=mean(fp[,2]))
  cat(rows[[length(rows)]]$design, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console()
}
for (rng in list(c(0,0.8))) {
  nrep <- 50; ex <- matrix(NA, nrep, 2); ms <- matrix(NA, nrep, 2)
  for (i in seq_len(nrep)) { Y <- gen(98000+i, rng, 0.06, c(100,200))
    for (v in 1:2) { K <- art_kink(Y, "hac_mc", nsim=500, scale=c("pooled","series")[v])$dates
      near <- sapply(c(100,200), function(b) sum(abs(K-b) <= 30)); ex[i,v] <- length(K) == 2 && all(near == 1); ms[i,v] <- any(near == 0) } }
  rows[[length(rows)+1]] <- data.frame(design="2 kinks s=0.06, phi~U(0,0.8): exact count", reps=nrep, pooled=mean(ex[,1]), series=mean(ex[,2]))
  rows[[length(rows)+1]] <- data.frame(design="2 kinks s=0.06, phi~U(0,0.8): missed window", reps=nrep, pooled=mean(ms[,1]), series=mean(ms[,2]))
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_dep4.csv", row.names=FALSE); options(width=140); print(r, row.names=FALSE, digits=3)
cat("sec", round(proc.time()[["elapsed"]]-t0), "\n")
