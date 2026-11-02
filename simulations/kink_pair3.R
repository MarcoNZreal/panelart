library(panelart)
nrep <- as.integer(Sys.getenv("NREP","100")); t0 <- proc.time()[["elapsed"]]
Tn <- 300; N <- 10; tt <- 1:Tn
meth <- list(alpha=list(), pair05=list(pair_test=TRUE, pair_level=0.05), pair01=list(pair_test=TRUE, pair_level=0.01))
rows <- list()
for (d in c(40, 25, 15)) for (s in c(0.04, 0.08)) {
  kn <- c(100, 100 + d)
  cnt <- matrix(NA, nrep, length(meth))
  for (i in seq_len(nrep)) {
    set.seed(65000 + i)
    sg <- sample(c(-1,1), N, TRUE)
    Y <- outer(pmax(tt-kn[1],0), s*sg) + outer(pmax(tt-kn[2],0), s*sg) + matrix(rnorm(Tn*N), Tn, N)
    for (v in seq_along(meth)) cnt[i, v] <- length(do.call(art_kink, c(list(Y=Y, penalty="mc", nsim=500), meth[[v]]))$dates)
  }
  for (v in seq_along(meth)) rows[[length(rows)+1]] <- data.frame(spacing=d, slope=s, method=names(meth)[v],
      exact_two=mean(cnt[,v]==2), one_knot=mean(cnt[,v]==1), more=mean(cnt[,v]>2))
  cat("spacing", d, "slope", s, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console()
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_pair3.csv", row.names=FALSE)
options(width=150); print(r, row.names=FALSE, digits=3)
