library(panelart)
nrep <- as.integer(Sys.getenv("NREP","120")); t0 <- proc.time()[["elapsed"]]
cfgs <- list(list(name="2 kinks, T=300, s=0.025", T=300, kn=c(100,200), s=0.025),
             list(name="2 kinks, T=300, s=0.03",  T=300, kn=c(100,200), s=0.03),
             list(name="2 kinks, T=300, s=0.04",  T=300, kn=c(100,200), s=0.04),
             list(name="3 kinks, T=400, s=0.05",  T=400, kn=c(100,200,300), s=0.05),
             list(name="null, T=300",             T=300, kn=numeric(0), s=0))
meth <- list(alpha=list(), twoN=list(alpha_thin="2N"),
             pair05=list(pair_test=TRUE, pair_level=0.05), pair01=list(pair_test=TRUE, pair_level=0.01))
sel <- as.integer(strsplit(Sys.getenv("CFG","1,2,3,4,5"), ",")[[1]]); tag <- Sys.getenv("TAG","all")
rows <- list()
for (cf in cfgs[sel]) {
  Tn <- cf$T; N <- 10; tt <- 1:Tn; m <- length(cf$kn)
  res <- array(NA, c(nrep, length(meth), 3))
  for (i in seq_len(nrep)) {
    set.seed(51000 + i)
    g <- cf$s * sample(c(-1,1), N, TRUE)
    mu <- if (m) Reduce(`+`, lapply(cf$kn, function(b) outer(pmax(tt-b,0), g))) else matrix(0, Tn, N)
    Y <- mu + matrix(rnorm(Tn*N), Tn, N)
    for (v in seq_along(meth)) {
      f <- do.call(art_kink, c(list(Y=Y, penalty="mc", nsim=500), meth[[v]]))
      d <- f$dates
      res[i, v, ] <- c(length(d) == m, length(d) > m, if (m) any(sapply(cf$kn, function(b) all(abs(d-b) > 30))) else FALSE)
    }
  }
  for (v in seq_along(meth)) rows[[length(rows)+1]] <- data.frame(design=cf$name, method=names(meth)[v],
      exact_count=mean(res[,v,1]), extra_runs=mean(res[,v,2]), missed_window=mean(res[,v,3]))
  cat(cf$name, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console()
}
r <- do.call(rbind, rows); write.csv(r, paste0("study_kink_pair1_", tag, ".csv"), row.names=FALSE)
options(width=150); print(r, row.names=FALSE, digits=3)
