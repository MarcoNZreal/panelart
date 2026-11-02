library(panelart)
nrep <- as.integer(Sys.getenv("NREP","100")); t0 <- proc.time()[["elapsed"]]
cfgs <- list(list(name="2 kinks, T=300, s=0.04", T=300, kn=c(100,200), s=0.04),
             list(name="2 kinks, T=300, s=0.06", T=300, kn=c(100,200), s=0.06),
             list(name="3 kinks, T=400, s=0.05", T=400, kn=c(100,200,300), s=0.05))
rows <- list()
for (cf in cfgs) for (pen in c("mc","iid")) {
  Tn <- cf$T; N <- 10; tt <- 1:Tn; m <- length(cf$kn)
  out <- matrix(NA, nrep, 5)
  for (i in seq_len(nrep)) {
    set.seed(43000 + i)
    g <- cf$s * sample(c(-1,1), N, TRUE)
    mu <- Reduce(`+`, lapply(cf$kn, function(b) outer(pmax(tt-b,0), g)))
    Y <- mu + matrix(rnorm(Tn*N), Tn, N)
    f <- art_kink(Y, pen, nsim=500)
    d <- f$dates
    near <- sapply(d, function(k) any(abs(k - cf$kn) <= 30))
    perwin <- sapply(cf$kn, function(b) sum(abs(d - b) <= 30))
    out[i,] <- c(length(d)==m, length(d)>m, any(!near), any(perwin>=2), any(perwin==0))
  }
  rows[[length(rows)+1]] <- data.frame(design=cf$name, penalty=pen, alpha=round(f$alpha,1), reps=nrep,
    exact_count=mean(out[,1]), extras=mean(out[,2]), far_knot=mean(out[,3]), duplicate=mean(out[,4]), window_missed=mean(out[,5]))
  cat(cf$name, pen, "done", round(proc.time()[["elapsed"]]-t0), "s\n")
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_extras.csv", row.names=FALSE)
options(width=170); print(r, row.names=FALSE, digits=3)
