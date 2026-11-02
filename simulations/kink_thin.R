library(panelart)
nrep <- as.integer(Sys.getenv("NREP","60")); Tn <- 300; N <- 10; tt <- 1:Tn; t0 <- proc.time()[["elapsed"]]
rows <- list()
for (sg in c(1,-1)) for (s in c(0.02, 0.03, 0.04)) {
  a <- matrix(NA_real_, nrep, 3); b <- a
  for (i in seq_len(nrep)) {
    set.seed(23000 + i)
    g1 <- s * sample(c(-1,1), N, TRUE); g2 <- sg * g1
    Y <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g2) + matrix(rnorm(Tn*N), Tn, N)
    f <- art_kink(Y, "mc", nsim=200)
    ev <- function(d) c(length(d), sum(abs(d-100)<=5)>0 && sum(abs(d-200)<=5)>0, mean(c(min(abs(d-100)), min(abs(d-200)))))
    a[i,] <- ev(f$dates_before_thin); b[i,] <- ev(f$dates)
  }
  rows[[length(rows)+1]] <- data.frame(pattern=if (sg==1) "same direction" else "tent", slope=s,
    exact_before=mean(a[,1]==2), exact_after=mean(b[,1]==2), knots_before=mean(a[,1]), knots_after=mean(b[,1]),
    both_found_before=mean(a[,2]), both_found_after=mean(b[,2]), mae_before=mean(a[,3]), mae_after=mean(b[,3]))
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_thin.csv", row.names=FALSE)
options(width=170); print(data.frame(r[,1:2], lapply(r[,-(1:2)], round, 2)), row.names=FALSE)
## null false-positive rate and one-kink designs
fp <- sapply(1:200, function(i) { set.seed(25000+i); length(art_kink(matrix(rnorm(200*10),200,10), "mc", nsim=200)$dates) })
cat("null panels (T=200, N=10): P(any knot) =", mean(fp>0), " (200 reps)\n")
cat("elapsed", round(proc.time()[["elapsed"]]-t0), "s\n")
