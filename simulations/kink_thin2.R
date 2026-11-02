library(panelart)
nrep <- as.integer(Sys.getenv("NREP","120")); Tn <- 300; N <- 10; tt <- 1:Tn; t0 <- proc.time()[["elapsed"]]
dupl <- function(d, truth=c(100,200), r=30) any(sapply(truth, function(b) sum(abs(d-b)<=r) >= 2))
far  <- function(d, truth=c(100,200), r=30) any(sapply(d, function(k) all(abs(k-truth)>r)))
rows <- list()
for (sg in c(1,-1)) for (s in c(0.03, 0.04)) {
  A <- matrix(NA_real_, nrep, 3); B <- A
  for (i in seq_len(nrep)) {
    set.seed(31000 + i)
    g1 <- s * sample(c(-1,1), N, TRUE); g2 <- sg * g1
    Y <- outer(pmax(tt-100,0), g1) + outer(pmax(tt-200,0), g2) + matrix(rnorm(Tn*N), Tn, N)
    f <- art_kink(Y, "mc", nsim=200)
    A[i,] <- c(length(f$dates_before_thin)==2, dupl(f$dates_before_thin), far(f$dates_before_thin))
    B[i,] <- c(length(f$dates)==2, dupl(f$dates), far(f$dates))
  }
  rows[[length(rows)+1]] <- data.frame(pattern=if (sg==1) "same direction" else "tent", slope=s, reps=nrep,
    exact_before=mean(A[,1]), exact_after=mean(B[,1]), dup_runs_before=mean(A[,2]), dup_runs_after=mean(B[,2]),
    far_knot_runs_before=mean(A[,3]), far_knot_runs_after=mean(B[,3]))
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_thin2.csv", row.names=FALSE)
options(width=170); print(data.frame(r[,1:3], lapply(r[,-(1:3)], round, 3)), row.names=FALSE)
cat("elapsed", round(proc.time()[["elapsed"]]-t0), "s\n")
