library(panelart)
nrep <- as.integer(Sys.getenv("NREP","400")); t0 <- proc.time()[["elapsed"]]; rows <- list()
for (cf in list(c(200,10,0.04), c(200,10,0.02), c(200,50,0.015))) {
  Tn <- cf[1]; N <- cf[2]; sl <- cf[3]; A <- rep(NA_real_, nrep); B <- A; ea <- A; eb <- A
  for (i in seq_len(nrep)) {
    set.seed(33000 + i)
    d <- art_sim_data(Tn, N, breaks=0.5, delta_bar=0, slope_bar=sl)
    f <- art_kink(d$Y, "mc", nsim=200)
    A[i] <- length(f$dates_before_thin)==1; B[i] <- length(f$dates)==1
    if (length(f$dates_before_thin)==1) ea[i] <- abs(f$dates_before_thin - d$dates)
    if (length(f$dates)==1) eb[i] <- abs(f$dates - d$dates)
  }
  rows[[length(rows)+1]] <- data.frame(T=Tn, N=N, slope=sl, reps=nrep, exact_before=mean(A), exact_after=mean(B), mae_before=mean(ea,na.rm=TRUE), mae_after=mean(eb,na.rm=TRUE))
}
r <- do.call(rbind, rows); write.csv(r, "study_kink_thin3.csv", row.names=FALSE); print(round(r,3), row.names=FALSE)
fp <- sapply(1:400, function(i) { set.seed(35000+i); f <- art_kink(matrix(rnorm(200*10),200,10), "mc", nsim=200); c(length(f$dates_before_thin)>0, length(f$dates)>0) })
cat("null panels T=200 N=10, 400 reps: P(any knot) before thinning", mean(fp[1,]), " after", mean(fp[2,]), "\n")
cat("elapsed", round(proc.time()[["elapsed"]]-t0), "s\n")
