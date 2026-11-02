library(panelart)
N <- 10; g <- 4; nrep <- as.integer(Sys.getenv("NREP","30"))
out <- list(); t0 <- proc.time()[["elapsed"]]
for (Tn in c(100, 200, 400, 800)) {
  e <- list(trend=numeric(0), kink=numeric(0))
  for (i in seq_len(nrep)) {
    set.seed(11000 + Tn + i)
    d <- art_sim_data(Tn, N, breaks=0.5, delta_bar=0, slope_bar=g/Tn)
    ft <- art_fit(d$Y, "mc", type="trend"); fk <- art_kink(d$Y, "mc")
    if (ft$n_breaks == 1) e$trend <- c(e$trend, abs(ft$dates - d$dates))
    if (fk$n_breaks == 1) e$kink  <- c(e$kink,  abs(fk$dates - d$dates))
  }
  out[[length(out)+1]] <- data.frame(T=Tn, found_trend=length(e$trend)/nrep, found_kink=length(e$kink)/nrep,
    trend_mean_err=mean(e$trend), kink_mean_err=mean(e$kink), trend_med=median(e$trend), kink_med=median(e$kink))
  cat("T", Tn, "done", round(proc.time()[["elapsed"]]-t0), "s\n")
}
r <- do.call(rbind, out); write.csv(r, "study_kink_rate.csv", row.names=FALSE)
options(width=150); print(round(r,2), row.names=FALSE)
sl <- function(y) round(unname(coef(lm(log(y) ~ log(r$T)))[2]), 2)
cat("log-log slope of mean error (dates) vs T: trend", sl(r$trend_mean_err), " kink", sl(r$kink_mean_err), "\n")
cat("=> slope for error as fraction of T: trend", sl(r$trend_mean_err)-1, " kink", sl(r$kink_mean_err)-1, "\n")
