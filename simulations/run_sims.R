library(panelart)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
which <- args[1]
cores <- as.integer(Sys.getenv("CORES", "1"))
nrep  <- as.integer(Sys.getenv("NREP", "200"))
t0 <- proc.time()[["elapsed"]]
wide <- function(res, grid) {
  s <- art_sim_summary(res)
  merge(s, cbind(scenario = seq_len(nrow(grid)), grid), by = "scenario")
}
if (which == "iid") {
  grid <- expand.grid(T = 200, N = c(1, 10, 50), delta_bar = c(0, 0.2, 0.4, 0.8))
  grid$m <- ifelse(grid$delta_bar == 0, 0, 2)
  meth <- art_methods()[c("classic_bic", "panel_iid", "panel_mc")]
  res <- art_sim_study(grid, meth, nrep = nrep, seed = 101, cores = cores)
  write.csv(wide(res, grid), "study_iid.csv", row.names = FALSE)
}
if (which == "ar1") {
  grid <- expand.grid(T = 300, N = 10, delta_bar = c(0, 0.7), phi = c(0, 0.3, 0.6, 0.8))
  grid$errors <- "ar1"; grid$m <- ifelse(grid$delta_bar == 0, 0, 2)
  meth <- art_methods()[c("classic_bic", "panel_mc", "panel_hac")]
  res <- art_sim_study(grid, meth, nrep = nrep, seed = 102, cores = cores)
  write.csv(wide(res, grid), "study_ar1.csv", row.names = FALSE)
}
if (which == "loc") {
  # Localization: fixed pooled jump energy Delta, growing N. Single break.
  Tn <- 300; Delta <- 3
  out <- list()
  for (N in c(3, 10, 30, 100, 300)) {
    dbar <- sqrt(Delta / N)
    err <- mclapply(seq_len(nrep), function(i) {
      set.seed(5000 + N * 7 + i)
      d <- art_sim_data(Tn, N, breaks = 0.5, delta_bar = dbar)
      f <- art_fit(d$Y, "mc")
      if (f$n_breaks == 1) f$dates - d$dates else NA_real_
    }, mc.cores = cores)
    e <- unlist(err)
    out[[length(out) + 1]] <- data.frame(N = N, delta_bar2_T = dbar^2 * Tn, Delta_T = Delta * Tn,
      found = mean(!is.na(e)), med_abs_err = median(abs(e), na.rm = TRUE),
      mean_abs_err = mean(abs(e), na.rm = TRUE), q90_abs_err = quantile(abs(e), 0.9, na.rm = TRUE),
      theory_med = 1.505 / Delta)
  }
  write.csv(do.call(rbind, out), "study_loc.csv", row.names = FALSE)
}
if (which == "ci") {
  out <- list()
  for (cfg in list(c(400, 10, 0.3), c(400, 10, 0.5), c(400, 10, 1.0), c(400, 50, 0.2), c(400, 50, 0.4))) {
    Tn <- cfg[1]; N <- cfg[2]; dbar <- cfg[3]
    res <- mclapply(seq_len(nrep), function(i) {
      set.seed(9000 + i + 17 * N + round(100 * dbar))
      d <- art_sim_data(Tn, N, breaks = 0.5, delta_bar = dbar)
      f <- art_fit(d$Y, "mc")
      if (f$n_breaks == 1) { ci <- art_ci(f); c(hit = as.numeric(d$dates >= ci$lower && d$dates <= ci$upper), w = ci$upper - ci$lower) }
      else c(hit = NA_real_, w = NA_real_)
    }, mc.cores = cores)
    m <- do.call(rbind, res)
    out[[length(out) + 1]] <- data.frame(T = Tn, N = N, delta_bar = dbar, Delta = N * dbar^2,
      found_one = mean(!is.na(m[, "hit"])), coverage95 = mean(m[, "hit"], na.rm = TRUE),
      mean_width = mean(m[, "w"], na.rm = TRUE))
  }
  write.csv(do.call(rbind, out), "study_ci.csv", row.names = FALSE)
}
if (which == "bw") {
  out <- list()
  for (L in c(2, 5, 10, 20, 40, 80)) {
    f <- function(Y) art_fit(Y, "hac_mc", L = L)$dates
    grid <- data.frame(T = 300, N = 10, m = 0, delta_bar = 0, errors = "ar1", phi = 0.5)
    grid2 <- data.frame(T = 300, N = 10, m = 2, delta_bar = 0.7, errors = "ar1", phi = 0.5)
    r0 <- art_sim_summary(art_sim_study(grid, list(hac = f), nrep = nrep, seed = 103, cores = cores))
    r1 <- art_sim_summary(art_sim_study(grid2, list(hac = f), nrep = nrep, seed = 104, cores = cores))
    out[[length(out) + 1]] <- data.frame(L = L, null_P_any_break = r0$p_any_spurious, null_n_est = r0$n_est,
      alt_exact_count = r1$p_exact_count, alt_found = r1$frac_found)
  }
  write.csv(do.call(rbind, out), "study_bw.csv", row.names = FALSE)
}
cat(which, "done in", round(proc.time()[["elapsed"]] - t0), "seconds\n")
