library(panelart); library(parallel)
which <- commandArgs(trailingOnly = TRUE)[1]
nrep <- as.integer(Sys.getenv("NREP", "100"))
t0 <- proc.time()[["elapsed"]]

if (which == "scale") {
  # Fixed pooled signal s = Delta * T / sigma^2, growing T and N. Two breaks (null when s = 0).
  grid <- expand.grid(T = c(100, 400, 1000), N = c(5, 50), s = c(0, 300, 600))
  grid$m <- ifelse(grid$s == 0, 0, 2)
  grid$delta_bar <- sqrt(grid$s / (grid$N * grid$T))
  meth <- art_methods()["panel_mc"]
  res <- art_sim_study(grid[, c("T", "N", "m", "delta_bar")], meth, nrep = nrep, seed = 201)
  s <- merge(art_sim_summary(res), cbind(scenario = seq_len(nrow(grid)), grid), by = "scenario")
  s$mae_frac <- s$mae / s$T
  write.csv(s, "study_scale.csv", row.names = FALSE)
}

if (which == "het") {
  # Three breaks with unequal jump sizes; per-break detection by method.
  profiles <- list(moderate = c(1, 0.6, 0.35), strong = c(1, 0.3, 0.15))
  out <- list()
  for (err in c("iid", "ar1")) for (pn in names(profiles)) for (base in c(0.5, 1.0)) {
    Tn <- 300; N <- 20; breaks <- c(0.25, 0.5, 0.75)
    db <- base * profiles[[pn]]
    methods <- list(classic_bic = function(Y) art_classic(Y, "bic")$dates,
                    panel_mc = function(Y) art_fit(Y, "mc")$dates,
                    panel_hac = function(Y) art_fit(Y, "hac_mc")$dates)
    rows <- list()
    for (i in seq_len(nrep)) {
      set.seed(3000 + i + 131 * base + 7 * match(pn, names(profiles)) + 11 * (err == "ar1"))
      d <- art_sim_data(Tn, N, breaks = breaks, delta_bar = db, errors = err, phi = if (err == "ar1") 0.3 else 0)
      tol <- max(5, ceiling(0.02 * Tn))
      for (nm in names(methods)) {
        est <- tryCatch(methods[[nm]](d$Y), error = function(e) NULL)
        if (is.null(est)) next
        fe <- sapply(d$dates, function(b) any(abs(est - b) <= tol))
        sp <- sum(vapply(est, function(e) all(abs(e - d$dates) > tol), logical(1)))
        rows[[length(rows) + 1]] <- data.frame(method = nm, large = fe[1], medium = fe[2], small = fe[3],
                                               n_est = length(est), spurious = sp)
      }
    }
    r <- do.call(rbind, rows)
    agg <- aggregate(cbind(large, medium, small, n_est, spurious > 0) ~ method, r, mean)
    names(agg)[ncol(agg)] <- "p_any_spurious"
    out[[length(out) + 1]] <- cbind(errors = err, profile = pn, base = base,
      sizes = paste(round(db, 2), collapse = "/"), agg)
  }
  write.csv(do.call(rbind, out), "study_het.csv", row.names = FALSE)
}

if (which == "factor") {
  grid <- data.frame(
    T = 300, N = c(20, 20, 50, 50, 20, 20, 20, 50),
    m = c(2, 2, 2, 2, 2, 2, 0, 0),
    delta_bar = c(0.7, 0.7, 0.45, 0.45, 0.7, 0.7, 0, 0),
    errors = "factor", phi = 0, r = c(1, 1, 1, 1, 2, 2, 1, 1), factor_var = 2,
    jump_in = c("factor", "orthogonal", "factor", "orthogonal", "factor", "orthogonal", "random", "random"),
    same_direction = TRUE, stringsAsFactors = FALSE)
  meth <- c(art_methods()[c("classic_bic", "panel_mc", "panel_hac", "panel_factor")],
            list(panel_factor_nowhite = function(Y) art_fit_factor(Y, whiten_factor = "none")$dates))
  res <- art_sim_study(grid, meth, nrep = nrep, seed = 202)
  s <- merge(art_sim_summary(res), cbind(scenario = seq_len(nrow(grid)), grid), by = "scenario")
  write.csv(s, "study_factor.csv", row.names = FALSE)
}
cat(which, "done in", round(proc.time()[["elapsed"]] - t0), "seconds\n")
