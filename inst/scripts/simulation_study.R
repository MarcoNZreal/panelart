## Example Monte Carlo studies with panelart.
## Run from R after installing the package: source(system.file("scripts/simulation_study.R", package = "panelart"))
## Adjust nrep and cores. Results are written to a CSV file in the working directory.
library(panelart)

nrep  <- 200
cores <- 1       # set to parallel::detectCores() - 1 on Linux/macOS

methods <- art_methods()[c("classic_bic", "panel_mc", "panel_hac")]

## --- Study 1: size and power with independent errors ----------------------
## Null (m = 0) and two breaks, over a grid of N and per-series jump sizes.
## The per-series signal that matters in the notes is delta_bar^2 * T / sigma^2.
grid1 <- expand.grid(T = 200, N = c(1, 10, 50), delta_bar = c(0, 0.2, 0.4, 0.8))
grid1$m <- ifelse(grid1$delta_bar == 0, 0, 2)
res1 <- art_sim_study(grid1, methods, nrep = nrep, seed = 1, cores = cores)
sum1 <- art_sim_summary(res1)
sum1 <- merge(sum1, cbind(scenario = seq_len(nrow(grid1)), grid1), by = "scenario")
write.csv(sum1, "study1_iid.csv", row.names = FALSE)

## --- Study 2: serial dependence (AR(1)) -----------------------------------
## Compares the iid-calibrated method with the long-run covariance method.
grid2 <- expand.grid(T = 300, N = 10, delta_bar = c(0, 0.7), phi = c(0, 0.3, 0.6, 0.8))
grid2$errors <- "ar1"
grid2$m <- ifelse(grid2$delta_bar == 0, 0, 2)
res2 <- art_sim_study(grid2, methods, nrep = nrep, seed = 2, cores = cores)
sum2 <- merge(art_sim_summary(res2), cbind(scenario = seq_len(nrow(grid2)), grid2), by = "scenario")
write.csv(sum2, "study2_ar1.csv", row.names = FALSE)

## --- Study 3: strong common factor ----------------------------------------
## Jumps inside versus orthogonal to the factor space.
grid3 <- expand.grid(T = 300, N = 20, delta_bar = 0.7, errors = "factor", r = 1,
                     factor_var = c(1, 4), jump_in = c("factor", "orthogonal"),
                     stringsAsFactors = FALSE)
grid3$m <- 2
methods3 <- art_methods()[c("classic_bic", "panel_mc", "panel_hac", "panel_factor")]
res3 <- art_sim_study(grid3, methods3, nrep = nrep / 4, seed = 3, cores = cores)
sum3 <- merge(art_sim_summary(res3), cbind(scenario = seq_len(nrow(grid3)), grid3), by = "scenario")
write.csv(sum3, "study3_factor.csv", row.names = FALSE)

## --- Study 4: coverage of the limit-law intervals --------------------------
## Single break, shrinking-jump regime.
cover <- function(T, N, delta_bar, nrep) {
  hit <- rep(NA, nrep)
  for (i in seq_len(nrep)) {
    d <- art_sim_data(T, N, breaks = 0.5, delta_bar = delta_bar)
    f <- art_fit(d$Y, "mc")
    if (f$n_breaks == 1) {
      ci <- art_ci(f)
      hit[i] <- (d$dates >= ci$lower) && (d$dates <= ci$upper)
    }
  }
  c(coverage = mean(hit, na.rm = TRUE), found = mean(!is.na(hit)))
}
print(cover(T = 400, N = 10, delta_bar = 0.4, nrep = 500))

## --- Study 5: the bandwidth window ----------------------------------------
## Fixed L from very short to long, AR(1) errors: look for false positives when L is long.
for (L in c(2, 5, 10, 20, 40)) {
  f <- function(Y) art_fit(Y, "hac_mc", L = L)$dates
  sc <- data.frame(T = 300, N = 10, m = 0, delta_bar = 0, errors = "ar1", phi = 0.5)
  r <- art_sim_summary(art_sim_study(sc, list(hac = f), nrep = nrep, seed = 4))
  cat("L =", L, " P(any false break) =", r$p_any_spurious, "\n")
}
