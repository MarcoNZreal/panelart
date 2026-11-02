library(panelart)
clr <- function() rm(list = ls(panelart:::.pkg_cache), envir = panelart:::.pkg_cache)
N <- 10; Tn <- 300; gap <- 6L; t0 <- proc.time()[["elapsed"]]
lam_or <- eigen(0.5*(1+.7)/(1-.7)*matrix(1,N,N) + 0.5*diag(N), symmetric=TRUE, only.values=TRUE)$values
cases <- list(unit = rep(1, N), oracle_factor = lam_or)
rows <- list()
for (nm in names(cases)) {
  lam <- cases[[nm]]
  clr(); ref <- art_penalty_kink_mc_w(lam, Tn, 0.05, 20000L, gap, seed = 1L)
  q500 <- sapply(1:20, function(s) { clr(); art_penalty_kink_mc_w(lam, Tn, 0.05, 500L, gap, seed = 100L + s) })
  q2000 <- sapply(1:12, function(s) { clr(); art_penalty_kink_mc_w(lam, Tn, 0.05, 2000L, gap, seed = 200L + s) })
  rows[[length(rows)+1]] <- data.frame(eigenvalues = nm, reference_20000 = ref, nsim500_mean_ratio = mean(q500)/ref, nsim500_sd_ratio = sd(q500)/ref, nsim500_min = min(q500)/ref, nsim500_max = max(q500)/ref,
        nsim2000_mean_ratio = mean(q2000)/ref, nsim2000_sd_ratio = sd(q2000)/ref)
  cat(nm, "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console() }
r <- do.call(rbind, rows); write.csv(r, "study_kink_calib3_mc_precision.csv", row.names=FALSE); options(width=170); print(r, row.names=FALSE, digits=3)
