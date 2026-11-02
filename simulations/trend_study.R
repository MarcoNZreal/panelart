library(panelart)
nrep <- as.integer(Sys.getenv("NREP","60"))
Tn <- 200; N <- 10
scen <- list(
  list(name="null, no trend",            m=0, delta=0,   slope=0,    base=0),
  list(name="null, common trend 0.01",   m=0, delta=0,   slope=0,    base=0.01),
  list(name="slope break 0.02",          m=1, delta=0,   slope=0.02, base=0),
  list(name="slope break 0.04",          m=1, delta=0,   slope=0.04, base=0),
  list(name="slope break 0.03 + trend",  m=1, delta=0,   slope=0.03, base=0.01),
  list(name="level jump 0.6 (no slope)", m=1, delta=0.6, slope=0,    base=0))
meth <- list(classic_bic=function(Y) art_classic(Y,"bic")$dates,
             mean_mc=function(Y) art_fit(Y,"mc")$dates,
             trend_mc=function(Y) art_fit(Y,"mc",type="trend")$dates)
rows <- list(); t0 <- proc.time()[["elapsed"]]
for (si in seq_along(scen)) {
  sc <- scen[[si]]
  for (i in seq_len(nrep)) {
    set.seed(7000 + 97*si + i)
    d <- art_sim_data(Tn, N, breaks = if (sc$m) 0.5 else numeric(0), delta_bar = sc$delta, slope_bar = sc$slope)
    Y <- d$Y
    if (sc$base != 0) Y <- Y + outer(seq_len(Tn), sc$base * sample(c(-1,1), N, TRUE))
    for (nm in names(meth)) {
      est <- tryCatch(meth[[nm]](Y), error=function(e) NULL)
      if (is.null(est)) next
      ev <- art_eval(est, d$dates, Tn)
      rows[[length(rows)+1]] <- data.frame(scenario=sc$name, method=nm, t(ev))
    }
  }
}
r <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(r, list(r$scenario, r$method), drop=TRUE), function(x)
  data.frame(scenario=x$scenario[1], method=x$method[1], n_est=mean(x$n_est), exact=mean(x$exact_count),
             any_break=mean(x$n_est>0), found=if (x$n_true[1]>0) mean(x$found/x$n_true) else NA, mae=mean(x$mae, na.rm=TRUE))))
agg <- agg[order(match(agg$scenario, sapply(scen, `[[`, "name")), agg$method),]
write.csv(agg, "study_trend.csv", row.names=FALSE)
options(width=150); print(data.frame(agg[,1:2], lapply(agg[,3:7], round, 2)), row.names=FALSE)
cat("elapsed", round(proc.time()[["elapsed"]]-t0), "s\n")
