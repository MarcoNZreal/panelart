source("lrv_ar_proto.R"); library(panelart)
Tn <- 300; N <- 10; tt <- 1:Tn; nrep <- 80; t0 <- proc.time()[["elapsed"]]
truth <- function(ar, ma) ((1 + sum(ma)) / (1 - sum(ar)))^2
cases <- list(list("white", numeric(0), numeric(0)), list("AR1 .5", .5, numeric(0)), list("AR1 .8", .8, numeric(0)),
              list("AR2 (.6,.2)", c(.6,.2), numeric(0)), list("AR2 (1.2,-.4)", c(1.2,-.4), numeric(0)),
              list("ARMA11 (.5,.5)", .5, .5), list("ARMA11 (.8,-.4)", .8, -.4), list("MA1 (.6)", numeric(0), .6), list("MA1 (-.5)", numeric(0), -.5))
rows <- list()
for (cs in cases) {
  ar <- cs[[2]]; ma <- cs[[3]]; tr <- truth(ar, ma); est <- matrix(NA, nrep, 4)
  for (i in 1:nrep) { set.seed(31000+i); g <- 0.06 * sample(c(-1,1), N, TRUE)
    Y <- outer(pmax(tt-100,0), g) + outer(pmax(tt-200,0), g) + sapply(1:N, function(j) as.numeric(arima.sim(list(ar = ar, ma = ma), n = Tn, n.start = 200)))
    for (p in 1:4) est[i,p] <- lrv_ar(Y, p, d = 1L)$omega2 }
  q <- function(x) quantile(x / tr, c(.1, .5, .9))
  rows[[length(rows)+1]] <- data.frame(errors = cs[[1]], true = tr, do.call(cbind, lapply(1:4, function(p) { z <- q(est[,p]); setNames(round(z,2), paste0("p", p, c("_q10","_med","_q90"))) }) |> (\(m) as.data.frame(t(unlist(as.list(m)))))()))
  cat(cs[[1]], "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console() }
r <- do.call(rbind, rows); write.csv(r, "study_lrv_ar_orders.csv", row.names = FALSE); options(width = 220); print(r, row.names = FALSE, digits = 3)
