source("lrv_ar_proto.R"); library(panelart)
Tn <- 300; N <- 10; tt <- 1:Tn; nrep <- 100
sim <- function(seed, ar = numeric(0), ma = numeric(0), s = 0.06, kn = c(100, 200)) {
  set.seed(seed); g <- s * sample(c(-1,1), N, TRUE)
  mu <- Reduce(`+`, lapply(kn, function(b) outer(pmax(tt-b,0), g)))
  E <- sapply(1:N, function(j) as.numeric(arima.sim(list(ar = ar, ma = ma), n = Tn, n.start = 200)))
  list(Y = mu + E, om2 = if (length(ar) || length(ma)) {   # true long-run variance relative to marginal variance 1 (we rescale below)
    1 } else 1) }
truth <- function(ar, ma) { th <- 1 + sum(ma); ph <- 1 - sum(ar); (th/ph)^2 }   # omega^2 for unit innovation variance
rows <- list(); t0 <- proc.time()[["elapsed"]]
cases <- list(list("white", numeric(0), numeric(0)), list("AR1 .3", .3, numeric(0)), list("AR1 .5", .5, numeric(0)), list("AR1 .7", .7, numeric(0)),
              list("AR1 .8", .8, numeric(0)), list("AR1 .9", .9, numeric(0)), list("AR2 (.6,.2)", c(.6,.2), numeric(0)),
              list("ARMA11 (.5,.5)", .5, .5), list("MA1 (.6)", numeric(0), .6))
for (cs in cases) {
  ar <- cs[[2]]; ma <- cs[[3]]; tr <- truth(ar, ma)
  est <- matrix(NA, nrep, 4); old <- numeric(nrep)
  for (i in 1:nrep) { Y <- sim(31000+i, ar, ma)$Y
    est[i,1] <- lrv_ar(Y, 1L, d=1L)$omega2; est[i,2] <- lrv_ar(Y, 2L, d=1L)$omega2
    est[i,3] <- lrv_ar(Y, 1L, d=2L)$omega2; est[i,4] <- median(c(est[i,1]))
    old[i] <- panelart:::.kink_lrv(Y, integer(0), floor(Tn^(1/3)))$omega2 }
  rows[[length(rows)+1]] <- data.frame(errors = cs[[1]], true_omega2 = tr, d1_ar1_mean = mean(est[,1]), d1_ar1_median = median(est[,1]), d1_ar1_sd = sd(est[,1]),
       d1_ar2_median = median(est[,2]), d2_ar1_median = median(est[,3]), old_pass0_mean = mean(old))
  cat(cs[[1]], "done", round(proc.time()[["elapsed"]]-t0), "s\n"); flush.console() }
r <- do.call(rbind, rows); write.csv(r, "study_lrv_ar_accuracy.csv", row.names = FALSE); options(width = 150); print(r, row.names = FALSE, digits = 3)
