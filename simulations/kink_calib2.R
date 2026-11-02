# Exact finite-T variance ratios v'Gamma v / omega^2 for unit residualized step and hinge vectors, AR(1) errors, versus the observed quantile ratios
obs <- data.frame(phi = c(0, .3, .5, .7, .8, .9), q300 = c(.982, .947, .927, .882, .839, .748))
rat <- function(Tn, phi) {
  tt <- 1:Tn; b <- round(Tn/2); om <- (1+phi)/(1-phi)
  Gam <- toeplitz(phi^(0:(Tn-1)))
  h <- qr.resid(qr(cbind(1,tt)), pmax(tt-b,0)); h <- h/sqrt(sum(h^2))
  st <- qr.resid(qr(cbind(1,tt,pmax(tt-b,0))), as.numeric(tt > b)); st <- st/sqrt(sum(st^2))
  c(hinge = drop(crossprod(h, Gam %*% h))/om, step = drop(crossprod(st, Gam %*% st))/om) }
r300 <- t(sapply(obs$phi, function(p) rat(300, p)))
out <- data.frame(phi = obs$phi, hinge_ratio = r300[,1], step_ratio = r300[,2], observed_quantile_ratio = obs$q300)
out$deficit_obs_over_step <- (1 - out$observed_quantile_ratio) / pmax(1 - out$step_ratio, 1e-9)
options(width = 130); print(round(out, 3), row.names = FALSE)
write.csv(out, "study_kink_calib2_exact_ratios.csv", row.names = FALSE)
