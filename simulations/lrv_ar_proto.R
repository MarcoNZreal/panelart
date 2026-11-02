# Prototype: long-run variance from an AR(p) fitted to second differences (robust to piecewise-linear mean)
# AR coefficients and autocorrelations (lags 0..lm) from partial autocorrelations r (Durbin-Levinson)
.pacf2acf <- function(r, lm) {
  p <- length(r); phi <- r[1]; rho <- numeric(max(lm, p) + 1); rho[1] <- 1; rho[2] <- r[1]
  if (p >= 2) for (k in 2:p) {
    pk <- phi
    rho[k + 1] <- r[k] * (1 - sum(pk * rho[2:k])) + sum(pk * rho[k:2])
    phi <- c(pk - r[k] * rev(pk), r[k]) }
  if (max(lm, p) > p) for (l in (p + 1):max(lm, p)) rho[l + 1] <- sum(phi * rho[l:(l - p + 1)])
  list(phi = phi, rho = rho)
}
.acov_d2 <- function(r, s2u, H, d = 2L) {
  p <- length(r); lm <- H + 2; a <- .pacf2acf(r, lm); phi <- a$phi; rho <- a$rho
  g0 <- s2u / (1 - sum(phi * rho[2:(p + 1)])); ge <- g0 * rho; cc <- if (d == 1L) c(1, -1) else c(1, -2, 1); q <- length(cc) - 1L
  vapply(0:H, function(h) { s <- 0; for (j in 0:q) for (k in 0:q) s <- s + cc[j+1]*cc[k+1]*ge[abs(h + j - k) + 1]; s }, numeric(1)) }
lrv_ar <- function(Y, p = 1L, H = 10L, d = 1L) {
  D <- diff(Y, differences = d); D <- sweep(D, 2L, colMeans(D), "-"); n <- nrow(D)
  gh <- vapply(0:H, function(h) mean(colSums(D[1:(n-h), , drop=FALSE] * D[(1+h):n, , drop=FALSE]) / n), numeric(1))
  obj <- function(th) { v <- sum((gh - .acov_d2(tanh(th[1:p]), exp(th[p+1]), H, d))^2); if (is.finite(v)) v else 1e12 }
  best <- NULL
  for (st in list(c(atanh(0.3), rep(0, p-1)), c(atanh(0.8), rep(0, p-1)))) {
    o <- try(optim(c(st, log(gh[1] / (if (d == 1L) 2 else 6))), obj, method = "BFGS", control = list(maxit = 200)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o }
  phi <- .pacf2acf(pmax(pmin(tanh(best$par[1:p]), 0.97), -0.97), 1L)$phi; s2u <- exp(best$par[p+1])
  list(omega2 = s2u / (1 - sum(phi))^2, phi = phi, s2u = s2u) }
