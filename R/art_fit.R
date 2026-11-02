#' Fit the extended ART estimator to a panel of time series
#'
#' Detects common mean breaks in the columns of \code{Y}. The estimator runs
#' (A) a forward-repartition-eliminate cycle at penalty \code{alpha}, (B) one
#' round of hull-merge thinning, and (C) repartition and elimination with
#' minimum spacing \code{eps * T}.
#'
#' Penalty modes:
#' \describe{
#'   \item{\code{"iid"}}{Independent errors with a common variance. The variance
#'     is \code{sigma2} or its first-difference estimate.}
#'   \item{\code{"hac"}}{Serially and cross-sectionally dependent errors. The
#'     penalty is built from a long-run covariance estimate computed iteratively
#'     from segment residuals: pass 0 uses the raw estimate with the breaks left
#'     in; later passes use the residuals of the current partition.}
#'   \item{\code{"mc"}}{Independent errors; the penalty is the Monte Carlo
#'     calibrated null quantile of \code{\link{art_penalty_mc}}, times the
#'     variance estimate. Much less conservative than \code{"iid"}.}
#'   \item{\code{"hac_mc"}}{As \code{"hac"}, but each penalty is the Monte Carlo
#'     null quantile for data drawn from the current long-run covariance estimate
#'     (\code{\link{art_penalty_mc_cov}}). Less conservative than \code{"hac"}.}
#'   \item{\code{"fixed"}}{User-supplied \code{alpha}.}
#' }
#'
#' Whitening replaces \code{Y} by \code{Y \%*\% Omega^(-1/2)}. With
#' \code{Omega} given, that matrix is used. \code{whiten = "diff"} estimates
#' the covariance by first differences (serially independent errors).
#' \code{whiten = "lrv"} first fits with \code{penalty = "hac"}, estimates the
#' long-run covariance from the residuals, and refits on whitened data.
#' After whitening the errors have identity covariance, so the penalty is
#' \code{N + 4 sqrt(N log T) + C0 log T} for \code{"iid"} and \code{"hac"} and the
#' Monte Carlo quantile of \code{\link{art_penalty_mc}} for \code{"mc"} and
#' \code{"hac_mc"}.
#'
#' @param Y Panel matrix (T x N), or a vector for a single series.
#' @param penalty One of \code{"iid"}, \code{"hac"}, \code{"fixed"}, \code{"mc"}, \code{"hac_mc"}.
#' @param sigma2 Known error variance for \code{penalty = "iid"}.
#' @param alpha Penalty for \code{penalty = "fixed"}.
#' @param whiten One of \code{"none"}, \code{"diff"}, \code{"lrv"}.
#' @param Omega Known N x N (long-run) covariance; if given, the data are
#'   whitened with it and the iid penalty with variance 1 is used.
#' @param L Bartlett bandwidth for the long-run covariance. Default for the
#'   \code{"hac"} and \code{"hac_mc"} penalties is \code{"auto"} (Andrews'
#'   AR(1) plug-in) from the second pass on; the first pass, which leaves the
#'   breaks in the data, uses \code{floor(T^(1/3))}. See \code{\link{art_lrcov}}.
#' @param eps Minimum spacing fraction of T between dates (phases B and C).
#' @param min_len Minimum distance of a split from the ends of its node in the
#'   forward cycle. Default 2 for \code{"iid"}, \code{max(2, ceiling(sqrt(N)))}
#'   for \code{"hac"}.
#' @param C0 Constant in the penalty; see \code{\link{art_penalty}}.
#' @param max_pass Maximum number of covariance passes for \code{"hac"}.
#' @param tol Stop the passes when the trace of the covariance estimate changes
#'   by less than \code{tol * fro * sqrt(log T)}.
#' @param thin Logical; apply hull-merge thinning.
#' @param level,nsim Null rejection level and number of simulations for
#'   \code{penalty = "mc"}.
#' @param nsim_cov Number of simulations per covariance pass for
#'   \code{penalty = "hac_mc"}.
#' @param verbose Print progress.
#' @param type \code{"mean"} (default): abrupt shifts in the mean, as in classic ART.
#'   \code{"trend"}: breaks in the intercept and/or slope; each segment gets its
#'   own line in every series (discontinuous piecewise-linear fit). The criterion
#'   then has two parameters per series and segment, so null values are about
#'   twice as large; the penalties account for this. The theory in the notes
#'   covers only \code{type = "mean"}; the trend version is an extension that
#'   has been checked by simulation only. Confidence intervals are not available.
#' @return An object of class \code{"art"}: a list with \code{dates},
#'   \code{n_breaks}, \code{segments}, \code{means}, \code{alpha},
#'   \code{penalty}, \code{sigma2}, \code{omega}, \code{history}, the data
#'   \code{Y} and the data used after whitening \code{Yw}.
#' @examples
#' set.seed(1)
#' sim <- art_sim_data(T = 200, N = 10, breaks = c(0.4, 0.7), delta_bar = 0.6)
#' fit <- art_fit(sim$Y)
#' fit
#' @export
art_fit <- function(Y, penalty = c("iid", "hac", "fixed", "mc", "hac_mc"), sigma2 = NULL,
                    alpha = NULL, whiten = c("none", "diff", "lrv"),
                    Omega = NULL, L = NULL, eps = 0.02, min_len = NULL,
                    C0 = 10, max_pass = 10L, tol = 0.125, thin = TRUE,
                    level = 0.05, nsim = 1000L, nsim_cov = 100L, verbose = FALSE,
                    type = c("mean", "trend")) {
  cl <- match.call()
  penalty <- match.arg(penalty); whiten <- match.arg(whiten)
  type <- match.arg(type); df <- if (type == "trend") 2 else 1
  dord <- if (type == "trend") 2L else 1L
  Y0 <- .as_panel(Y); Tn <- nrow(Y0); N <- ncol(Y0)
  if (Tn < 20L) stop("need at least 20 dates")
  lt <- log(Tn)
  Yw <- Y0; Wm <- NULL; whitened <- FALSE
  hac_pen <- if (penalty == "hac") "hac" else "hac_mc"
  if (!is.null(Omega)) {
    Wm <- art_whiten_matrix(Omega); Yw <- Y0 %*% Wm; whitened <- TRUE
  } else if (whiten == "diff") {
    Wm <- art_whiten_matrix(art_sigma_diff(Y0, order = dord), bias_n = Tn)
    Yw <- Y0 %*% Wm; whitened <- TRUE
  } else if (whiten == "lrv") {
    f0 <- art_fit(Y0, penalty = hac_pen, L = L, eps = eps, min_len = min_len,
                  C0 = C0, max_pass = max_pass, tol = tol, thin = thin,
                  level = level, nsim_cov = nsim_cov, type = type)
    om <- art_lrcov(Y0, f0$dates, L, full = TRUE, type = type)
    Wm <- art_whiten_matrix(om$Omega); Yw <- Y0 %*% Wm; whitened <- TRUE
  }
  if (whitened && penalty != "fixed") {
    penalty <- if (penalty %in% c("mc", "hac_mc")) "mc" else "iid"
    sigma2 <- 1
  }
  is_hac <- penalty %in% c("hac", "hac_mc")
  if (is_hac && is.null(L)) L <- "auto"
  min_floor <- if (type == "trend") 3L else 2L
  if (is.null(min_len)) min_len <- if (is_hac) max(min_floor, as.integer(ceiling(sqrt(N)))) else min_floor
  pen_fun <- function(om) {
    if (penalty == "hac_mc")
      art_penalty_mc_cov(if (!is.null(om$Omega)) om$Omega else om$diag, Tn, level, nsim_cov, min_len, type = type)
    else art_penalty(df * om$tr, sqrt(df) * om$fro, om$op, Tn, C0, hat = TRUE)
  }
  S0 <- .cumsum0(Yw, type)
  history <- NULL; omega <- NULL; start <- integer(0)

  if (penalty == "fixed") {
    if (is.null(alpha)) stop("alpha must be given when penalty = 'fixed'")
  } else if (penalty == "iid") {
    if (is.null(sigma2)) sigma2 <- art_sigma2_diff(Yw, dord)
    alpha <- art_penalty_iid(N, Tn, sigma2, C0, df)
  } else if (penalty == "mc") {
    if (is.null(sigma2)) sigma2 <- art_sigma2_diff(Yw, dord)
    alpha <- art_penalty_mc(N, Tn, sigma2, level, nsim, min_len, type)
  } else {
    # pass 0 leaves the breaks in the data, so an automatic bandwidth would read
    # the steps as strong autocorrelation and pick a long window; use a short one
    L0 <- if (identical(L, "auto")) max(2L, as.integer(floor(Tn^(1 / 3)))) else L
    om <- art_lrcov(Yw, integer(0), L0, type = type)
    alpha <- pen_fun(om)
    history <- data.frame(pass = 0L, alpha = alpha, tr = om$tr, n_dates = 0L)
    dates <- integer(0)
    for (pass in seq_len(max_pass)) {
      dates <- art_forward(S0, alpha, dates, min_len)
      om_new <- art_lrcov(Yw, dates, L, type = type)
      alpha_new <- pen_fun(om_new)
      history <- rbind(history, data.frame(pass = pass, alpha = alpha_new,
                                           tr = om_new$tr, n_dates = length(dates)))
      if (verbose) message(sprintf("pass %d: %d dates, alpha %.2f -> %.2f",
                                   pass, length(dates), alpha, alpha_new))
      done <- abs(om_new$tr - om$tr) <= tol * om_new$fro * sqrt(lt)
      om <- om_new; alpha <- alpha_new
      if (done) break
    }
    start <- dates
  }
  dates <- art_forward(S0, alpha, start, min_len)
  if (thin && length(dates) > 1L) dates <- art_thin(S0, dates, eps)
  dates <- .phaseC(S0, dates, alpha, eps, min_len)
  if (is_hac) omega <- art_lrcov(Yw, dates, L, type = type)
  b <- c(0L, dates, Tn)
  structure(list(
    dates = as.integer(dates), n_breaks = length(dates),
    segments = data.frame(start = b[-length(b)] + 1L, end = b[-1L]),
    means = .segment_means(Y0, dates), type = type,
    slopes = if (type == "trend") .segment_slopes(Y0, dates) else NULL,
    alpha = alpha, penalty = penalty, sigma2 = sigma2, omega = omega,
    history = history, whitened = whitened, W = Wm, Y = Y0, Yw = Yw,
    eps = eps, min_len = min_len, C0 = C0, L = L, call = cl),
    class = "art")
}

#' Number of factors by the Bai-Ng criterion
#'
#' Information criterion IC_p2 of Bai and Ng (2002) applied to the principal
#' components of the centred panel.
#'
#' @param Y Panel matrix (T x N).
#' @param rmax Largest number of factors considered.
#' @return Selected number of factors (possibly zero).
#' @export
art_nfactors <- function(Y, rmax = 8L) {
  Y <- .as_panel(Y); Tn <- nrow(Y); N <- ncol(Y)
  rmax <- min(rmax, N - 1L, Tn - 1L)
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  d <- svd(Yc, nu = 0L, nv = 0L)$d
  tot <- sum(Yc^2)
  V <- c(tot, tot - cumsum(d[seq_len(rmax)]^2)) / (N * Tn)
  pen <- (N + Tn) / (N * Tn) * log(min(N, Tn))
  ic <- log(pmax(V, 1e-300)) + (0:rmax) * pen
  which.min(ic) - 1L
}

#' ART with factor projection
#'
#' Two pipelines. The factor pipeline fits \code{art_fit} on the \code{rbar}
#' estimated factor series \code{Y \%*\% Lambda / N}, whitened with a long-run
#' covariance estimate. The orthogonal pipeline fits \code{art_fit} with
#' \code{penalty = "hac"} on the data projected orthogonally to the estimated
#' loadings. The two sets of dates are merged: dates closer than
#' \code{3 * eps * T} are linked and the factor-pipeline date is kept.
#'
#' Loadings are estimated by principal components of the centred panel with
#' \code{rbar = r + extra} components. This implementation does not model the
#' dependence between the loading estimate and the data; see the notes.
#'
#' @param Y Panel matrix (T x N).
#' @param r Number of factors; default from \code{\link{art_nfactors}}.
#' @param rmax Largest number considered when estimating \code{r}.
#' @param extra Number of extra components added to \code{r} (over-selection
#'   is cheap, under-selection is not).
#' @param whiten_factor Whitening of the factor series: \code{"lrv"} (long-run
#'   covariance, as in the notes) or \code{"none"}. Whitening adds one noise
#'   coordinate per extra component, which costs power when \code{extra > 0}.
#' @param ... Passed to \code{\link{art_fit}} for both pipelines (default
#'   \code{penalty = "hac_mc"}).
#' @return A list of class \code{"art_factor"} with \code{dates},
#'   \code{factor_fit}, \code{orth_fit}, \code{r} and \code{rbar}.
#' @export
art_fit_factor <- function(Y, r = NULL, rmax = 8L, extra = 1L,
                           whiten_factor = c("lrv", "none"), ...) {
  whiten_factor <- match.arg(whiten_factor)
  Y <- .as_panel(Y); Tn <- nrow(Y); N <- ncol(Y)
  if (is.null(r)) r <- art_nfactors(Y, rmax)
  rbar <- min(r + extra, N - 1L, Tn - 1L)
  args <- list(...)
  eps <- if (!is.null(args$eps)) args$eps else 0.02
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  V <- svd(Yc, nu = 0L, nv = rbar)$v[, seq_len(rbar), drop = FALSE]
  Fz <- Y %*% V / sqrt(N)
  Yperp <- Y - Y %*% V %*% t(V)
  pen <- if (!is.null(args$penalty)) args$penalty else "hac_mc"
  args$penalty <- NULL
  ff <- do.call(art_fit, c(list(Fz, penalty = pen, whiten = whiten_factor), args))
  fo <- do.call(art_fit, c(list(Yperp, penalty = pen), args))
  d <- sort(c(ff$dates, fo$dates))
  if (length(d) > 1L) {
    gap <- ceiling(3 * eps * Tn)
    comp <- cumsum(c(TRUE, diff(d) >= gap))
    keep <- integer(0)
    for (cc in unique(comp)) {
      idx <- which(comp == cc)
      fd <- intersect(d[idx], ff$dates)
      keep <- c(keep, if (length(fd)) fd[1L] else d[idx][1L])
    }
    d <- sort(keep)
  }
  structure(list(dates = as.integer(d), factor_fit = ff, orth_fit = fo,
                 r = r, rbar = rbar, Y = Y), class = "art_factor")
}
