#' panelart: common breaks in panels with regression trees
#'
#' The package works on a \eqn{T \times N} matrix \code{Y} (rows are dates,
#' columns are series). A break date \code{d} means the mean changes between
#' row \code{d} and row \code{d + 1}.
#'
#' Main entry points: \code{\link{art_fit}}, \code{\link{art_classic}},
#' \code{\link{art_fit_factor}}, \code{\link{art_kink}} (continuous
#' piecewise-linear breaks), \code{\link{art_sim_data}},
#' \code{\link{art_sim_study}}.
#'
#' Several help pages refer to \dQuote{the notes}: theory notes with proofs,
#' each result carrying a status label, which accompany a paper in preparation
#' and are not part of the package. Statements about simulation behaviour in the
#' help pages come from the simulation studies in the source repository.
#'
#' @importFrom stats rnorm sd integrate uniroot pnorm setNames aggregate runif quantile qnorm
#' @importFrom graphics matplot abline lines par plot legend mtext
#' @importFrom grDevices adjustcolor
#' @importFrom parallel mclapply
#' @keywords internal
"_PACKAGE"
