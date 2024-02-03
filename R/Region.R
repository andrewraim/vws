#' Region
#'
#' A region contains all of the problem-specific logic for the vws sampler.
#' This is an "abstract" region that shows which functions need to be
#' implemented.
#'
#' The following methods should be implemented:
#' \describe{
#' \item{\code{d_base}}{Density function \eqn{g} for the base distribution.}
#' \item{\code{w}}{Weight function \eqn{w}.}
#' \item{\code{r}}{Generate a draw from \eqn{g_j} specific to this region.}
#' \item{\code{d}}{Density of \eqn{g_j} specific to this region.}
#' \item{\code{in_support}}{Test if given \code{x} is in the support for the
#' \eqn{g_j} specific to this region.}
#' \item{\code{w_major}}{Majorized weight function \eqn{\overline{w}_j} for
#' this region.}
#' \item{\code{bifurcate}}{	Bifurcate this region into two regions. Use
#' \\code{x} as the bifurcation point if it is not \code{NULL}. Otherwise,
#' select a point for bifurcation.}
#' \item{\code{xi_upper}}{The quantity \eqn{\overline{\xi}_j} for this
#' region.}
#' \item{\code{xi_lower}}{The quantity \eqn{\underline{\xi}_j} for this
#' region.}
#' }
#'
#' @name Region
NULL

Region$methods(d_base = function(x, log = FALSE)
{
	stop("d_base: please implement for your Region")
})

Region$methods(w = function(x, log = TRUE)
{
	stop("w: please implement for your Region")
})

Region$methods(r = function(n) {
	stop("r: please implement for your Region")
})

Region$methods(d = function(x) {
	stop("d: please implement for your Region")
})

Region$methods(in_support = function(x) {
	stop("in_support: please implement for your Region")
})

Region$methods(w_major = function(x, log = TRUE)
{
	stop("w_major: please implement for your Region")
})

Region$methods(bifurcate = function(x = NULL)
{
	stop("bifurcate: please implement for your Region")
})

Region$methods(xi_upper = function(log = TRUE)
{
	stop("xi_upper: please implement for your Region")
})

Region$methods(xi_lower = function(log = TRUE)
{
	stop("xi_lower: please implement for your Region")
})
