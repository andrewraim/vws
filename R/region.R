#' Region
#'
#' A region contains all of the problem-specific logic for the vws sampler.
#' This is an "abstract" region that shows which functions need to be
#' implemented.
#'
#' @name Region
NULL

#' @name Region
#' @export
Region$methods(d_base = function(x, log = FALSE)
{
	stop("d_base: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(w = function(x, log = TRUE)
{
	stop("w: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(r = function(n) {
	stop("r: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(d = function(x) {
	stop("d: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(in_support = function(x) {
	stop("in_support: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(w_major = function(x, log = TRUE)
{
	stop("w_major: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(bifurcate = function(x = NULL)
{
	stop("bifurcate: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(xi_upper = function(log = TRUE)
{
	stop("xi_upper: please implement for your Region")
})

#' @name Region
#' @export
Region$methods(xi_lower = function(log = TRUE)
{
	stop("xi_lower: please implement for your Region")
})
