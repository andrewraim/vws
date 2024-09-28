#' Equally Spaced Univariate Knots
#'
#' Prepare \eqn{N} equally-spaced intervals between finite endpoints `lo` and
#' `hi`.
#'
#' @param N TBD
#' @param lo TBD
#' @param hi TBD
#'
#' @return TBD
#'
#' @examples
#' print("TBD")
#'
#' @export
equal_spaced_knots = function(N, lo, hi)
{
	stopifnot(N > 0)
	x = seq(lo, hi, length.out = N + 1)
	x |> head(n = -1) |> tail(n = -1)
}
