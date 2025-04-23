#' Rejection Sampling with VWS Proposal
#'
#' Rejection sampling for univariate targets using a constant majorizer.
#' Constants are obtained using numerical optimization. The proposal is created
#' internally, used, and then discarded. See the C++ framework if access to the
#' proposal is needed or other variations of VWS are desired.
#'
#' @param n Number of desired draws.
#' @param lo Lower bound of support.
#' @param hi Upper bound of support.
#' @param w The weight function \eqn{w}.
#' @param helper A helper object constructed from [univariate_helper] that
#' characterizes the base distribution \eqn{g}.
#' @param control An object constructed using [rejection_control].
#'
#' @details
#' VWS requires decomposition of the target \eqn{f} into
#' \eqn{f(x) \propto w(x) g(x)}.
#'
#' The functions `rejection_numeric` and `rejection_int` are for targets with
#' continuous and integer supports, respectively.
#'
#' @return A list with the following elements.
#' - `draws`: Vector with accepted draws.
#' - `rejects`: Vector with counts of rejections needed to obtain each draw.
#' - `lbdd`: Value of bound for the log-probability of rejection over
#'   refinement of the proposal. The maximum length of this vector is \eqn{N},
#'   but it may be shorter if tolerance `tol` was achieved before partitioning
#'   into \eqn{N} regions.
#'
#' The length of `draws` and `rejects` will be \eqn{n} unless the maximum
#' number of rejections was exceeded.
#'
#' @examples
#' library(tidyverse)
#' set.seed(1234)
#'
#' z = 66
#' mu = 5
#' sigma = 0.25
#' lambda2 = 100
#' tol = 0.2
#'
#' helper = normal_helper(mean = z, sd = sqrt(lambda2))
#' w = function(y, log = TRUE) {
#'   dlnorm(y, mu, sigma, log = log)
#' }
#'
#' out = rejection_numeric(n = 10000, lo = 0, hi = Inf, N = 50, tol = tol, w,
#'   helper, control = rejection_control(max_rejects = 10000, action = "stop"))
#' sum(out$rejects)
#'
#' # Some plots of the results
#' data.frame(x = out$draws) %>%
#'   ggplot() +
#'   geom_histogram(aes(x = x), col = "black", bins = 20) +
#'   theme_light()
#'
#' data.frame(draw = seq_along(out$lbdd), lbdd = out$lbdd) %>%
#'   ggplot() +
#'   geom_line(aes(x = draw, y = lbdd)) +
#'   geom_hline(yintercept = log(tol), lty = 2) +
#'   xlab("Draw") +
#'   ylab("Log of Bound for Rejection Probability") +
#'   theme_light()
#'
#' data.frame(draw = seq_along(out$rejects), crej = cumsum(out$rejects)) %>%
#'   ggplot() +
#'   geom_line(aes(x = draw, y = crej)) +
#'   xlab("Draw") +
#'   ylab("Cumulative Rejection Count") +
#'   theme_light()
#'
#' print("TBD: add an example with integer support")
#'
#' @name rejection
NULL

#' @name rejection
#' @export
rejection_numeric = function(n, lo, hi, w, helper,
	control = rejection_control())
{
	# Require for now that maxopt and minopt are either both specified or both
	# unspecified.
	stopifnot(is.null(control$maxopt) == is.null(control$minopt))

	if (is.null(control$opt)) {
		out = rejection_numeric_rcpp(n = n, lo = lo, hi = hi, w = w,
			d_base = helper$d, p_base = helper$p, q_base = helper$q,
			N = control$N, tol = control$tol, control = control)
	} else {
		out = rejection_numeric_opt_rcpp(n = n, lo = lo, hi = hi, w = w,
			d_base = helper$d, p_base = helper$p, q_base = helper$q,
			maxopt = control$maxopt, minopt = control$minopt,
			N = control$N, tol = control$tol, control = control)
	}

	return(out)
}

#' @name rejection
#' @export
rejection_int = function(n, lo, hi, w, helper,
	control = rejection_control())
{
	# Require for now that maxopt and minopt are either both specified or both
	# unspecified.
	stopifnot(is.null(control$maxopt) == is.null(control$minopt))

	if (is.null(control$maxopt)) {
		out = rejection_int_rcpp(n = n, lo = lo, hi = hi, w = w,
			d_base = helper$d, p_base = helper$p, q_base = helper$q,
			N = control$N, tol = control$tol, control = control)
	} else {
		out = rejection_int_opt_rcpp(n = n, lo = lo, hi = hi, w = w,
			d_base = helper$d, p_base = helper$p, q_base = helper$q,
			maxopt = control$maxopt, minopt = control$minopt,
			N = control$N, tol = control$tol, control = control)
	}

	return(out)
}
