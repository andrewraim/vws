#' int_univariate_const_region
#'
#' A more friendly constructor for \code{IntUnivariateConstRegion}.
#'
#' @param a Lower knot of region.
#' @param b Upper knot of region.
#' @param w Weight function.
#' @param g Object that encapsulates base distribution.
#'
#' @examples
#' # Define base distribution and weight function
#' g = poisson_univariate_helper(lambda = 5)
#' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
#'
#' reg = int_univariate_const_region(-Inf, 10, w, g)
#' print(reg)
#'
#' out = reg$bifurcate(5)
#' print(out[[1]])
#' print(out[[2]])
#'
#' out[[1]]$r(100) |> as.numeric()
#' out[[2]]$r(100) |> as.numeric()
#'
#' @export
int_univariate_const_region = function(a, b, w, g)
{
	IntUnivariateConstRegion$new(a = a, b = b, w = w, g = g)
}

#' Integer Univariate Region with Constant Majorizer
#'
#' A region based on univariate intervals with a constant majorizer for the
#' weight function. This version is for integer supports.
#'
#' @field w Weight function for the target distribution.
#'
#' @export
IntUnivariateConstRegion = R6::R6Class(

classname = "IntUnivariateConstRegion",
inherit = UnivariateConstRegion,
portable = TRUE,
lock_class = FALSE,
private = list(
	a = NULL,
	b = NULL,
	g = NULL,
	log_w_max = NULL,
	log_w_min = NULL,
	log_prob = NULL
),
public = list(

w = NULL,

#' @param a Lower limit of interval.
#' @param b Upper limit of interval.
#' @param g An object created by \code{univariate_helper}.
#' @param w Weight function for the target distribution.
initialize = function(a, b, w, g)
{
	stopifnot(a <= b)
	stopifnot("univariate_helper" %in% class(g))

	private$a = a
	private$b = b
	private$g = g
	self$w = w

	private$log_w_max = self$optimize(maximize = TRUE)
	private$log_w_min = self$optimize(maximize = FALSE)

	# Compute g$p(b) - g$p(a) on the log scale
	if (a < b) {
		# Not a point mass
		private$log_prob = log_sub2_exp(g$p(b, log.p = TRUE), g$p(a, log.p = TRUE))
	} else {
		# Point mass
		private$log_prob = g$d(b, log = TRUE)
	}
},

#' @description
#' Generate a draw from \eqn{g_j} specific to this region.
#' @param n Number of draws to generate.
#' @return A list of draws, with one draw per list element.
r = function(n)
{
	# Compute g$q((pb - pa) * u + pa) on the log scale.
	# Some of the quantile functions seem to produce NaN if the probability
	# is is numerically larger than one: or zero on the log scale. We
	# replace these values of >= 0 with a negative of machine epsilon.
	u = runif(n)
	log_pa = private$g$p(private$a, log.p = TRUE)
	log_p = log_add2_exp(private$log_prob + log(u), rep(log_pa, n))
	log_p_adj = pmin(log_p, -.Machine$double.eps)
	x = private$g$q(log_p_adj, log.p = TRUE)
	as.list(x)
},

#' @description
#' Bifurcate this region into two regions. Use \code{x} as the bifurcation
#' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
#' @param x An optional bifurcation point.
bifurcate = function(x = NULL)
{
	a = private$a
	b = private$b

	if (is.null(x)) {
		if (is.infinite(a) && is.infinite(a) && a < 0 && b > 0) {
			# In this case, we have an interval (-Inf, Inf). Make a split at zero.
			x = 0
		} else if (is.infinite(a) && a < 0) {
			# Left endpoint is -Inf. Split based on right endpoint.
			x = b - abs(b) - 1
		} else if (is.infinite(b) && b > 0) {
			# Right endpoint is Inf. Split based on left endpoint.
			x = a + abs(a) + 1
		} else {
			# Both endpoints are finite. Take the midpoint.
			# For discrete intervals, only bifurcate at integers
			x = ceiling((a + b) / 2)
		}
	}

	s1 = IntUnivariateConstRegion$new(a = a, b = x, w = self$w, g = private$g)
	s2 = IntUnivariateConstRegion$new(a = x, b = b, w = self$w, g = private$g)
	list(s1, s2)
},

#' @description
#' Return a logical value indicating whether this region is bifurcatable.
is_bifurcatable = function()
{
	a = private$a
	b = private$b

	# A discrete interval is bifurcatable if there is at least one integer
	# between the limits
	if (is.infinite(a) || is.infinite(b)) {
		out = TRUE
	} else {
		out = ceiling((a + b) / 2) < b
	}

	return(out)
},

#' @description
#' Print a description of the region.
print = function()
{
	printf("Integer Univariate Const Region (%g, %g]\n", private$a, private$b)
},

#' @description
#' Maximize or minimize the function \eqn{w(x)} over this region. Optimization
#' is carried out with \code{optim} using arguments
#' \code{method = "Nelder-Mead"} and
#' \code{control = list(maxit = 100000, warn.1d.NelderMead = FALSE)}.
#' @param maximize logical; if \code{TRUE} do maximization. Otherwise do
#' minimization.
#' @param log logical; if \code{TRUE} return optimized value of \eqn{\log w(x)}.
#' Otherwise return optimized value of \eqn{w(x)}.
optimize = function(maximize = TRUE, log = TRUE)
{
	a = private$a
	b = private$b
	w = self$w

	method = "Nelder-Mead"
	control = list(maxit = 100000, warn.1d.NelderMead = FALSE)

	if (is.finite(a) && is.finite(b)) {
		# If both endpoints are finite, do brute force optimization. This may
		# not be desirable in some situations but should mostly work for
		# textbook problems ...
		x_seq = seq.int(ceiling(a + 0.0001), floor(b), by = 1)
		log_w_seq = w(x_seq, log = TRUE)
		if (maximize) {
			out = max(log_w_seq)
		} else {
			out = min(log_w_seq)
		}
	} else {
		# Transform to bounded interval, if necessary
		tx = function(x) {
			if (is.infinite(a) && is.infinite(b) && a < 0 && b > 0) {
				out = x
			} else if (is.infinite(a) && a < 0) {
				out = b*plogis(x)
			} else if (is.infinite(b) && b > 0) {
				out = exp(x) + a
			} else {
				out = (b-a) * plogis(x) + a
			}
			return(out)
		}

		f_opt = function(x) {
			w(tx(x), log = TRUE)
		}

		init = 0
		log_w_endpoints = c(w(a, log = TRUE), w(b, log = TRUE))
		log_w_endpoints = log_w_endpoints[!is.na(log_w_endpoints)]

		endpoint_pos_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints > 0)
		endpoint_neg_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints < 0)

		if (maximize && endpoint_pos_inf) {
			out = Inf
		} else if (maximize) {
			control$fnscale = -1
			opt_out = optim(init, f_opt, method = method, control = control)
			if (opt_out$convergence != 0) {
				warning("opt_out: convergence status was ", opt_out$convergence)
				browser()
			}
			out = max(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints)
		}

		if (!maximize && endpoint_neg_inf) {
			out = -Inf
		} else if (!maximize) {
			control$fnscale = 1
			opt_out = optim(init, f_opt, method = method, control = control)
			if (opt_out$convergence != 0) {
				warning("opt_out: convergence status was ", opt_out$convergence)
				browser()
			}
			out = min(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints)
		}
	}

	if (log) { return(out) } else { return(exp(out)) }
}

) # Close public
) # Close class
