#' Univariate Region with Constant Majorizer
#'
#' A region based on univariate intervals with a constant majorizer for the
#' weight function.
#'
#' @param a Lower knot of region.
#' @param b Upper knot of region.
#' @param w Weight function.
#' @param g Object that encapsulates base distribution.
#'
#' @examples
#' # Define base distribution and weight function
#' g = normal_univariate_helper(mean = 0, sd = 5)
#' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
#'
#' reg = univariate_const_region(-Inf, 10, w, g)
#' print(reg)
#'
#' out = reg$bifurcate(0)
#' print(out[[1]])
#' print(out[[2]])
#'
#' @name UnivariateConstRegion
#' @export
univariate_const_region = function(a, b, w, g)
{
	stopifnot(a <= b)
	stopifnot("univariate_helper" %in% class(g))

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

	control1 = list(fnscale = -1, maxit = 100000, warn.1d.NelderMead = FALSE)
	opt1_out = optim(init, f_opt, method = "Nelder-Mead", control = control1)
	if (opt1_out$convergence != 0) {
		warning("opt1_out: convergence status was ", opt1_out$convergence)
	}
	log_w_max = max(f_opt(floor(opt1_out$par)), f_opt(ceiling(opt1_out$par)), log_w_endpoints)

	control2 = list(fnscale = 1, maxit = 100000, warn.1d.NelderMead = FALSE)
	opt2_out = optim(init, f_opt, method = "Nelder-Mead", control = control2)
	if (opt2_out$convergence != 0) {
		warning("opt2_out: convergence status was ", opt2_out$convergence)
	}
	log_w_min = min(f_opt(floor(opt2_out$par)), f_opt(ceiling(opt2_out$par)), log_w_endpoints)

	# Compute g$p(b) - g$p(a) on the log scale
	if (a < b) {
		# Not a point mass
		log_prob = log_sub2_exp(g$p(b, log.p = TRUE), g$p(a, log.p = TRUE))
	} else {
		# Point mass
		log_prob = g$d(b, log = TRUE)
	}

	UnivariateConstRegion$new(a = a, b = b, w = w, g = g, log_w_max = log_w_max,
		log_w_min = log_w_min, log_prob = log_prob)
}

UnivariateConstRegion = R6::R6Class(
	classname = "UnivariateConstRegion",
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
		w = NULL
	)
)

UnivariateConstRegion$set("public", "initialize", function(a, b, w, g,
	log_w_max, log_w_min, log_prob)
{
	private$a = a
	private$b = b
	private$g = g
	private$log_w_max = log_w_max
	private$log_w_min = log_w_min
	private$log_prob = log_prob
	self$w = w
})

UnivariateConstRegion$set("public", "d_base", function(x, log = FALSE)
{
	private$g$d(x, log = log)
})

UnivariateConstRegion$set("public", "r", function(n)
{
	# Generate a draw from $g_j$; i.e., the density $g$ truncated to this region.
	# Compute g$q((pb - pa) * u + pa) on the log scale
	u = runif(n)
	log_pa = private$g$p(private$a, log.p = TRUE)
	log_p = log_add2_exp(private$log_prob + log(u), rep(log_pa, n))
	x = private$g$q(log_p, log.p = TRUE)
	as.list(x)
})

UnivariateConstRegion$set("public", "d", function(x)
{
	n = length(x)
	out = rep(-Inf, n)
	idx = which(a < x & x <= b)
	out[idx] = private$g$d(x[idx], log = TRUE) -
		log_sub2_exp(private$g$p(private$b, log.p = TRUE), private$g$p(private$a, log.p = TRUE))

	if (log) { return(out) } else { return(exp(out)) }
})

UnivariateConstRegion$set("public", "s", function(x)
{
	private$a < x & x <= private$b & private$g$s(x)
})

UnivariateConstRegion$set("public", "w_major", function(x, log = TRUE)
{
	if (!private$g$s(x)) {
		out = ifelse(log, -Inf, 0)
		return(out)
	}
	out = private$log_w_max
	if (log) { return(out) } else { return(exp(out)) }
})

UnivariateConstRegion$set("public", "bifurcate", function(x = NULL)
{
	a = private$a
	b = private$b

	if (is.null(x)) {
		if (is.infinite(a) && is.infinite(b) && a < 0 && b > 0) {
			# In this case, we have an interval (-Inf, Inf). Make a split at zero.
			x = 0
		} else if (is.infinite(a) && a < 0) {
			# Left endpoint is -Inf. Split based on right endpoint.
			x = b - abs(b) - 1
		} else if (is.infinite(b) && b > 0) {
			# Right endpoint is Inf. Split based on left endpoint.
			x = a + abs(a) + 1
		} else {
			x = (a + b) / 2
		}
	}

	s1 = univariate_const_region(a, x, self$w, private$g)
	s2 = univariate_const_region(x, b, self$w, private$g)
	list(s1, s2)
})

UnivariateConstRegion$set("public", "is_bifurcatable", function(x)
{
	return(TRUE)
})

UnivariateConstRegion$set("public", "xi_upper", function(log = TRUE)
{
	log_xi_upper = private$log_w_max + private$log_prob
	ifelse(log, log_xi_upper, exp(log_xi_upper))
})

UnivariateConstRegion$set("public", "xi_lower", function(log = TRUE)
{
	log_xi_lower = private$log_w_min + private$log_prob
	ifelse(log, log_xi_lower, exp(log_xi_lower))
})

UnivariateConstRegion$set("public", "description", function()
{
	sprintf("(%g, %g]", private$a, private$b)
})

UnivariateConstRegion$set("public", "print", function()
{
	printf("Univariate Const Region (%g, %g]\n", private$a, private$b)
})

UnivariateConstRegion$lock()
