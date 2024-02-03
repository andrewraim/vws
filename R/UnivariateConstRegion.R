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

UnivariateConstRegion$methods(d_base = function(x, log = FALSE)
{
	g$d(x, log = log)
})

UnivariateConstRegion$methods(w = function(x, log = TRUE)
{
	w(x, log = log)
})

UnivariateConstRegion$methods(r = function(n)
{
	# Generate a draw from $g_j$; i.e., the density $g$ truncated to this region.
	# Compute g$q((pb - pa) * u + pa) on the log scale
	u = runif(n)
	log_pa = g$p(a, log.p = TRUE)
	log_p = log_add2_exp(log_prob + log(u), rep(log_pa,n))
	x = g$q(log_p, log.p = TRUE)
	return(as.list(x))
})

UnivariateConstRegion$methods(d = function(x)
{
	n = length(x)
	out = rep(-Inf, n)
	idx = which(a < x & x <= b)
	out[idx] = d(x[idx], log = TRUE) -
		log_sub2_exp(p(b, log.p = TRUE), p(a, log.p = TRUE))

	if (log) { return(out) } else { return(exp(out)) }
})

UnivariateConstRegion$methods(in_support = function(x)
{
	a < x & x <= b & g$in_support(x)
})

UnivariateConstRegion$methods(w_major = function(x, log = TRUE)
{
	if (!g$in_support(x)) {
		out = ifelse(log, -Inf, 0)
		return(out)
	}
	if (log) { return(log_w_max) } else { return(exp(log_w_max)) }
})

UnivariateConstRegion$methods(bifurcate = function(x = NULL)
{
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

	s1 = univariate_const_region(a, x, w, g)
	s2 = univariate_const_region(x, b, w, g)
	list(s1, s2)
})

UnivariateConstRegion$methods(is_bifurcatable = function(x)
{
	return(TRUE)
})

UnivariateConstRegion$methods(xi_upper = function(log = TRUE)
{
	log_xi_upper = log_w_max + log_prob
	ifelse(log, log_xi_upper, exp(log_xi_upper))
})

UnivariateConstRegion$methods(xi_lower = function(log = TRUE)
{
	log_xi_lower = log_w_min + log_prob
	ifelse(log, log_xi_lower, exp(log_xi_lower))
})

UnivariateConstRegion$methods(description = function()
{
	sprintf("(%g, %g]", a, b)
})

UnivariateConstRegion$methods(show = function()
{
	printf("Univariate Const Region (%g, %g]\n", a, b)
})
