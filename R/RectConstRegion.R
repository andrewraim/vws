#' Univariate Region with Constant Majorizer
#'
#' An R6 class which represents a region based on univariate intervals with a
#' constant majorizer for the weight function. This version is for continuous
#' supports.
#'
#' @field w Weight function for the target distribution. Its expected interface
#' is \code{w(x, log = TRUE)} so that results are returned on the log-scale by
#' default.
#'
#' @examples
#' # Define base distribution and weight function
#' g = normal_univariate_helper(mean = 0, sd = 5)
#' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
#'
#' reg = RectConstRegion$new(-Inf, 10, w, g)
#' print(reg)
#'
#' out = reg$bifurcate(0)
#' print(out[[1]])
#' print(out[[2]])
#'
#' @export
RectConstRegion = R6::R6Class(

classname = "RectConstRegion",
portable = TRUE,
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

#' @param a Lower limits of intervals.
#' @param b Upper limits of intervals.
#' @param w Weight function for the target distribution.
#' @param g An list of objects created by \code{univariate_helper}.
initialize = function(a, b, w, g)
{
	k = length(g)
	stopifnot(length(a) == length(b))
	stopifnot(all(a <= b))
	stopifnot(all(unlist(Map(function(x) { "univariate_helper" %in% class(x) }, g))))
	stopifnot(k == length(a))

	private$a = a
	private$b = b
	private$g = g
	self$w = w

	private$log_w_max = self$optimize(maximize = TRUE)
	private$log_w_min = self$optimize(maximize = FALSE)

	# Compute g$p(b) - g$p(a) on the log scale
	lp = numeric(k)
	for (j in 1:k) {
		lp[j] = log_sub2_exp(g[[j]]$p(b[j], log.p = TRUE), g[[j]]$p(a[j], log.p = TRUE))
	}
	private$log_prob = sum(lp)
},

#' @description
#' Density function \eqn{g} for the base distribution.
#' @param x Density argument.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
d_base = function(x, log = FALSE)
{
	k = length(private$g)
	stopifnot(length(x) == k)

	out = numeric(k)
	for (j in 1:k) {
		out[j] = private$g[[j]]$d(x[j], log = log)
	}
},

#' @description
#' Generate a draw from \eqn{g_j} specific to this region.
#' @param n Number of draws to generate.
#' @return A list of draws, with one draw per list element.
r = function(n)
{
	k = length(private$g)
	x = matrix(NA, n, k)

	for (j in 1:k) {
		u = runif(n)
		log_pa = private$g[[j]]$p(private$a[j], log.p = TRUE)
		log_p = log_add2_exp(private$log_prob + log(u), rep(log_pa, n))
		x[,j] = private$g[[j]]$q(log_p, log.p = TRUE)
	}

	out = lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
},

# TBD: do we need to be able to handle a list of input values, or just one at a
# time? Make sure C++ and R interfaces are consistent!!

#' @description
#' Density of \eqn{g_j} specific to this region.
#' @param x Density argument.
d = function(x)
{
	k = length(private$g)
	n = length(x)
	xx = matrix(unlist(x), n, k, byrow = TRUE)
	lp = matrix(-Inf, n, k)

	for (j in 1:k) {
		idx = which(a[j] < xx[,j] & xx[,j] <= b[j])
		lp[,j] = private$g[[j]]$d(x[idx], log = TRUE) -
			log_sub2_exp(private$g[[j]]$p(private$b[j], log.p = TRUE),
				private$g[[j]]$p(private$a[j], log.p = TRUE))
	}

	out = apply(lp, 1, log_sum_exp)
	if (log) { return(out) } else { return(exp(out)) }
},

#' @description
#' Test if given \code{x} is in the support for the \eqn{g_j} specific to this
#' region.
#' @param x Density argument.
s = function(x)
{
	k = length(private$g)
	n = length(x)
	xx = matrix(unlist(x), n, k, byrow = TRUE)
	pp = matrix(-Inf, n, k)

	for (j in 1:k) {
		pp[,j] = private$a[j] < xx[,j] & xx[,j] <= private$b[j] & private$g[[j]]$s(xx[,j])
	}

	out = apply(pp, 1, prod)
	if (log) { return(out) } else { return(exp(out)) }
},

#' @description
#' Majorized weight function \eqn{\overline{w}_j} for this region.
#' @param x Argument to weight function.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
w_major = function(x, log = TRUE)
{
	if (!private$s(x)) {
		out = ifelse(log, -Inf, 0)
		return(out)
	}
	out = private$log_w_max
	if (log) { return(out) } else { return(exp(out)) }
},

midpoint = function()
{
	k = length(private$g)
	a = private$a
	b = private$b

	k = length(private$g)
	x = numeric(k)

	for (j in 1:k) {
		if (is.infinite(a[j]) && is.infinite(a[j]) && a[j] < 0 && b[j] > 0) {
			# Here we have an interval (-Inf, Inf). Make a split at zero.
			x[j] = 0
		} else if (is.infinite(a[j]) && a[j] < 0) {
			# Left endpoint is -Inf. Split based on right endpoint.
			x[j] = b[j] - abs(b[j]) - 1
		} else if (is.infinite(b[j]) && b[j] > 0) {
			# Right endpoint is Inf. Split based on left endpoint.
			x[j] = a[j] + abs(a[j]) + 1
		} else {
			x[j] = (a[j] + b[j]) / 2
		}
	}

	return(x)
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
		x = self$midpoint()
	}

	s1 = RectConstRegion$new(a = a, b = x, w = self$w, g = private$g)
	s2 = RectConstRegion$new(a = x, b = b, w = self$w, g = private$g)
	list(s1, s2)
},

#' @description
#' Return a logical value indicating whether this region is bifurcatable.
is_bifurcatable = function()
{
	return(TRUE)
},

#' @description
#' The quantity \eqn{\overline{\xi}_j} for this region.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
xi_upper = function(log = TRUE)
{
	log_xi_upper = private$log_w_max + private$log_prob
	ifelse(log, log_xi_upper, exp(log_xi_upper))
},

#' @description
#' The quantity \eqn{\underline{\xi}_j} for this region.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
xi_lower = function(log = TRUE)
{
	log_xi_lower = private$log_w_min + private$log_prob
	ifelse(log, log_xi_lower, exp(log_xi_lower))
},

#' @description
#' The upper limits \eqn{\alpha_j} for this region.
upper = function()
{
	private$b
},

#' @description
#' The lower limits \eqn{\alpha_{j-1}} for this region.
lower = function()
{
	private$a
},

#' @description
#' A string that describes the region.
description = function()
{
	comps = sprintf("(%g, %g]", private$a, private$b)
	paste(comps, collapse = " x ")
},

#' @description
#' Print a description of the region.
print = function()
{
	comps = sprintf("(%g, %g]", private$a, private$b)
	printf("Rect Const Region %s\n", paste(comps, collapse = " x "))
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
	k = length(private$g)
	a = private$a
	b = private$b
	w = self$w

	method = "L-BFGS-B"
	control = list(maxit = 100000, trace = 6)

	f_opt = function(x) { w(x, log = TRUE) }

	init = self$midpoint()
	#log_w_endpoints = c(w(a, log = TRUE), w(b, log = TRUE))
	#log_w_endpoints = log_w_endpoints[!is.na(log_w_endpoints)]

	# endpoint_pos_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints > 0)
	# endpoint_neg_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints < 0)

	# if (maximize && endpoint_pos_inf) {
	#	out = Inf
	# } else if (maximize) {

	if (maximize) {
		control$fnscale = -1
		opt_out = optim(init, f_opt, method = method, lower = a, upper = b, control = control)
		if (opt_out$convergence != 0) {
			warning("opt_out: convergence status was ", opt_out$convergence)
			browser()
		}
		# out = max(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints)
		out = opt_out$value
	}

	# if (!maximize && endpoint_neg_inf) {
	#	out = -Inf
	# } else if (!maximize) {
	if (!maximize) {
		control$fnscale = 1
		opt_out = optim(init, f_opt, method = method, lower = a, upper = b, control = control)
		if (opt_out$convergence != 0) {
			warning("opt_out: convergence status was ", opt_out$convergence)
			browser()
		}
		# out = min(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints)
		out = opt_out$value
	}

	if (log) { return(out) } else { return(exp(out)) }
}

) # Close public
) # Close class
