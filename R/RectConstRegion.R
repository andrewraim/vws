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
	stopifnot(all(a <= b))
	stopifnot(all(unlist(Map(function(x) { "univariate_helper" %in% class(x) }, g))))
	stopifnot(k == length(a))
	stopifnot(k == length(b))

	private$a = a
	private$b = b
	private$g = g
	self$w = w

	private$log_w_max = self$optimize(maximize = TRUE)
	private$log_w_min = self$optimize(maximize = FALSE)

	# Compute Prob((a,b]) under dist g on the log scale
	lp = numeric(k)
	for (l in 1:k) {
		lpa = g[[l]]$p(a[l], log.p = TRUE)
		lpb = g[[l]]$p(b[l], log.p = TRUE)
		lp[l] = log_sub2_exp(lpb, lpa)
	}
	private$log_prob = sum(lp)
},

#' @description
#' Density function \eqn{g} for the base distribution.
#' @param x Density argument.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
d_base = function(x, log = FALSE)
{
	g = private$g
	k = length(g)
	stopifnot(length(x) == k)

	lp = numeric(k)
	for (l in 1:k) {
		lp[l] = g[[l]]$d(x[l], log = TRUE)
	}

	out = sum(lp)
	if (log) { return(out) } else { return(exp(out)) }
},

#' @description
#' Generate a draw from \eqn{g_j} specific to this region.
#' @param n Number of draws to generate.
#' @return A list of draws, with one draw per list element.
r = function(n)
{
	g = private$g
	a = private$a
	b = private$b
	k = length(g)
	x = matrix(NA, n, k)

	for (l in 1:k) {
		u = runif(n)
		lpa = g[[l]]$p(a[l], log.p = TRUE)
		lpb = g[[l]]$p(b[l], log.p = TRUE)
		log_prob_l = log_sub2_exp(lpb, lpa)
		log_p = log_add2_exp(log_prob_l + log(u), rep(lpa, n))
		x[,l] = g[[l]]$q(log_p, log.p = TRUE)
	}

	out = lapply(seq_len(n), function(i) x[i,])
	return(out)
},

# TBD: do we need to be able to handle a list of input values, or just one at a
# time? Make sure C++ and R interfaces are consistent!!

#' @description
#' Density of \eqn{g_j} specific to this region.
#' @param x Density argument.
d = function(x)
{
	g = private$g
	a = private$a
	b = private$b
	k = length(g)
	stopifnot(k == length(x))
	lp = rep(-Inf, k)

	for (l in 1:k) {
		if (a[l] < x[l] & x[l] <= b[l] & g[[l]]$s(x[l])) {
			lpa = g[[l]]$p(a[l], log.p = TRUE)
			lpb = g[[l]]$p(b[l], log.p = TRUE)
			lpx = g[[l]]$d(x[l], log = TRUE)
			lp[l] = lpx - log_sub2_exp(lpb, lpa)
		}
	}

	out = sum(lp)
	if (log) { return(out) } else { return(exp(out)) }
},

#' @description
#' Test if given \code{x} is in the support for the \eqn{g_j} specific to this
#' region.
#' @param x Density argument.
s = function(x)
{
	g = private$g
	k = length(g)
	a = private$a
	b = private$b
	stopifnot(k == length(x))
	pp = numeric(k)

	for (l in 1:k) {
		pp[l] = (a[l] < x[l] & x[l] <= b[l]) & g[[l]]$s(x[l])
	}

	return(prod(pp))
},

#' @description
#' Majorized weight function \eqn{\overline{w}_j} for this region.
#' @param x Argument to weight function.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
w_major = function(x, log = TRUE)
{
	if (!self$s(x)) {
		out = ifelse(log, -Inf, 0)
		return(out)
	}
	out = private$log_w_max
	if (log) { return(out) } else { return(exp(out)) }
},

midpoint = function()
{
	a = private$a
	b = private$b
	g = private$g
	k = length(g)
	x = numeric(k)

	for (l in 1:k) {
		if (is.infinite(a[l]) && is.infinite(a[l]) && a[l] < 0 && b[l] > 0) {
			# Here we have an interval (-Inf, Inf). Make a split at zero.
			x[l] = 0
		} else if (is.infinite(a[l]) && a[l] < 0) {
			# Left endpoint is -Inf. Split based on right endpoint.
			x[l] = b[l] - abs(b[l]) - 1
		} else if (is.infinite(b[l]) && b[l] > 0) {
			# Right endpoint is Inf. Split based on left endpoint.
			x[l] = a[l] + abs(a[l]) + 1
		} else {
			x[l] = (a[l] + b[l]) / 2
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
	g = private$g
	k = length(g)

	if (is.null(x)) {
		# cuts = self$midpoint()
		# pairs_left = list()
		# pairs_right = list()
		# reduction = rep(NA, k)

		# for (j in 1:k) {
		# 	a1 = a
		# 	a2 = a
		# 	b1 = b
		# 	b2 = b
		# 	a2[j] = cuts[j]
		# 	b1[j] = cuts[j]

		# 	s1 = RectConstRegion$new(a = a1, b = b1, w = self$w, g = private$g)
		# 	s2 = RectConstRegion$new(a = a2, b = b2, w = self$w, g = private$g)
		# 	pairs_left[[j]] = s1
		# 	pairs_right[[j]] = s2

		# 	lp0 = log_sub2_exp(self$xi_upper(), self$xi_lower())
		# 	lp1 = log_sub2_exp(s1$xi_upper(), s1$xi_lower())
		# 	lp2 = log_sub2_exp(s2$xi_upper(), s2$xi_lower())
		# 	reduction[j] = log_sub2_exp(lp0, log_add2_exp(lp1, lp2))
		# }

		# browser()
		# idx = r_categ(n = 1, p = reduction, log_p = TRUE)
		# out = list(pairs_left[[idx]], pairs_right[[idx]])
		# return(out)

		# Sample a cut orientation randomly
		cuts = self$midpoint()
		idx = sample(1:k, size = 1)
		a1 = a
		a2 = a
		b1 = b
		b2 = b
		a2[idx] = cuts[idx]
		b1[idx] = cuts[idx]
		s1 = RectConstRegion$new(a = a1, b = b1, w = self$w, g = g)
		s2 = RectConstRegion$new(a = a2, b = b2, w = self$w, g = g)
		out = list(s1, s2)
	} else {
		s1 = RectConstRegion$new(a = a, b = x, w = self$w, g = g)
		s2 = RectConstRegion$new(a = x, b = b, w = self$w, g = g)
		out = list(s1, s2)
	}

	return(out)
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
	out = private$log_w_max + private$log_prob
	ifelse(log, out, exp(out))
},

#' @description
#' The quantity \eqn{\underline{\xi}_j} for this region.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
xi_lower = function(log = TRUE)
{
	out = private$log_w_min + private$log_prob
	ifelse(log, out, exp(out))
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
	a = private$a
	b = private$b
	k = length(a)
	w = self$w

	method = "L-BFGS-B"
	f_opt = function(x) {
		w(x, log = TRUE)
	}

	control = list(maxit = 100000, trace = 0)
	control$fnscale = ifelse(maximize, -1, 1)

	init = self$midpoint()
	opt_out = optim(init, f_opt, method = method, lower = a, upper = b, control = control)
	if (opt_out$convergence != 0 && opt_out$convergence != 52) {
		warning("opt_out: convergence status was ", opt_out$convergence)
		browser()
	}
	# out = min(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints)
	out = opt_out$value


	if (FALSE) {
		method = "L-BFGS-B"
		# method = "Nelder-Mead"
		control = list(maxit = 100000, trace = 0)

		f_opt = function(x) {
			tx = inv_rect(x, a, b)
			w(tx, log = TRUE)
		}

		init = numeric(k)

		#init = self$midpoint()
		#log_w_endpoints = c(w(a, log = TRUE), w(b, log = TRUE))
		#log_w_endpoints = log_w_endpoints[!is.na(log_w_endpoints)]

		# endpoint_pos_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints > 0)
		# endpoint_neg_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints < 0)

		# if (maximize && endpoint_pos_inf) {
		#	out = Inf
		# } else if (maximize) {

		if (maximize) {
			control$fnscale = -1
			opt_out = optim(init, f_opt, method = method, control = control)
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
			opt_out = optim(init, f_opt, method = method, control = control)
			if (opt_out$convergence != 0) {
				warning("opt_out: convergence status was ", opt_out$convergence)
				browser()
			}
			# out = min(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints)
			out = opt_out$value
		}
	}

	if (log) { return(out) } else { return(exp(out)) }
}

) # Close public
) # Close class
