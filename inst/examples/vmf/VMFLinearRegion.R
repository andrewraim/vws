# Define Class

VMFLinearRegion = R6::R6Class(
	classname = "VMFLinearRegion",
	portable = TRUE,
	lock_class = FALSE,
	private = list(
		a = NULL,
		b = NULL,
		dim = NULL,
		kappa = NULL,
		beta0_min = NULL,
		beta1_min = NULL,
		beta0_max = NULL,
		beta1_max = NULL
	)
)

## Initializer
VMFLinearRegion$set("public", "initialize", function(a, b, kappa, d)
{
	stopifnot(a <= b)

	private$a = a
	private$b = b
	private$dim = d
	private$kappa = kappa

	w = self$w

	# First derivative of log-weight function
	d_log_w = function(x) {
		-(d-3) * x / (1 - x^2)
	}

	# MGF the truncated and reweighted g
	mgf = function(s, log = FALSE)  {
		kappa_s = s + kappa
		if (kappa_s < 0) {
			out = log(kappa) - log(-kappa_s) +
				log_sub2_exp(kappa_s * a, kappa_s * b) -
				log_sub2_exp(kappa * b, kappa * a)
		} else {
			out = log(kappa) - log(kappa_s) +
				log_sub2_exp(kappa_s * b, kappa_s * a) -
				log_sub2_exp(kappa * b, kappa * a)
		}
		if (log) { return(out) } else { return(exp(out)) }
	}

	obj_line = function(x) {
		w(x, log = TRUE) - x * d_log_w(x) + mgf(x, log = TRUE)
	}

	# Determine if log(w(x)) is convex, concave, or a constant. Majorization
	# and minorization choices will depend on this.
	if (d > 3) {
		# log(w(x)) is concave. Any tangent line will lie above the curve.

		# For the majorizer
		optim_out = optimize(f = obj_line, interval = c(a, b), maximum = FALSE)
		c_star = optim_out$minimum
		beta0_max = w(c_star) - c_star*d_log_w(c_star)
		beta1_max = d_log_w(c_star)

		# For the minorizer
		A = matrix(c(1,1,a,b), 2, 2)
		c = c(w(a), w(b))
		x = solve(A, c)
		beta0_min = x[1]
		beta1_min = x[2]
	} else if (d < 3) {
		# log(w(x)) is convex. Make a line that passes through a and b.

		# For the majorizer
		A = matrix(c(1,1,a,b), 2, 2)
		c = c(w(a), w(b))
		x = solve(A, c)
		beta0_max = x[1]
		beta1_max = x[2]

		# For the minorizer
		optim_out = optimize(f = obj_line, interval = c(a, b), maximum = TRUE)
		c_star = optim_out$maximum
		beta0_min = w(c_star) - c_star*d_log_w(c_star)
		beta1_min = d_log_w(c_star)
	} else {
		# log(w(x)) is constant with value zero
		beta0_max = 0
		beta1_max = 0
		beta0_min = 0
		beta1_min = 0
	}

	private$beta0_max = beta0_max
	private$beta1_max = beta1_max
	private$beta0_min = beta0_min
	private$beta1_min = beta1_min
})

## Base Density
## Define the base distribution density as a method.
VMFLinearRegion$set("public", "d_base", function(x, log = FALSE)
{
	out = private$kappa*x + log(private$kappa) -
		vws::log_sub2_exp(private$kappa,-private$kappa) + log(-1 < x & x < 1)
	if (log) { return(out) } else { return(exp(out)) }
})

## Weight Function
## Define the weight function as a method.
VMFLinearRegion$set("public", "w", function(x, log = TRUE)
{
	out = (private$dim - 3) / 2 * log(1 - x^2)
	if (log) { return(out) } else { return(exp(out)) }
})

## Generate from Region
## Generate random variables from the truncated & reweighted base distribution for
## this region.
VMFLinearRegion$set("public", "r", function(n)
{
	a = private$a
	b = private$b
	kappa = private$kappa
	beta1_max = private$beta1_max

	u = runif(n)
	kappa_j = kappa + beta1_max
	x = vws::q_truncated(u, lo = a, hi = b, pf = p_base, qf = q_base,
		kappa = kappa_j)
	as.list(x)
})

## Density for Region
## Density for the truncated & reweighted base distribution on this region.
VMFLinearRegion$set("public", "d", function(x, log = FALSE)
{
	a = private$a
	b = private$b
	kappa = private$kappa
	beta1_max = private$beta1_max

	kappa_j = beta1_max + kappa
	if (kappa_j > 0) {
		out = log(kappa_j) + kappa_j*x -
			vws:::log_sub2_exp(kappa_j*b, kappa_j*a) +
			log(a < x & x <= b)
	} else {
		out = log(-kappa_j) + kappa_j*x -
			vws:::log_sub2_exp(kappa_j*a, kappa_j*b) +
			log(a < x & x <= b)
	}
	if (log) { return(out) } else { return(exp(out)) }
})

## Support Indicator
## Define the base distribution density as a method.

VMFLinearRegion$set("public", "s", function(x)
{
	private$a < x & x <= private$b
})

## Majorized Weight Function
VMFLinearRegion$set("public", "w_major", function(x, log = TRUE)
{
	out = rep(-Inf, length(x))
	idx = which(self$s(x))
	out[idx] = private$beta0_max + private$beta1_max*x[idx]
	if (log) { return(out) } else { return(exp(out)) }
})

## Bifurcation
VMFLinearRegion$set("public", "is_bifurcatable", function()
{
	return(TRUE)
})

VMFLinearRegion$set("public", "bifurcate", function(x = NULL, ...)
{
	a = private$a
	b = private$b

	if (is.null(x)) {
		if (is.infinite(a) && is.infinite(a) && a < 0 && b > 0) {
			# Here we have an interval (-Inf, Inf). Make a split at zero.
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

	s1 = VMFLinearRegion$new(a, x, private$kappa, private$dim)
	s2 = VMFLinearRegion$new(x, b, private$kappa, private$dim)
	list(s1, s2)
})

## Upper and Lower Integrals
VMFLinearRegion$set("public", "xi_upper", function(log = TRUE)
{
	a = private$a
	b = private$b
	kappa = private$kappa
	beta0_max = private$beta0_max
	beta1_max = private$beta1_max

	kappa_j = beta1_max + kappa
	if (kappa_j > 0) {
		out = log(kappa) + beta0_max - log(kappa_j) -
			vws::log_sub2_exp(kappa,-kappa) +
			vws::log_sub2_exp(kappa_j*b, kappa_j*a)
	} else {
		out = log(kappa) + beta0_max - log(-kappa_j) -
			vws::log_sub2_exp(kappa,-kappa) +
			vws::log_sub2_exp(kappa_j*a, kappa_j*b)
	}
	if (log) { return(out) } else { return(exp(out)) }
})

VMFLinearRegion$set("public", "xi_lower", function(log = TRUE)
{
	a = private$a
	b = private$b
	dim = private$dim
	kappa = private$kappa
	beta0_max = private$beta0_max
	beta1_max = private$beta1_max

	lb = p_target(b, kappa, dim, log.p = TRUE)
	la = p_target(a, kappa, dim, log.p = TRUE)
	lnc = nc_target(kappa, dim, log = TRUE)
	out = vws::log_sub2_exp(lb, la) + lnc

	log_xi_upper = self$xi_upper(log = TRUE)

	if (log_xi_upper < out) {
		# This condition can happen numerically. If it occurs, just take lower
		# to be equal to upper.
		out = log_xi_upper
	}
	if (log) { return(out) } else { return(exp(out)) }
})

## Print Method
VMFLinearRegion$set("public", "description", function()
{
	sprintf("(%g, %g]", private$a, private$b)
})

VMFLinearRegion$set("public", "print", function()
{
	printf("VMF Linear Region (%g, %g]\n", private$a, private$b)
})
