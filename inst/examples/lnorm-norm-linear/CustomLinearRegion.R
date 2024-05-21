# Define Class

CustomLinearRegion = R6::R6Class(
	classname = "CustomLinearRegion",
	portable = TRUE,
	lock_class = FALSE,
	private = list(
		a = NULL,
		b = NULL,
		mu = NULL,
		sigma2 = NULL,
		z = NULL,
		lambda2 = NULL,
		beta0_min = NULL,
		beta1_min = NULL,
		beta0_max = NULL,
		beta1_max = NULL
	)
)

## Initializer
CustomLinearRegion$set("public", "initialize", function(a, b, mu, sigma2, z, lambda2)
{
	stopifnot(a <= b)

	private$a = a
	private$b = b
	private$mu = mu
	private$sigma2 = sigma2
	private$z = z
	private$lambda2 = lambda2

	w = self$w

	# First derivative of log w(x)
	d_log_w = function(x) { -1/x * (1 + (log(x) - mu) / sigma2) }

	# MGF the truncated and reweighted g
	mgf = function(s, log = FALSE, tol = 1e-6) {
		# If we are truncating way into the upper tail of a distribution, working
		# with the complement of the CDF helps to retain precision. Otherwise,
		# work with the CDF function.

		lp_num_a = pnorm(a, mean = z + s*lambda2, sd = sqrt(lambda2), log.p = TRUE)
		lp_num_b = pnorm(b, mean = z + s*lambda2, sd = sqrt(lambda2), log.p = TRUE)
		clp_num_a = pnorm(a, mean = z + s*lambda2, sd = sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		clp_num_b = pnorm(b, mean = z + s*lambda2, sd = sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		lp_num = ifelse(lp_num_a > log1p(-tol),
			log_sub2_exp(clp_num_a, clp_num_b),
			log_sub2_exp(lp_num_b, lp_num_a)
		)

		lp_den_a = pnorm(a, mean = z, sd = sqrt(lambda2), log.p = TRUE)
		lp_den_b = pnorm(b, mean = z, sd = sqrt(lambda2), log.p = TRUE)
		clp_den_a = pnorm(a, mean = z, sd = sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		clp_den_b = pnorm(b, mean = z, sd = sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		lp_den = ifelse(lp_den_a > log1p(-tol),
			log_sub2_exp(clp_den_a, clp_den_b),
			log_sub2_exp(lp_den_b, lp_den_a)
		)

		out = lp_num - lp_den + s*z + s^2 * lambda2 / 2
		if (log) { return(out) } else { return(exp(out)) }
	}

	obj_line = function(x) {
		gr = d_log_w(x)
		w(x, log = TRUE) - x * gr + mgf(gr, log = TRUE)
	}

	l_concave = log(a) < mu - sigma2 + 1
	r_convex = log(b) > mu - sigma2 + 1
	if (l_concave && r_convex) {
		msg = sprintf("%s y = %g %s\n",
			"Partition your region so that",
			exp(mu - sigma2 + 1),
			"is not in the interior")
		stop(msg)
	}

	is_concave = l_concave
	is_convex = r_convex

	optim_out = optimize(f = obj_line, interval = c(a, b), maximum = FALSE)
	c_star = optim_out$minimum

	if (is_concave) {
		# log w(x) is concave

		# For the majorizer
		beta0_max = w(c_star) - c_star*d_log_w(c_star)
		beta1_max = d_log_w(c_star)

		# For the minorizer
		c = c(w(a), w(b))
		if (all(is.finite(c))) {
			A = matrix(c(1,1,a,b), 2, 2)
			x = solve(A, c)
			beta0_min = x[1]
			beta1_min = x[2]
		} else {
			# Special handling when log w(x) is not finite at one of the endpoints.
			beta0_min = -Inf
			beta1_min = 0
		}
	} else if (is_convex) {
		# log w(x) is convex

		# For the majorizer
		beta0_min = w(c_star) - c_star*d_log_w(c_star)
		beta1_min = d_log_w(c_star)

		# For the majorizer
		c = c(w(a), w(b))
		if (all(is.finite(c))) {
			A = matrix(c(1,1,a,b), 2, 2)
			x = solve(A, c)
			beta0_max = x[1]
			beta1_max = x[2]
		} else {
			# Special handling when log w(x) is not finite at one of the endpoints.
			beta0_max = max(c)
			beta1_max = 0
		}
	} else {
		# log w(x) is constant
		beta0_min = 0
		beta1_min = 0
		beta0_max = 0
		beta1_max = 0
	}

	if (FALSE) {
		# Plot on log-scale
		curve(w(x, log = TRUE), xlim = c(a,b))
		abline(coef = c(beta0_max, beta1_max), lty = 2, col = "blue")
		abline(coef = c(beta0_min, beta1_min), lty = 2, col = "red")

		# Plot on original scale
		w_maj = function(x) { exp(beta0_max + beta1_max * x) }
		w_min = function(x) { exp(beta0_min + beta1_min * x) }
		curve(w(x, log = FALSE), xlim = c(a, 2*a))
		curve(w_maj, add = TRUE, lty = 2)
		curve(w_min, add = TRUE, lty = 2)
	}

	private$beta0_max = beta0_max
	private$beta1_max = beta1_max
	private$beta0_min = beta0_min
	private$beta1_min = beta1_min

	if (is.na(self$xi_upper())) browser()
	if (is.na(self$xi_lower())) browser()
	if (self$xi_upper() < self$xi_lower()) browser()
})

## Base Density
## Define the base distribution density as a method.
CustomLinearRegion$set("public", "d_base", function(x, log = FALSE)
{
	dnorm(unlist(x), mean = private$z, sd = sqrt(private$lambda2), log = log)
})

## Weight Function
## Define the weight function as a method.
CustomLinearRegion$set("public", "w", function(x, log = TRUE)
{
	mu = private$mu
	sigma2 = private$sigma2
	out = -log(x) - (log(x) - mu)^2 / (2*sigma2) + log(x > 0)
	out[x == 0] = -Inf
	if (log) { return(out) } else { return(exp(out)) }
})

## Generate from Region
## Generate random variables from the truncated & reweighted base distribution for
## this region.
CustomLinearRegion$set("public", "r", function(n)
{
	mean = private$z + private$beta1_max * private$lambda2
	sd = sqrt(private$lambda2)
	r_truncated(n, lo = private$a, hi = private$b, pf = pnorm, qf = qnorm,
		mean = mean, sd = sd) |> as.list()
})

## Density for Region
## Density for the truncated & reweighted base distribution on this region.
CustomLinearRegion$set("public", "d", function(x, log = FALSE)
{
	mean = private$z + private$beta1_max * private$lambda2
	sd = sqrt(private$lambda2)
	d_truncated(unlist(x), lo = private$a, hi = private$b, df = dnorm,
		pf = pnorm, mean = mean, sd = sd, log = log)
})

CustomLinearRegion$set("public", "s", function(x)
{
	private$a < x & x <= private$b
})

CustomLinearRegion$set("public", "w_major", function(x, log = TRUE)
{
	out = rep(-Inf, length(x))
	idx = which(self$s(x))
	out[idx] = private$beta0_max + private$beta1_max*x[idx]
	if (log) { return(out) } else { return(exp(out)) }
})

CustomLinearRegion$set("public", "is_bifurcatable", function()
{
	return(TRUE)
})

CustomLinearRegion$set("public", "midpoint", function()
{
	a = private$a
	b = private$b

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

	return(x)
})

CustomLinearRegion$set("public", "bifurcate", function(x = NULL, ...)
{
	a = private$a
	b = private$b
	mu = private$mu
	sigma2 = private$sigma2
	z = private$z
	lambda2 = private$lambda2

	if (is.null(x)) {
		x = self$midpoint()
	}

	s1 = CustomLinearRegion$new(a, x, mu, sigma2, z, lambda2)
	s2 = CustomLinearRegion$new(x, b, mu, sigma2, z, lambda2)
	list(s1, s2)
})

## Upper and Lower Integrals
CustomLinearRegion$set("public", "xi_upper", function(log = TRUE, tol = 1e-6)
{
	a = private$a
	b = private$b
	mu = private$mu
	sigma2 = private$sigma2
	z = private$z
	lambda2 = private$lambda2
	beta0_max = private$beta0_max
	beta1_max = private$beta1_max

	lp_a = pnorm(a, z + beta1_max*lambda2, sqrt(lambda2), log.p = TRUE)
	if (lp_a > log1p(-tol)) {
		# If we are truncating way into the upper tail of a distribution, working
		# with the complement of the CDF helps to retain precision.
		clp_a = pnorm(a, z + beta1_max*lambda2, sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		clp_b = pnorm(b, z + beta1_max*lambda2, sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		out = beta0_max + beta1_max * z + beta1_max^2 * lambda2 / 2 +
			log_sub2_exp(clp_a, clp_b)
	} else {
		# Otherwise, work with the CDF function.
		lp_b = pnorm(b, z + beta1_max*lambda2, sqrt(lambda2), log.p = TRUE)
		out = beta0_max + beta1_max * z + beta1_max^2 * lambda2 / 2 +
			log_sub2_exp(lp_b, lp_a)
	}

	if (log) { return(out) } else { return(exp(out)) }
})

CustomLinearRegion$set("public", "xi_lower", function(log = TRUE, tol = 1e-6)
{
	a = private$a
	b = private$b
	mu = private$mu
	sigma2 = private$sigma2
	z = private$z
	lambda2 = private$lambda2
	beta0_min = private$beta0_min
	beta1_min = private$beta1_min

	lp_a = pnorm(a, z + beta1_min*lambda2, sqrt(lambda2), log.p = TRUE)
	if (lp_a > log1p(-tol)) {
		# If we are truncating way into the upper tail of a distribution, working
		# with the complement of the CDF helps to retain precision.
		clp_a = pnorm(a, z + beta1_min*lambda2, sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		clp_b = pnorm(b, z + beta1_min*lambda2, sqrt(lambda2), log.p = TRUE, lower.tail = FALSE)
		out = beta0_min + beta1_min * z + beta1_min^2 * lambda2 / 2 +
			log_sub2_exp(clp_a, clp_b)
	} else {
		# Otherwise, work with the CDF function.
		lp_b = pnorm(b, z + beta1_min*lambda2, sqrt(lambda2), log.p = TRUE)
		out = beta0_min + beta1_min * z + beta1_min^2 * lambda2 / 2 +
			log_sub2_exp(lp_b, lp_a)
	}

	# If we hopelessly run out of precision and the result is larger than
	# xi_upper, set xi_lower to xi_upper.
	out = min(out, self$xi_upper(log = TRUE))

	if (log) { return(out) } else { return(exp(out)) }
})

#' @description
#' The upper limit \eqn{\alpha_j} for this region.
CustomLinearRegion$set("public", "upper", function()
{
	private$b
})

#' @description
#' The upper limit \eqn{\alpha_j} for this region.
CustomLinearRegion$set("public", "lower", function()
{
	private$a
})

## Print Method
CustomLinearRegion$set("public", "description", function()
{
	sprintf("(%g, %g]", private$a, private$b)
})

CustomLinearRegion$set("public", "print", function()
{
	printf("Custom Linear Region (%g, %g]\n", private$a, private$b)
})
