Rcpp::sourceCpp("draw-ln-norm.cpp")
Rcpp::sourceCpp("ln-norm-proposal.cpp")

#' r_mvnorm_prec
#'
#' Generate one multivariate normal variate using the specified mean vector and
#' precision matrix.
#'
#' @param mu Vector which represents the mean.
#' @param Omega Matrix which represents the precision.
#'
#' @returns A vector
#' @export
r_mvnorm_prec = function(mu, Omega)
{
	stopifnot(length(mu) == nrow(Omega))
	stopifnot(nrow(Omega) == ncol(Omega))

    k = length(mu)
    z = rnorm(k)
    A = chol(Omega)
    out = solve(A, z) + mu

	return(out)
}

#' printf
#'
#' Print formatted output to the console.
#'
#' @param fmt A format string as in [base::sprintf].
#' @param ... Values to be passed to fmt as in [base::sprintf].
#'
#' @returns None
#'
#' @export
printf = function (fmt, ...)
{
    cat(sprintf(fmt, ...))
}

#' control_inner
#'
#' List of arguments to VWS sampler within Gibbs.
#'
#' @param tol_suff The tolerance \eqn{\epsilon_1}.
#' @param tol_merge The tolerance \eqn{\epsilon_2}.
#' @param max_rejects Maximum number of rejects allowed before halting.
#' @param tune Number of Gibbs sampling steps where self-tuning will be used.
#' Only has an effect when `method = "vws-tune"`.
#' @param N Maximum number of regions allowed in a proposal. Only has an effect
#' when `method = "vws-basic"`.
#' @param method The method to use: either `"vws-tune"` or `"vws-basic"`.
#'
#' @returns
#' A list with the data.
#'
#' @export
control_inner = function(tol_suff = 1e-2, tol_merge = 1e-4,
	max_rejects = 1e6, tune = 1e6, N = 50, method = c("vws-tune", "vws-basic"))
{
	ret = list(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = max_rejects, tune = tune, N = N,
		method = match.arg(method))
	class(ret) = "control_inner"
	return(ret)
}

#' control_gibbs
#'
#' List of arguments to Gibbs sampler
#'
#' @param R Length of MCMC chain.
#' @param burn Number of draws to burn.
#' @param thin The period for saving samples.
#' @param report The period for printing progress to the console.
#' @param save_latent An integer vector whose elements are a subset of `1:n`.
#' Draws of \eqn{y} will be saved for these indices.
#' @param inner An object from [control_inner].
#'
#' @returns
#' A list with the data.
#'
#' @export
control_gibbs = function(R = 1000, burn = 0, thin = 1, report = R+1,
	save_latent = integer(0), inner = control_inner())
{
	ret = list(R = R, burn = burn, thin = thin, report = report,
		save_latent = save_latent, inner = inner)
	class(ret) = "control_gibbs"
	return(ret)
}

#' fixed_gibbs
#'
#' Specify components as fixed in the Gibbs sampler so that they will be kept
#' at their initial value.
#'
#' @param beta Logical; if `TRUE`, fix \eqn{\beta} at initial value.
#' @param sigma2 Logical; if `TRUE`, fix \eqn{\sigma^2} at initial value.
#' @param y Logical; if `TRUE`, fix \eqn{y} at initial value.
#'
#' @returns
#' A list with the data.
#'
#' @export
fixed_gibbs = function(beta = FALSE, sigma2 = FALSE, y = FALSE)
{
	ret = list(beta = beta, sigma2 = sigma2, y = y)
	class(ret) = "fixed_gibbs"
	return(ret)
}

#' prior_gibbs
#'
#' Specify hyperparameters for prior in Gibbs sampler.
#'
#' @param sigma_beta Prior standard deviation for \eqn{\beta}.
#' @param a_sigma Prior shape parameter for \eqn{\sigma^2}.
#' @param b_sigma Prior rate parameter for \eqn{\sigma^2}.
#'
#' @returns
#' A list with the data.
#'
#' @export
prior_gibbs = function(sigma_beta = 10, a_sigma = 5, b_sigma = 5)
{
	ret = list(sigma_beta = sigma_beta, a_sigma = a_sigma, b_sigma = b_sigma)
	class(ret) = "prior_gibbs"
	return(ret)
}

#' init_gibbs
#'
#' Specify initial values for Gibbs sampler.
#'
#' @param n Number of observations.
#' @param d Number of covariates.
#' @param beta A vector of length \eqn{d} or `NULL`.
#' @param sigma2 A positive scalar or `NULL`.
#' @param y A vector of length \eqn{n} or `NULL`.
#'
#' @details
#' When a `NULL` value is specified, we select an initial value internally.
#'
#' @returns
#' A list with the data.
#'
#' @export
init_gibbs = function(n, d, beta = NULL, sigma2 = NULL, y = NULL)
{
	if (is.null(beta)) { beta = numeric(d)	}
	if (is.null(sigma2)) { sigma2 = 1 }
	if (is.null(y)) { y = rep(1, n) }

	stopifnot(length(beta) == d)
	stopifnot(length(sigma2) == 1)
	stopifnot(length(y) == n)

	ret = list(beta = beta, sigma2 = sigma2, y = y)
	class(ret) = "init_gibbs"
	return(ret)
}

#' gibbs
#'
#' Gibbs sampler.
#'
#' @param z Vector of \eqn{n} observations.
#' @param lambda Vector of \eqn{n} standard deviations.
#' @param X An \eqn{n \times d} matrix with covariates.
#' @param init A object from [init_gibbs].
#' @param prior A object from [prior_gibbs].
#' @param control A object from [control_gibbs].
#' @param fixed A object from [fixed_gibbs].
#'
#' @returns
#' A list with draws of \eqn{\beta}, \eqn{\sigma^2}, \eqn{y}, and other results.
#'
#' @export
gibbs = function(z, lambda, X,
	init = init_gibbs(n = length(y), d = ncol(X)),
	prior = prior_gibbs(),
	control = control_gibbs(),
	fixed = fixed_gibbs())
{
	n = length(z)
	d = ncol(X)
	stopifnot(n == length(lambda))
	stopifnot(n == nrow(X))

	XtX = crossprod(X)

	stopifnot(class(control) == "control_gibbs")
	R = control$R
	burn = control$burn
	thin = control$thin
	report = control$report
	save_latent = control$save_latent
	inner_ctrl = control$inner

	## Ensure several conditions for indices:
	## - No duplicate values,
	## - All values are in 1:n.
	## Convert to 0-based indices to pass to C++
	stopifnot(table(save_latent) == 1)
	stopifnot(all(save_latent %in% 1:n))

	inner_method = inner_ctrl$method
	max_rejects = inner_ctrl$max_rejects
	tol_suff = inner_ctrl$tol_suff
	tol_merge = inner_ctrl$tol_merge
	tune = inner_ctrl$tune
	N = inner_ctrl$N

	rep_keep = 0
	R_keep = ceiling((R - burn) / thin)

	# Set up histories
	beta_hist = matrix(NA, R_keep, d)
	sigma2_hist = numeric(R_keep)
	y_hist = matrix(NA, R_keep, length(save_latent))

	# Some statistics on tuning
	rejects_hist = integer(R)
	comps_hist = integer(R)
	tunes_hist = integer(R)
	tuned_hist = integer(R)

	# Set up initial values
	stopifnot(class(init) == "init_gibbs")
	beta = init$beta
	y = init$y
	sigma2 = init$sigma2

	Xbeta = X %*% beta

	# Initialize self-tuned VWS proposals
	proposals = list()
	for (i in 1:n) {
	 	proposals[[i]] = ln_norm_proposal$new(z[i], lambda[i], Xbeta[i], sqrt(sigma2))
	}

	# Set up fixed parameters
	stopifnot(class(fixed) == "fixed_gibbs")

	# Set up prior
	stopifnot(class(prior) == "prior_gibbs")
	sigma_beta = prior$sigma_beta
	a_sigma = prior$a_sigma
	b_sigma = prior$b_sigma
	sigma2_beta = sigma_beta*sigma_beta

	# Set up timers
	elapsed = list(beta = 0, sigma2 = 0, y = 0)

	for (rep in 1:R)
	{
		# Draw [y | rest]
		if (!fixed$y)
		{
			st = Sys.time()

			if (inner_method == "vws-tune")
			{
				# Self-tuned VWS
				for (i in 1:n) {
					proposals[[i]]$update(Xbeta[i], sqrt(sigma2))

					if (rep < tune) {
						vws_out = proposals[[i]]$draw_tune(tol_suff, tol_merge, max_rejects)
						y[i] = vws_out$draws
						rejects_hist[rep] = rejects_hist[rep] + vws_out$rejects
						tunes_hist[rep] = tunes_hist[rep] + vws_out$tunes
						tuned_hist[rep] = tuned_hist[rep] + (vws_out$tunes > 0)
					} else {
						vws_out = proposals[[i]]$draw(max_rejects)
						y[i] = vws_out$draws
						rejects_hist[rep] = rejects_hist[rep] + vws_out$rejects
					}
				}
			}
			else if (inner_method == "vws-basic")
			{
				# VWS without tuning
				for (i in 1:n) {
					vws_out = draw_ln_norm(z[i], lambda[i], Xbeta[i],
						sqrt(sigma2), tol_suff, N, max_rejects)
					y[i] = vws_out$draws
					rejects_hist[rep] = rejects_hist[rep] + vws_out$rejects
                }
            }
			else
            {
                stop("Unrecognized method in inner_ctrl")
            }

			elapsed$y = elapsed$y + as.numeric(Sys.time() - st, units = "secs")
		}

		# Draw [beta | rest]
		if (!fixed$beta) {
			st = Sys.time()
			Omega = (1 / sigma2) * XtX + (1 / sigma2_beta) * diag(d)
			mm = solve(Omega, crossprod(X, log(y) / sigma2))
			beta = r_mvnorm_prec(mm, Omega)
			Xbeta = X %*% beta
			elapsed$beta = elapsed$beta + as.numeric(Sys.time() - st, units = "secs")
		}

		# Draw [sigma2 | rest]
		if (!fixed$sigma2) {
			st = Sys.time()
			e = log(y) - Xbeta
			aa = a_sigma + n / 2
			bb = b_sigma + 1 / 2 * crossprod(e)
			sigma2 = 1 / rgamma(1, aa, bb)
			elapsed$sigma2 = elapsed$sigma2 + as.numeric(Sys.time() - st, units = "secs")
		}

		# Save total number of mixture components at this point
		comps_hist[rep] = 0
		for (i in 1:n) {
		 	comps_hist[rep] = comps_hist[rep] + proposals[[i]]$size()
		}
		avg_comps = comps_hist[rep] / n

		if (rep > burn && rep %% thin == 0) {
			rep_keep = rep_keep + 1
			beta_hist[rep_keep,] = beta
			sigma2_hist[rep_keep] = sigma2

			for (l in 1:length(save_latent)) {
				i = save_latent[l]
				y_hist[rep_keep, l] = y[i]
			}
		}

		if (rep %% report == 0) {
			if (inner_method == "vws-tune") {
				s = ifelse(rep > report, rep - report, 1)
				idx = seq(s, rep)
				rejects = sum(rejects_hist[idx])
				tunes = sum(tunes_hist[idx])
				vws::logger("[%d] avg-N: %0.4f  tunes: %d  rejects: %d\n", rep,
					avg_comps, tunes, rejects)
			} else {
				s = ifelse(rep > report, rep - report, 1)
				idx = seq(s, rep)
				rejects = sum(rejects_hist[idx])
				vws::logger("[%d] rejects: %d\n", rep, rejects)
			}
		}
	}

	out = list(
		beta = beta_hist,
		sigma2 = sigma2_hist,
		y = y_hist,
		R_keep = R_keep,
		elapsed = elapsed,
		R = R,
		burn = burn,
		thin = thin,
		rejects = rejects_hist,
		comps = comps_hist,
		tunes = tunes_hist,
		tuned = tuned_hist,
		n = n,
		fixed = fixed,
		method = control$inner$method
	)
	class(out) = "gibbs"
	return(out)
}

#' summary
#'
#' S3 method to compute summary from results of [gibbs].
#'
#' @param object A result from [gibbs].
#' @param pr A vector of quantiles to report.
#' @param ... Not used.
#'
#' @returns
#' A data frame.
#'
#' @export
summary.gibbs = function(object, pr = c(0.05, 0.95), ...)
{
	d = ncol(object$beta)

	df_beta = as.data.frame(cbind(
		apply(object$beta, 2, mean),
		apply(object$beta, 2, sd),
		t(apply(object$beta, 2, quantile, probs = pr))
	))
	rownames(df_beta) = sprintf("beta%d", 1:d)

	df_sigma2 = as.data.frame(cbind(
		mean(object$sigma2),
		sd(object$sigma2),
		t(quantile(object$sigma2, probs = pr))
	))
	rownames(df_sigma2) = sprintf("sigma2")

	df = rbind(df_beta, df_sigma2)
	quantile_names = sprintf("%g%%", 100 * pr)
	colnames(df) = c("mean", "sd", quantile_names)
	df = round(df, 4)
	return(df)
}

#' print
#'
#' S3 method to print summary from results of [gibbs].
#'
#' @param x A result from [gibbs].
#' @param pr A vector of quantiles to report.
#' @param ... Not used.
#'
#' @returns None
#' @export
print.gibbs = function(x, pr = c(0.05, 0.95), ...)
{
	cat("Summary of fit\n")
	print(summary(x, pr))

	cat("----\n")
	printf("Iterations: %d  Burn: %d  Thin: %d  Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)

	if (!x$fixed$y && x$method == "vws-tune") {
		cat("----\n")
		printf("y step\n")
		printf("   Proposed: %d  Rejected: %d  Rate: %0.2f%%\n",
			sum(x$rejects) + x$R * x$n, sum(x$rejects),
			100 * sum(x$rejects) / (sum(x$rejects) + x$R * x$n))
		printf("   Avg regions: %g\n", sum(x$comps) / (x$R * x$n))
	} else if (!x$fixed$y && x$method == "vws-basic") {
		cat("----\n")
		printf("y step\n")
		printf("   Proposed: %d  Rejected: %d\n", sum(x$rejects) + x$R * x$n,
			sum(x$rejects))
		printf("   Rejection rate: %g%%\n",
			100 * sum(x$rejects) / (sum(x$rejects) + x$R * x$n))
	}

	cat("----\n")
	printf("Elapsed time (sec):\n")
	x$elapsed$total = sum(unlist(x$elapsed))
	tab = round(as.data.frame(x$elapsed), 4)
	rownames(tab) = ""
	print(tab)
}
