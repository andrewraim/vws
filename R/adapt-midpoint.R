#' Adapt
#'
#' Adapt an FMM proposal using a midpoint rule-of-thumb.
#'
#' @param h An FMM proposal
#' @param N Number of additional mixture components after adaptation.
#' @param control A control object from \code{rejection_control}
#' @param method Either \code{greedy} or \code{pps}. The \code{greedy} method
#' always chooses to bifrucate the region with the largest contribution to the
#' probability of rejection. The \code{pps} method draws randomly, with
#' probability proportion to that contribution.
#'
#' @examples
#' #  Define base distribution
#' g = univariate_helper(
#'     r = function(n) rnorm(n, 0, 5),
#'     d = function(x, log = FALSE) dnorm(x, 0, 5, log),
#'     p = function(q, lower.tail = TRUE, log.p = FALSE) {
#'         pnorm(q, 0, 5, lower.tail, log.p)
#'     },
#'     q = function(p, lower.tail = TRUE, log.p = FALSE) {
#'         qnorm(p, 0, 5, lower.tail, log.p)
#'     },
#'     in_support = function(x) { TRUE }
#' )
#'
#' # Define weight function
#' w = function(x, log = FALSE) {
#'     dlnorm(10 - x, meanlog = 5, sdlog = 2, log)
#' }
#'
#' # Set up support
#' support = univariate_const_region(-Inf, 10, w, g)
#' regions = list(support)
#'
#' # Adapt the proposal
#' h_init = fmm_proposal(regions)
#' adapt_out = adapt_midpoint(h_init, N = 100, control = rejection_control(report = 1))
#'
#' # Create a finite mixture proposal
#' h = adapt_out$h
#' h$rejection_bound()
#' h$rejection_bound(byregion = TRUE)
#'
#' out = rejection(h, n = 1000)
#' print(out |> unlist())
#'
#' out = rejection(h, n = 1000, rejection_control(extra_outputs = TRUE))
#' print(out$draws)
#' print(out$rejects)
#'
#' @name adapt
#' @export
adapt_midpoint = function(h, N, method = "pps", control = rejection_control())
{
	log_ub_hist = numeric(N+1)
	log_lb_hist = numeric(N+1)
	log_bdd_hist = numeric(N+1)
	log_ub_hist[1] = log_sum_exp(h$log_xi_upper)
	log_lb_hist[1] = log_sum_exp(h$log_xi_lower)
	log_bdd_hist[1] = h$rejection_bound(log = TRUE)

	report = control$report

	if (report <= N) {
		logger("After 0 steps, attained a = %g ", log_ub_hist[1])
		printf("rejection prob <= %g\n", exp(log_bdd_hist[1]))
	}

	for (j in seq_len(N)) {
		# Recall that region volumes reflect where there mixture is further
		# from the target: it takes into account both the weight difference
		# and the probability of being in that region.
		L = length(h$regions)
		log_xi_upper = h$log_xi_upper
		log_xi_lower = h$log_xi_lower

		# Each region's contribution to the rejection rate
		log_volume = h$rejection_bound(byregion = TRUE, log = TRUE)

		bifurcatable = h$bifurcatable
		idx = which(bifurcatable)
		if (length(idx) == 0) {
			warning("No regions left to bifurcate")
			break
		}

		if (method == "greedy") {
			jdx = which.max(log_volume[idx])
		} else if (method == "pps") {
			jdx = r_categ(n = 1, log_volume[idx], log_p = TRUE)
		} else {
			stopifnot("method must be greedy or pps")
		}
		reg = h$regions[[idx[jdx]]]

		# Split the target region and make another proposal with it
		bif_out = reg$bifurcate()
		regions_new = append(h$regions[-idx[jdx]], bif_out)
		h = fmm_proposal(regions_new)

		log_ub_hist[j+1] = log_sum_exp(h$log_xi_upper)
		log_lb_hist[j+1] = log_sum_exp(h$log_xi_lower)
		log_bdd_hist[j+1] = h$rejection_bound(log = TRUE)

		if (j %% report == 0) {
			logger("After %d steps, attained a = %g ", j, log_ub_hist[j+1])
			printf("rejection prob <= %g\n", exp(log_bdd_hist[j+1]))
		}
	}

	list(h = h, log_ub_hist = log_ub_hist, log_lb_hist = log_lb_hist,
		log_bdd_hist = log_bdd_hist)
}
