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
#' @name adapt
#' @export
adapt_midpoint = function(h, N, control = rejection_control(), method = "greedy")
{
	log_ub_hist = numeric(N+1)
	log_lb_hist = numeric(N+1)
	log_bdd_hist = numeric(N+1)
	log_ub_hist[1] = log_sum_exp(h$log_xi_upper)
	log_lb_hist[1] = log_sum_exp(h$log_xi_lower)
	log_bdd_hist[1] = rejection_bound(h$log_xi_upper, h$log_xi_lower, log = TRUE)

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
		log_volume = log_sub2_exp(log_xi_upper, log_xi_lower) -
			log_sum_exp(log_xi_upper)

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
		bif_out = bifurcate(reg)
		regions_new = append(h$regions[-idx[jdx]], bif_out)
		h = get_fmm_proposal(regions_new)

		log_ub_hist[j+1] = log_sum_exp(h$log_xi_upper)
		log_lb_hist[j+1] = log_sum_exp(h$log_xi_lower)
		log_bdd_hist[j+1] = rejection_bound(h$log_xi_upper, h$log_xi_lower, log = TRUE)

		if (j %% report == 0) {
			logger("After %d steps, attained a = %g ",
				j, log_ub_hist[j+1])
			printf("rejection prob <= %g\n", exp(log_bdd_hist[j+1]))
		}
	}

	list(h = h, log_ub_hist = log_ub_hist, log_lb_hist = log_lb_hist,
		log_bdd_hist = log_bdd_hist)
}
