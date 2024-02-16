#' Adapt Using Midpoint
#'
#' Adapt an FMM proposal using a midpoint rule-of-thumb.
#'
#' @param h An object of class \code{FMMProposal}.
#' @param N Number of additional mixture components after adaptation.
#' @param report Report progress each time this many candidates are proposed.
#'
#' @examples
#' g = normal_univariate_helper(0, 5)
#' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
#' support = UnivariateConstRegion$new(-Inf, 10, w, g)
#' regions = list(support)
#'
#' # Adapt the proposal
#' h_init = FMMProposal$new(regions)
#' adapt_out = adapt(h_init, N = 100)
#' h = adapt_out$h
#'
#' # Create a finite mixture proposal
#' h$rejection_bound()
#' h$rejection_bound(byregion = TRUE)
#'
#' out = rejection(h, n = 200)
#' print(out |> unlist())
#'
#' out = rejection(h, n = 200, rejection_control(extra_outputs = TRUE))
#' print(out$draws |> unlist())
#' print(out$rejects)
#'
#' @name adapt
#' @export
adapt = function(h, N, report = N+1)
{
	log_ub_hist = numeric(N+1)
	log_lb_hist = numeric(N+1)
	log_bdd_hist = numeric(N+1)
	log_ub_hist[1] = log_sum_exp(h$get_xi_upper(log = TRUE))
	log_lb_hist[1] = log_sum_exp(h$get_xi_lower(log = TRUE))
	log_bdd_hist[1] = h$rejection_bound(log = TRUE)

	if (report <= N) {
		logger("Initial log Pr{rejection} <= %g\n", log_bdd_hist[1])
	}

	for (j in seq_len(N)) {
		# Recall that region volumes reflect where there mixture is further
		# from the target: it takes into account both the weight difference
		# and the probability of being in that region.
		L = length(h$get_regions())
		log_xi_upper = h$get_xi_upper(log = TRUE)
		log_xi_lower = h$get_xi_lower(log = TRUE)

		# Each region's contribution to the rejection rate
		log_volume = h$rejection_bound(byregion = TRUE, log = TRUE)

		idx = which(h$get_bifurcatable())
		if (length(idx) == 0) {
			warning("No regions left to bifurcate")
			break
		}

		jdx = r_categ(n = 1, log_volume[idx], log_p = TRUE)
		reg = h$get_regions()[[idx[jdx]]]

		# Split the target region and make another proposal with it
		bif_out = reg$bifurcate()
		regions_new = append(h$get_regions()[-idx[jdx]], bif_out)
		h = FMMProposal$new(regions_new)

		log_ub_hist[j+1] = log_sum_exp(h$get_xi_upper(log = TRUE))
		log_lb_hist[j+1] = log_sum_exp(h$get_xi_lower(log = TRUE))
		log_bdd_hist[j+1] = h$rejection_bound(log = TRUE)

		if (j %% report == 0) {
			logger("After %d steps log Pr{rejection} <= %g\n", j, log_bdd_hist[j+1])
		}
	}

	list(h = h, log_ub_hist = log_ub_hist, log_lb_hist = log_lb_hist,
		log_bdd_hist = log_bdd_hist)
}
