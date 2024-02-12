#' Vertical Weighted Strips Rejection Sampler
#'
#' Accept-reject algorithm using our proposal for weighted distributions.
#'
#' @param h An \code{fmm_proposal} object
#' @param n Number of desired draws
#' @param control A control object from \code{rejection_control}
#'
#' @return
#' A list whose structure depends on
#' \code{extra_outputs}. If \code{extra_outputs = FALSE}, the list is of length
#' \code{n} where each element represents one draw. If
#' \code{extra_outputs = TRUE}, the list contains the following named elements:
#'
#' \item{draws}{is list is of length \code{n} where each element represents one
#' draw.}
#' \item{rejects}{is a vector of counts. The \eqn{i}th element is the number of
#' rejections before the \eqn{i}th successful draw.}
#' \item{h}{returns the \code{h} that was input to the sampler.}
#'
#' Note that \code{h} was originally intended to show adaptation during
#' sampling, but that is currently not done.
#'
#' @examples
#' # Define base distribution and weight function
#' g = normal_univariate_helper(0, 5)
#' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
#'
#' # Set up support
#' support = univariate_const_region(-Inf, 10, w, g)
#' regions = support$bifurcate()
#'
#' # Create a finite mixture proposal
#' h = fmm_proposal(regions)
#' h$rejection_bound()
#' h$rejection_bound(byregion = TRUE)
#'
#' out = rejection(h, n = 1000)
#' print(out |> unlist())
#'
#' out = rejection(h, n = 1000, rejection_control(extra_outputs = TRUE))
#' print(out$draws |> unlist())
#' print(out$rejects)
#'
#' @name rejection
#' @export
rejection = function(h, n = 1, control = rejection_control())
{
	out = list()
	N_rejects = 0
	rejects = integer(n)
	accept = FALSE

	stopifnot(is(h, "FMMProposal"))

	stopifnot("rejection_control" %in% class(control))
	max_rejects = control$max_rejects
	report = control$report
	extra_outputs = control$extra_outputs
	action_incomplete = control$action_incomplete

	# The constant M in the acceptance ratio is always M = 1.
	log_M = 0

	for (i in 1:n) {
		accept = FALSE
		while (!accept && N_rejects < max_rejects) {
			v = runif(1)
			x = h$r(n = 1)
			log_fx = h$log_target_pdf_unnorm(x)
			log_hx = h$d(x, normalize = FALSE, log = TRUE)
			log_ratio = log_fx - log_hx - log_M

			if (log(v) < log_ratio) {
				# Accept x as a draw from f(x)
				out[[i]] = x |> unlist()
				accept = TRUE
			} else {
				# Reject x and adapt the proposal
				N_rejects = N_rejects + 1
				rejects[i] = rejects[i] + 1
			}

			# Report progress after `report` candidates
			N_accepts = i - 1 + accept
			if ((N_rejects + N_accepts) %% report == 0) {
				logger("After %d candidates, %d accepts and %d rejects\n",
					N_accepts + N_rejects, N_accepts, N_rejects)
			}
		}
	}

	if (N_rejects == max_rejects) {
		msg = sprintf("Reached maximum number of rejects: %d", max_rejects)
		action_incomplete(msg)
	}

	if (extra_outputs) {
		res = list(draws = out, rejects = rejects, h = h)
	} else {
		res = out
	}

	return(res)
}
