#' Rejection Control
#'
#' Control object for rejection sampler.
#'
#' @param N Maximum number of regions to create in the proposal.
#' @param tol Desired bound for probability of rejection in the proposal.
#' @param max_rejects Maximum number of rejections to tolerate before bailing out.
#' @param report Report progress each time this many candidates are proposed.
#' @param ratio_ub TBD
#' @param action What should happen if sampler halts with `max_rejects`
#' rejections: one of `"stop"`,  `"warning"`, or `"message"`.
#'
#' @return
#' A control object to be passed to the `rejection` function.
#'
#' @name rejection_control
#' @export
rejection_control = function(N = 30, tol = 0.1,
	max_rejects = .Machine$integer.max,
	report = .Machine$integer.max, ratio_ub = exp(1e-5),
	action = c("stop", "warning", "message", "none"))
{
	action = switch(match.arg(action),
		stop = 0L,
		warning = 1L,
		message = 2L,
		none = 3L
	)

	out = list(
		N = N,
		tol = tol,
		max_rejects = max_rejects,
		report = report,
		ratio_ub = ratio_ub,
		action = action
	)

	structure(out, class = "rejection_control")
}
