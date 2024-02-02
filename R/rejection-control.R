#' Rejection Control
#'
#' Control object for rejection sampler.
#'
#' @param max_rejects Maximum number of rejections to tolerate before bailing out.
#' @param report Report progress each time this many candidates are proposed.
#' @param extra_outputs If \code{TRUE}, return a list with extended output
#' in addition to the accepted draws. Otherwise only return accepted draws.
#' @param action_incomplete What should happen if sampler halts with
#' \code{max_rejects} rejections: ne of \code{"stop"},  \code{"warning"}, or
#' \code{"message"}.
#'
#' @return
#' A control object to be passed to the \code{rejection} function.
#'
#' @name rejection_control
#' @export
rejection_control = function(max_rejects = 1000, report = 100,
	extra_outputs = FALSE, action_incomplete = c("stop", "warning", "message"))
{
	action_incomplete = switch(match.arg(action_incomplete),
		stop = stop,
		warning = warning,
		message = message)

	out = list(
		max_rejects = max_rejects,
		report = report,
		extra_outputs = extra_outputs,
		action_incomplete = action_incomplete
	)
	structure(out, class = "rejection_control")
}
