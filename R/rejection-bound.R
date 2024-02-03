#' Bound for probability of rejection
#'
#' Compute the expression
#' \deqn{
#' 1 - \sum_{j=1}^N \underline{\xi}_j \Bigg/ \sum_{j=1}^N \overline{\xi}_j
#' }
#' in a numerically stable way.
#'
#' @param log_xi_upper The vector \eqn{\log \overline{\xi}_1, \ldots, \log \overline{\xi}_N}
#' @param log_xi_lower The vector \eqn{\log \underline{\xi}_1, \ldots, \log \underline{\xi}_N}
#' @param log If \code{TRUE}, return result on the log-scale
#'
#' @export
rejection_bound = function(log_xi_upper, log_xi_lower, log = FALSE) {
	# Do the following calculation, but in a more numerically stable way.
	# out = 1 - sum(exp(log_xi_lower)) / sum(exp(log_xi_upper))

	log_num = log_sum_exp(log_xi_lower)
	log_den = log_sum_exp(log_xi_upper)
	out = log_sub2_exp(0, log_num - log_den)
	if (log) { return(out) } else { return(exp(out)) }
}

