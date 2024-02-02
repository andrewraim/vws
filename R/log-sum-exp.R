#' Log-Sum-Exp
#'
#' Compute \code{log(sum(exp(x)))} but in a more stable way.
#'
#' @param x A numeric vector
#' @param y A numeric vector
#'
#' @details Computed using the method described by user Ben in StackExchange
#' thread \url{https://stats.stackexchange.com/questions/381936/vectorised-computation-of-logsumexp}.
#' A faster C version (requiring a dependency and possibly compilation) is
#' provided in \link[matrixStats]{logSumExp}.
#' 
#' The function \code{log_sub2_exp} expects that each element of \code{x} is
#' larger than or equal to its corresponding element in \code{y}. Otherwise,
#' \code{NaN} will be returned with a warning.
#' 
#' The function \code{log_sub2_exp_signed} can handle inputs where elements of
#' \code{x} are smaller than corresponding elements of \code{y}. The return
#' value \code{modulus} contains \code{log(abs(exp(x) - exp(y)))}, and
#' \code{sign} contains the sign of \code{exp(x) - exp(y)}. Therefore, the
#' difference on the original scale can be reconstituted as
#' \code{sign * exp(modulus)}.
#' 
#' @examples
#' pi = 1:6 / sum(1:6)
#' x = log(2*pi)
#' log(sum(exp(x)))
#' log_sum_exp(x)
#' 
#' #' # Result should be 0
#' x = c(-Inf -Inf, 0)
#' log_sum_exp(x)
#' 
#' # Result should be -Inf
#' x = c(-Inf -Inf, -Inf)
#' log_sum_exp(x)
#' 
#' # Result should be -Inf
#' x = c(-Inf -Inf, Inf)
#' log_sum_exp(x)
#' 
#' # Result should be 5 on the original scale
#' out = log_add2_exp(log(3), log(2))
#' exp(out)
#' 
#' # Result should be 7 on the original scale
#' out = log_sub2_exp(log(12), log(5))
#' exp(out)
#' 
#' # Results should be 7 and -7 on the original scale, respectively
#' out1 = log_sub2_exp_signed(log(12), log(5))
#' out2 = log_sub2_exp_signed(log(5), log(12))
#' out1$sign * exp(out1$modulus)
#' out2$sign * exp(out2$modulus)
#'
#' @name Log-Sum-Exp
NULL

#' @name Log-Sum-Exp
#' @export
log_sum_exp = function(x) {
	k = length(x)
	v = sort(x, decreasing = TRUE)
	s = numeric(k)

	s[1] = v[1]
	for (j in setdiff(seq_len(k), 1)) {
		dd = ifelse(j > 1 && s[j-1] > -Inf, v[j] - s[j-1], v[j])
		s[j] = max(v[j], s[j-1]) + log1p(exp(-abs(dd)))
	}

	return(s[k])
}

#' @name Log-Sum-Exp
#' @export
log_add2_exp = function(x, y)
{
	s = pmin(x,y)
	t = pmax(x,y)
	t + log1p(exp(s - t))
}

#' @name Log-Sum-Exp
#' @export
log_sub2_exp = function(x, y)
{
	x + log1p(-exp(y - x))
}

#' @name Log-Sum-Exp
#' @export
log_sub2_exp_signed = function(x, y)
{
	s = sign(x - y)
	m = log_sub2_exp(ifelse(s >= 0, x, y), ifelse(s >= 0, y, x))
	list(modulus = m, sign = s)
}

