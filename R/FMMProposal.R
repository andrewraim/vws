#' fmm_proposal
#'
#' A more friendly constructor for \code{FMMProposal}.
#'
#' @param regions A list of regions that form a partition of the support.
#'
#' @examples
#' # Define base distribution and weight function
#' g = normal_univariate_helper(mean = 0, sd = 5)
#' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
#'
#' # Set up support
#' support = univariate_const_region(-Inf, 10, w, g)
#' regions = support$bifurcate()
#'
#' # Create a finite mixture with one component
#' fmm = fmm_proposal(regions)
#' print(fmm)
#'
#' @export
fmm_proposal = function(regions)
{
	FMMProposal$new(regions)
}

#' FMMProposal
#'
#' An R6 class which represents a VWS proposal: a finite mixture that some
#' specific certain operations.
#'
#' @param regions A list of \code{N} regions that form a partition of the support.
#' @param log_xi_upper The numeric vector
#' \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
#' @param log_xi_lower The numeric vector
#' \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
#' @param bifurcatable A vector of logical values which indicates whether the
#' corresponding regions may be bifurcated further.
#'
#' @details
#' \itemize{
#' \item The list \code{regions} represents
#' (\eqn{\mathcal{D}_1, \ldots, \mathcal{D}_N}).
#' Each entry should be an object based on an R6  subclass of \code{Region}.
#' \item The numeric vectors \code{log_xi_upper} and \code{log_xi_lower}
#' represent
#' \eqn{(\log \overline{\xi}_1, \ldots, \log \overline{\xi}_N)}
#' and
#' \eqn{(\log \underline{\xi}_1, \ldots, \log \underline{\xi}_N)},
#' respectively.
#' \item The logical vector \code{bifurcatable} indicates whether each
#' region can be bifurcated or not.
#' }
FMMProposal = R6::R6Class(

classname = "FMMProposal",
portable = TRUE,
private = list(
	regions = NULL,
	log_xi_upper = NULL,
	log_xi_lower = NULL,
	bifurcatable = NULL
),

public = list(

#' @description
#' Constructor for FMMProposal.
initialize = function(regions)
{
	private$regions = regions
	private$log_xi_upper = Map(function(reg) { reg$xi_upper(log = TRUE) }, regions) |> unlist()
	private$log_xi_lower = Map(function(reg) { reg$xi_lower(log = TRUE) }, regions) |> unlist()
	private$bifurcatable = Map(function(reg) { reg$is_bifurcatable() }, regions) |> unlist()
},

#' @description
#' Access the vector \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
#' @param log If \code{TRUE} compute result on log-scale.
get_xi_upper = function(log = FALSE)
{
	out = private$log_xi_upper
	if (log) return(out) else { return(exp(out)) }
},

#' @description
#' Access the vector \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
#' @param log If \code{TRUE} compute result on log-scale.
get_xi_lower = function(log = FALSE)
{
	out = private$log_xi_lower
	if (log) return(out) else { return(exp(out)) }
},

#' @description
#' Access the vector \code{bifurcatable}.
get_bifurcatable = function()
{
	private$bifurcatable
},

#' @description
#' Access the vector \code{regions}.
get_regions = function()
{
	private$regions
},

#' @description
#' Upper bound for rejection probability.
#' @param byregion If \code{TRUE}, compute bound by region. Otherwise compute
#' total.
#' @param log If \code{TRUE} compute result on log-scale.
#' @return A vector of size \code{N} or a single scalar.
rejection_bound = function(byregion = FALSE, log = FALSE)
{
	# Each region's contribution to the rejection rate bound
	out = log_sub2_exp(private$log_xi_upper, private$log_xi_lower) -
		log_sum_exp(private$log_xi_upper)

	if (!byregion) {
		# Overall rejection rate bound
		out = vws::log_sum_exp(out)
	}
if (any(is.na(out))) browser()
	if (log) { return (out) } else { return(exp(out)) }
},

#' @description
#' Normalizing constant for proposal distribution.
#' @param log If \code{TRUE} compute result on log-scale
nc = function(log = FALSE)
{
	out = vws::log_sum_exp(private$log_xi_upper)
	if (log) { return (out) } else { return(exp(out)) }
},

#' @description
#' Generate draws from proposal distribution.
#' @param n Number of draws to generate.
#' @param indices If \code{TRUE}, return indices of mixture components selected
#' during draws.
#' @return A list which each element is a saved draw.
r = function(n = 1, indices = FALSE)
{
	N = length(private$regions)

	# Draw from the mixing weights, which are given on the log scale and not
	# normalized.
	idx = r_categ(n, p = private$log_xi_upper, log_p = TRUE)

	# Draw the values from the respective mixture components.
	x = list()
	for (i in 1:n) {
		j = idx[i]
		x[[i]] = private$regions[[j]]$r(n = 1) |> unlist()
	}

	if (indices) {
		out = list(x = x, idx = idx)
	} else {
		out = x
	}

	return(out)
},

#' @description
#' Compute density of proposal distribution.
#' @param x A list where each element is a density value to evaluate.
#' @param normalize If \code{TRUE}, apply the normalizing constant; otherwise
#' do not.
#' @param log If \code{TRUE}, return density results on the log-scale.
#' @return A vector or scalar of density values corresponding to \code{x}.
d = function(x, normalize = TRUE, log = FALSE)
{
	stopifnot(is.list(x))
	n = length(x)

	log_nc = 0
	if (normalize) {
		log_nc = self$nc(log = TRUE)
	}

	log_wg = rep(-Inf, n)

	N = length(private$regions)
	for (i in 1:n) {
		# This search could be more efficient, but would need to be done in a
		# way that can support any kind of region. For example, if we can
		# define a "<" operator for region objects, we could consider a binary
		# search.
		for (j in 1:N) {
			reg = private$regions[[j]]
			if (reg$s(x[[i]])) {
				log_wg[i] = reg$w_major(x[[i]], log = TRUE) +
					reg$d_base(x[[i]], log = TRUE)
			}
		}
	}

	out = log_wg - log_nc
	if (log) { return(out) } else { return(exp(out)) }
},

#' @description
#' Compute \eqn{\log w(x) + \log g(x)} for each element of the given \code{x}
#' which is a list whose elements are density values.
#' @param x a list where each element represents one argument.
#' @return A vector or scalar of density values corresponding to \code{x}.
log_target_pdf_unnorm = function(x)
{
	stopifnot(is.list(x))
	reg = private$regions[[1]]
	out = Map(function(z) { reg$w(z, log = TRUE) + reg$d_base(z, log = TRUE) }, x)
	return(unlist(out))
},

#' @description
#' Summary table of regions which compose the proposal distribution. Returns a
#' data frame.
summary = function()
{
	tbl = data.frame(
		Region = Map(function(x) { x$description() }, private$regions) |> unlist(),
		log_xi_upper = Map(function(x) { x$xi_upper(log = TRUE) }, private$regions) |> unlist(),
		log_xi_lower = Map(function(x) { x$xi_lower(log = TRUE) }, private$regions) |> unlist()
	)
	return(tbl)
},

#' @description
#' Display summary table of regions which compose the proposal distribution.
#' Limit the display to \code{n} regions.
#' @param n Number of regions to print.
print = function(n = 5)
{
	N = length(private$regions)
	printf("FMM Proposal with %d regions (display is unsorted)\n", N)

	tbl = self$summary()
	print(head(tbl, n))

	if (N > n) {
		printf("There are %d more regions not displayed\n", N - n)
	}
}

) # Close public
) # Close class

