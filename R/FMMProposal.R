#' FMMProposal
#'
#' An R6 class which represents a VWS proposal: a finite mixture that some
#' specific certain operations.
#'
#' @details
#' \itemize{
#' \item The list `regions` represents
#' (\eqn{\mathcal{D}_1, \ldots, \mathcal{D}_N}).
#' Each entry should be an object based on an R6  subclass of `Region`.
#' \item The numeric vectors `log_xi_upper` and `log_xi_lower`
#' represent
#' \eqn{(\log \overline{\xi}_1, \ldots, \log \overline{\xi}_N)}
#' and
#' \eqn{(\log \underline{\xi}_1, \ldots, \log \underline{\xi}_N)},
#' respectively.
#' \item The logical vector `bifurcatable` indicates whether each
#' region can be bifurcated or not.
#' }
#'
#' @export
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
#' @param regions A list of objects whose class derives from `Region`.
initialize = function(regions)
{
	private$regions = regions
	private$log_xi_upper = Map(function(reg) { reg$xi_upper(log = TRUE) }, regions) |> unlist()
	private$log_xi_lower = Map(function(reg) { reg$xi_lower(log = TRUE) }, regions) |> unlist()
	private$bifurcatable = Map(function(reg) { reg$is_bifurcatable() }, regions) |> unlist()
},

#' @description
#' Access the vector \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
#' @param log If `TRUE` compute result on log-scale.
get_xi_upper = function(log = FALSE)
{
	out = private$log_xi_upper
	if (log) return(out) else { return(exp(out)) }
},

#' @description
#' Access the vector \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
#' @param log If `TRUE` compute result on log-scale.
get_xi_lower = function(log = FALSE)
{
	out = private$log_xi_lower
	if (log) return(out) else { return(exp(out)) }
},

#' @description
#' Access the vector `bifurcatable`.
get_bifurcatable = function()
{
	private$bifurcatable
},

#' @description
#' Access the vector `regions`.
get_regions = function()
{
	private$regions
},

#' @description
#' Upper bound for rejection probability.
#' @param log If `TRUE` compute result on log-scale.
#' @return A single scalar.
rejection_bound = function(log = FALSE)
{
	# Each region's contribution to the rejection rate bound

    lxl = private$log_xi_lower
    lxu = private$log_xi_upper
	log_bound = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu)
    out = log_sum_exp(log_bound)
	if (log) { return (out) } else { return(exp(out)) }
},

#' @description
#' Upper bound for rejection probability; contribution per region.
#' @param log If `TRUE` compute result on log-scale.
#' @return A vector of size `N`.

rejection_bound_regions = function(log = FALSE)
{
    # Each region's contribution to the rejection rate bound.
    lxl = private$log_xi_lower
    lxu = private$log_xi_upper
    out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu)
	if (log) { return (out) } else { return(exp(out)) }
},

#' @description
#' Normalizing constant for proposal distribution.
#' @param log If `TRUE` compute result on log-scale
nc = function(log = FALSE)
{
	out = vws::log_sum_exp(private$log_xi_upper)
	if (log) { return (out) } else { return(exp(out)) }
},

#' @description
#' Generate draws from proposal distribution.
#' @param n Number of draws to generate.
#' @param indices If `TRUE`, return indices of mixture components selected
#' during draws.
#' @return A list which each element is a saved draw.
r = function(n = 1, indices = FALSE)
{
	N = length(private$regions)

	# Draw from the mixing weights, which are given on the log scale and not
	# normalized.
	idx = r_categ(n, p = private$log_xi_upper, log = TRUE,one_based = TRUE)

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
#' @param normalize If `TRUE`, apply the normalizing constant; otherwise
#' do not.
#' @param log If `TRUE`, return density results on the log-scale.
#' @return A vector or scalar of density values corresponding to `x`.
d = function(x, normalize = TRUE, log = FALSE)
{
	log_nc = 0
	if (normalize) {
		log_nc = self$nc(log = TRUE)
	}

	log_wg = -Inf

	N = length(private$regions)

	# This search could be more efficient, but would need to be done in a
	# way that can support any kind of region. For example, if we can
	# define a "<" operator for region objects, we could consider a binary
	# search.
	for (j in 1:N) {
		reg = private$regions[[j]]
		if (reg$s(x)) {
			log_wg = reg$w_major(x, log = TRUE) + reg$d_base(x, log = TRUE)
		}
	}

	out = log_wg - log_nc
	if (log) { return(out) } else { return(exp(out)) }
},

#' @description
#' Compute \eqn{f_0(x) = w(x) g(x)} for each element of the given `x`
#' which is a list whose elements are density values.
#' @param x a list where each element represents one argument.
#' @param log If `TRUE`, return density results on the log-scale.
#' @return A vector or scalar of density values corresponding to `x`.
d_target_unnorm = function(x, log = TRUE)
{
	reg = private$regions[[1]]
	out = reg$w(x, log = TRUE) + reg$d_base(x, log = TRUE)
	if (log) { return(unlist(out)) } else { return(exp(unlist(out))) }
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
#' Limit the display to `n` regions.
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

