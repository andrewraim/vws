#' FMM Proposal
#'
#' @param regions A list of regions that form a partition of the support.
#'
#' @field regions A list of regions that form a partition of the support.
#' @field log_xi_upper The numeric vector
#' \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
#' @field log_xi_lower The numeric vector
#' \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
#' @field bifurcatable A vector of logical values which indicates whether the
#' corresponding regions may be bifurcated further.
#'
#' @details
#' The \code{regions} list represents \eqn{\mathcal{D}_1, \ldots, \mathcal{D}_N}.
#' They should be objects based on a subclass of \code{Region}, which uses
#' R's reference class construct.
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
#'     dlnorm(10 - x, 5, 2, log)
#' }
#'
#' # Set up support
#' support = univariate_const_region(-Inf, 10, w, g)
#' regions = support$bifurcate()
#'
#' # Create a finite mixture with one component
#' fmm = fmm_proposal(regions)
#' print(fmm)
#'
#' @import methods
#' @name FMMProposal
NULL

#' FMMProposal
#'
#' @name FMMProposal
#' @export
FMMProposal = setRefClass("FMMProposal",
	fields = c("regions", "log_xi_upper", "log_xi_lower", "bifurcatable"))

#' @name FMMProposal
#' @export
fmm_proposal = function(regions)
{
	fmm = FMMProposal$new()

	fmm$regions = regions
	fmm$log_xi_upper = Map(function(reg) { reg$xi_upper(log = TRUE) }, regions) |> unlist()
	fmm$log_xi_lower = Map(function(reg) { reg$xi_lower(log = TRUE) }, regions) |> unlist()
	fmm$bifurcatable = Map(function(reg) { reg$is_bifurcatable() }, regions) |> unlist()

	return(fmm)
}

FMMProposal$methods(rejection_bound = function(byregion = FALSE, log = FALSE)
{
	'
	Upper bound for rejection probability.
	\\subsection{Arguments}{\\itemize{
		\\item{\\code{byregion} If \\code{TRUE}, compute bound by region.
			Otherwise compute total.}
		\\item{\\code{log} If \\code{TRUE} compute result on log-scale.}
	}}
	\\subsection{Value}{A vector of size \\code{N} or a single scalar.}
	'

	# Each region's contribution to the rejection rate bound
	out = log_sub2_exp(log_xi_upper, log_xi_lower) - log_sum_exp(log_xi_upper)

	if (!byregion) {
		# Overall rejection rate bound
		out = vws::log_sum_exp(out)
	}

	if (log) { return (out) } else { return(exp(out)) }
})

FMMProposal$methods(nc = function(log = FALSE)
{
	'
	Normalizing constant for proposal distribution.
	\\subsection{Arguments}{\\itemize{
		\\item{\\code{log} If \\code{TRUE} compute result on log-scale.}
	}}
	\\subsection{Value}{A scalar.}
	'

	out = vws::log_sum_exp(log_xi_upper)
	if (log) { return (out) } else { return(exp(out)) }
})

FMMProposal$methods(r = function(n = 1, indices = FALSE)
{
	'
	Generate draws from proposal distribution.
	\\subsection{Arguments}{\\itemize{
		\\item{\\code{n} Number of draws to generate.}
		\\item{\\code{indices} If \\code{TRUE}, return indices of mixture
			components selected during draws.}
	}}
	\\subsection{Value}{A list which each element is a saved draw.}
	'

	N = length(regions)

	# Draw from the mixing weights, which are given on the log scale and not
	# normalized.
	idx = r_categ(n, p = log_xi_upper, log_p = TRUE)

	# Draw the values from the respective mixture components.
	x = list()
	for (i in 1:n) {
		j = idx[i]
		x[[i]] = regions[[j]]$r(n = 1) |> unlist()
	}

	if (indices) {
		out = list(x = x, idx = idx)
	} else {
		out = x
	}

	return(out)
})

FMMProposal$methods(d = function(x, normalize = TRUE, log = FALSE)
{
	'
	Compute density of proposal distribution.
	\\subsection{Arguments}{\\itemize{
		\\item{\\code{x} A list where each element is a density value to
			evaluate.}
		\\item{\\code{normalize} If \\code{TRUE}, apply the normalizing
			constant; otherwise do not.}
		\\item{\\code{log} If \\code{TRUE}, return density results on the
			log-scale.}
	}}
	\\subsection{Value}{A vector of density values corresponing to \\code{x}.}
	'

	stopifnot(is.list(x))
	n = length(x)

	log_nc = 0
	if (normalize) {
		log_nc = nc(log = TRUE)
	}

	log_wg = rep(-Inf, n)

	N = length(regions)
	for (i in 1:n) {
		for (j in 1:N) {
			reg = regions[[j]]
			if (reg$in_support(x[[i]])) {
				log_wg[i] = reg$w_major(x[[i]], log = TRUE) +
					reg$d_base(x[[i]], log = TRUE)
			}
		}
	}

	out = log_wg - log_nc
	if (log) { return(out) } else { return(exp(out)) }
})

# The Region class knows how to compute the unnormalized target pdf.
FMMProposal$methods(log_target_pdf_unnorm = function(x)
{
	'
	Compute \\eqn{\\log w(x) + \\log g(x)} for each element of the given
	\\code{x} which is a list whose elements are density values.
	'

	if (!is.list(x)) { x = Map(identity, x) }
	reg = regions[[1]]
	out = Map(function(x) { reg$w(x, log = TRUE) + reg$d_base(x, log = TRUE) }, x)
	return(unlist(out))
})

FMMProposal$methods(summary = function()
{
	'
	Summary table of regions which compose the proposal distribution. Returns a
	data frame.
	'

	tbl = data.frame(
		Region = Map(function(x) { x$description() }, regions) |> unlist(),
		log_xi_upper = Map(function(x) { x$xi_upper(log = TRUE) }, regions) |> unlist(),
		log_xi_lower = Map(function(x) { x$xi_lower(log = TRUE) }, regions) |> unlist()
	)
	return(tbl)
})

FMMProposal$methods(show = function(n = 5)
{
	'
	Display summary table of regions which compose the proposal distribution.
	Limit the display to \\code{n} regions.
	'

	N = length(regions)
	printf("FMM Proposal with %d regions (display is unsorted)\n", N)

	tbl = summary()
	print(head(tbl, n))

	if (N > n) {
		printf("There are %d more regions not displayed\n", N - n)
	}
})
