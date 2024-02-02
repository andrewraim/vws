#' FMM Proposal
#'
#' @param regions A list of regions that form a partition of the support
#' @param x An FMM Proposal to operate on
#' @param object An FMM Proposal to operate on
#' @param n Number of regions to
#' @param ... Additional arguments
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
#' @name FMMProposal
NULL

#' FMMProposal
#'
#' Use reference classes here to allow proposal to be mutable.
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

FMMProposal$methods(rejection_bound = function(log = FALSE)
{
	vws::rejection_bound(.self$log_xi_upper, .self$log_xi_lower, log = log)
})

FMMProposal$methods(nc = function(log = FALSE)
{
	log_nc = vws::log_sum_exp(.self$log_xi_upper)
	if (log) { return (log_nc) } else { return(exp(log_nc)) }
})

FMMProposal$methods(r = function(n = 1, indices = FALSE)
{
	N = length(.self$regions)

	# Draw from the mixing weights, which are given on the log scale and not
	# normalized.
	idx = r_categ(n, p = .self$log_xi_upper, log_p = TRUE)

	# Draw the values from the respective mixture components.
	x = list()
	for (i in 1:n) {
		j = idx[i]
		x[[i]] = .self$regions[[j]]$r(n = 1) |> unlist()
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
	stopifnot(is.list(x))
	n = length(x)

	log_nc = 0
	if (normalize) {
		log_nc = .self$nc(log = TRUE)
	}

	log_wg = rep(-Inf, n)

	N = length(.self$regions)
	for (i in 1:n) {
		for (j in 1:N) {
			reg = .self$regions[[j]]
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
	if (!is.list(x)) { x = Map(identity, x) }
	reg = .self$regions[[1]]
	out = Map(function(x) { reg$w(x, log = TRUE) + reg$d_base(x, log = TRUE) }, x)
	return(unlist(out))
})

FMMProposal$methods(summary = function()
{
	tbl = data.frame(
		Region = Map(function(x) { x$description() }, .self$regions) |> unlist(),
		log_xi_upper = Map(function(x) { x$xi_upper(log = TRUE) }, .self$regions) |> unlist(),
		log_xi_lower = Map(function(x) { x$xi_lower(log = TRUE) }, .self$regions) |> unlist()
	)
	return(tbl)
})

FMMProposal$methods(show = function(n = 5)
{
	N = length(.self$regions)
	printf("FMM Proposal with %d regions (display is unsorted)\n", N)

	tbl = .self$summary()
	print(head(tbl, n))

	if (N > n) {
		printf("There are %d more regions not displayed\n", N - n)
	}
})
