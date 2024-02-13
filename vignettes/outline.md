# Overview

We have some examples of VWS sampling already worked out. Let's try to
implement some of them in a vignette. Also we want to give an overview of how
the package works. The examples should also serve as test cases for the
package.

For each example, we can give a bit of background, show the target, show the
decomposition, show or point to the code for functions, demonstrate the sampler
in action, plot the adaptation rate achieved, and plot the result.

# Outline

Review of VWS Method

- Avoid rehashing everything from the other paper. Just try to make this one
  self-contained.
- Briefly explain constant and linear VWS. 

Overview of Package

- We use R6 in some places to promote more formal object-orientation.
- We don't use Rcpp at the moment because interoperability between C++ classes
  and R becomes more complicated. It is possible with Rcpp Modules though.
- The `FMMProposal` class encapsulates the finite mixture proposal.
- A `Region` contains all of the problem-specific logic.
- `UnivariateConstRegion` implements constant VWS for univariate regions that
  are partitioned into intervals. It has a built in optimization (may want to
  make this more configurable!)
- `IntUnivariateConstRegion` is a special version for integer supports. E.g.,
  it avoids bifurcating into regions that contain no integers.
- `univariate_helper` for `UnivariateConstRegion` and `IntUnivariateConstRegion`.
- The `adapt_midpoint` function partitions the proposal to make it better
  approximate the target.
- The `rejection` function does rejection sampling with a proposal.
- Maybe explain all of the pieces first, then show how to use them in more
  detail.
- Remark: many of the calculations are done on the log-scale. Describe the
  log-sum-exp functions.
- Remark: Gumbel trick to draw from the proposal. (Does this preclude doing a
  binary search though?)

Possible Examples

- Lognormal-Normal
- Conway-Maxwell Poisson
- Bessel count distribution
- Skew Normal
- VMF: the nontrivial marginal
	+ Const
	+ Linear
- VMF: the posterior in iid Bayesian analysis
- Gaussian process regression

