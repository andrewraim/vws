# Overview

We have some examples of VWS sampling already worked out. Let's try to
implement some of them in a vignette. Also we want to give an overview of how
the package works. The examples should also serve as test cases for the
package.

For each example, we can give a bit of background, show the target, show the
decomposition, show or point to the code for functions, demonstrate the sampler
in action, plot the adaptation rate achieved, and plot the result.

Some misc items

- Give figures captions
- Make figures render nicely. This might be doable with a global option.
- See if there's a nice way to let the user determine how optimization works in
  UnivariateConstRegion and Int version. What if they can do it in closed form?
  Maybe make the min and max functions into methods that can be overridden.

# Outline

1. Introduction
2. Preliminaries
	- Review of VWS
	- Lambdas
	- fntl Package
	- The model of coding in C++ and exposing functions to R
	- Explain that computations are kept on log-scale in the package and why we
	  do it.
3. Usage
	a. Overview of API
	b. Quick start in R
		- Quickly recap a sampling problem and constant majorizer approach
		- Give R code to draw a sample
		- Describe important lines of the code.
	c. First example in C++
		- Redo the quick R section with C++.
	d. Custom optimization
		- Mention that we use a numerical method by default
		- R version
		- C++ version
	e. Refining the proposal
		- Our rule of thumb (and greedy option)
		- At specified knots
	f. User-defined regions
4. R API
5. C++ API
	- FMMProposal
	- Regions
		+ Region Base Class
		+ RealConstRegion
		+ IntConstRegion
	- UnivariateHelper
	- Rejection Sampling
	- TypeDefs (including weight_optimization)
	- Utility Functions
		+ optimize_hybrid
		+ Gumbel Distribution
		+ Categorical Distribution
		+ log-sum-exp and related
		+ which
6. Examples
	a. Disclosure avoidance
	b. Bessel
	c. VMF variate generation


**Max and min**: do we want to make this two separate functions, or is that too
much to deal with, having two pass around two objects where one or both might
be missing?

**Review of VWS**: just review enough so that the rest of the vignette is
self-contained. May want to revisit this toward the end of writing.

**Overview of API**: mention both Real and Int versions of ConstRegion. Mention
ability to customize the optimization.

**Review of VWS Method**

- Avoid rehashing everything from the other paper. Just try to make this one
  self-contained.
- Briefly explain constant and linear VWS. 

Candidate C++ code for removal

- [ ] RectConstRegion.h
- [ ] Region.h
- [ ] invgamma.h
- [ ] logit.h
- [ ] rect.h
- [ ] texp.h
- [ ] unif.h
- [ ] which.h

Candidate R code for removal

- [ ] mvnorm.R
- [ ] polar.R

# Older (Possibly Stale) Notes

Overview of Package

- The `FMMProposal` class encapsulates the finite mixture proposal.
- A `Region` contains all of the problem-specific logic.
- `RealConstRegion` implements constant VWS for univariate regions that
  are partitioned into intervals.
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

