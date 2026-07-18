#ifndef VWS_H
#define VWS_H

/*
* Export classes, functions, etc to package users.
*
* Note that `result.h` must be listed first because it uses a technique
* described in [1] to serialize objects between R and C++.
*
* [1] Dirk Eddelbuettel and Romain Francois. Rcpp Extending. Vignette from
* Rcpp version 1.0.13.
*/

#include "result.h"

#include "categ.h"
#include "fmm-proposal.h"
#include "gumbel.h"
#include "int-const-region.h"
#include "log-sum-exp.h"
#include "logit.h"
#include "optimize-hybrid.h"
#include "region.h"
#include "real-const-region.h"
#include "real-const-region-defaults.h"
#include "rect.h"
#include "rejection.h"
#include "rejection-args.h"
#include "rejection-tune.h"
#include "seq.h"
#include "typedefs.h"
#include "univariate-helper.h"

/*
* Aliases for original camel case class names.
*/

namespace vws {
	template<class T, class R>
	using FMMProposal = fmm_proposal<T,R>;

	template<class T>
	using Region = region<T>;

	using RealConstRegion = real_const_region;
	using IntConstRegion = int_const_region;
	using UnivariateHelper = univariate_helper;
}

#endif

