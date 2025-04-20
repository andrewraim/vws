#ifndef VWS_H
#define VWS_H

/*
* Export classes, functions, etc to package users.
*
* Note that `result.h` must be listed first because it uses a technique
* described in [1].
*
* [1] Dirk Eddelbuettel and Romain Francois. Rcpp Extending. Vignette from
* Rcpp version 1.0.13.
*/

#include "result.h"

#include "categ.h"
#include "FMMProposal.h"
#include "gamma-trunc.h"
#include "gumbel.h"
#include "IntConstRegion.h"
#include "invgamma.h"
#include "log-sum-exp.h"
#include "logger.h"
#include "logit.h"
#include "optimize-hybrid.h"
#include "Region.h"
#include "RealConstRegion.h"
#include "rect.h"
#include "rejection.h"
#include "rejection-args.h"
#include "texp.h"
#include "typedefs.h"
#include "UnivariateHelper.h"
#include "unif.h"

#endif
