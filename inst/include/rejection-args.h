#ifndef VWS_REJECTION_ARGS_H
#define VWS_REJECTION_ARGS_H

/*
* This code follows a specific structure so that we can use the `as` and `wrap`
* constructs to serialize between structs and Rcpp Lists. The structs are
* defined first without Rcpp included yet, then Rcpp is included and the
* implementations for serialization operations are give.
*
* The pattern we follow here is referred to as "intrusive" (rather than
* "non-intrusive") because `wrap` and `as` are defined via member functions.
*
* See the following articles:
* <https://gallery.rcpp.org/articles/custom-templated-wrap-and-as-for-seamingless-interfaces>
* <https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-extending.pdf>
*
* The issue also came up in Stack\ Overflow threads such as:
* <https://stackoverflow.com/questions/51110244/in-rcpp-how-to-get-a-user-defined-structure-from-c-into-r>
* <https://stackoverflow.com/questions/74887786/specialising-rcppas-for-stdarray>
*/

#include <RcppCommon.h>
#include "typedefs.h"
#include "fntl.h"

namespace vws {

/*
* Structure of optional arguments for rejection sampler.
*
*  - `max_rejects`: Maximum number of rejections to tolerate before bailing out.
*
*  - `report_period`: specifies the period in which progress should be reported
*    (printed to the screen as a log message).
*
*  - `max_rejects_action`: what should happen if `max_rejects` rejections have
*    been obtained during sampling. The default action `STOP` results in an
*    exception being thrown; here, any successful draws that may have been
*    obtained are not returned.
*
*  - `log_ratio_ub`: it is possible numerically for log-ratio
*    $\log[f_0(x) / h_0(x)]$ to be greater than zero. This condition should
*    not occur otherwise, and usually indicates a mistake in user code. This
*    argument is the maximum value allowed where an exception will not be
*    thrown.
*/
struct rejection_args
{
	unsigned int max_rejects = std::numeric_limits<unsigned int>::max();
	unsigned int report = std::numeric_limits<unsigned int>::max();
	double ratio_ub = std::exp(1e-5);
	error_action action = error_action::STOP;

	rejection_args() { };
	rejection_args(SEXP obj);
	operator SEXP() const;
};

}

#include <Rcpp.h>

namespace vws {

/*
* Constructor from SEXP objects
*/
inline rejection_args::rejection_args(SEXP obj)
{
	const Rcpp::List& x = Rcpp::as<Rcpp::List>(obj);

	const Rcpp::StringVector& ex_names = { "max_rejects", "report",
		"ratio_ub", "action", "N", "tol"};
	const Rcpp::StringVector& ac_names = x.names();
	const auto& diff = Rcpp::setdiff(ac_names, ex_names);
	if (diff.size() > 0) {
		Rcpp::stop("Unexpected list entries: %s", fntl::paste(diff, ", "));
	}

	if (x.containsElementNamed("action")) {
		unsigned int ac = x["action"];
		action = error_action(ac);
	}

	max_rejects = x.containsElementNamed("max_rejects") ? x["max_rejects"] : max_rejects;
	report = x.containsElementNamed("report") ? x["report"] : report;
	ratio_ub = x.containsElementNamed("ratio_ub") ? x["ratio_ub"] : ratio_ub;
}

/*
* Conversion operators to SEXP objects
*/

inline rejection_args::operator SEXP() const
{
	return Rcpp::List::create(
		Rcpp::Named("max_rejects") = max_rejects,
		Rcpp::Named("report") = report,
		Rcpp::Named("ratio_ub") = ratio_ub,
		Rcpp::Named("action") = fntl::to_underlying(action)
	);
}

}

#endif
