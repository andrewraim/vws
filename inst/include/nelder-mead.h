#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#include <Rcpp.h>
#include <R_ext/Applic.h>

namespace vws {

typedef std::function<double(const Rcpp::NumericVector&)> mv_function;

class NelderMeadFunctional
{
protected:
	const mv_function& _f;
	double _fnscale;

public:
	NelderMeadFunctional(const mv_function& f, double fnscale)
		: _f(f), _fnscale(fnscale)
	{
	}

	double get_fnscale() const { return _fnscale; }
	const mv_function& get_f() const { return _f; }

	// User should override this with their objective function
	static double eval(int n, double* par, void* ex) {
		const NelderMeadFunctional* p = static_cast<const NelderMeadFunctional*>(ex);
		Rcpp::NumericVector x(par, par + n);
		const mv_function& f = p->get_f();
		return p->get_fnscale() * f(x);
	}
};

// Template Result and Control??!!
// Put fnscale and std::function into ex?
// Put fnscale into Control?

struct NelderMeadResult
{
	Rcpp::NumericVector par;
	double value;
	int fail;
	int fncount;
};

struct NelderMeadControl
{
	double alpha = 1.0;
	double beta = 0.5;
	double gamma = 2.0;
	unsigned int trace = 0;
	double abstol = R_NegInf;
	double intol = sqrt(std::numeric_limits<double>::epsilon());
	unsigned int maxit = 500;
	double fnscale = 1.0;
};

/*
* Nelder-Mead algorithm from R.
* See <https://cran.r-project.org/doc/manuals/R-exts.html#Optimization>
* and <https://stackoverflow.com/questions/12765304/calling-r-function-optim-from-c>
*/
NelderMeadResult nelder_mead(const Rcpp::NumericVector& init,
	mv_function f,
	const NelderMeadControl& control)
{
	NelderMeadResult out;
	unsigned int n = init.length();

	double xin[n];
	for (unsigned int i = 0; i < n; i++) {
		xin[i] = init(i);
	}

	double xout[n];

	NelderMeadFunctional nmf(f, control.fnscale);

	nmmin(
		n,               // In:  int n [number of parameters]
		xin,             // In:  double *xin [initial value]
		xout,            // Out: double *x [point at which optimum is found]
		&out.value,      // Out: double *Fmin [objective value at which optimum is found]
		NelderMeadFunctional::eval, // In:  optimfn fn [objective function]
		&out.fail,       // Out: int *fail [true if the function failed]
		R_NegInf,        // In:  double abstol [absolute tolerance]
		control.intol,   // In:  double intol [user-initialized conversion tolerance]
		&nmf,            // In:  void *ex [external data to pass to the objective function]
		control.alpha,   // In:  double alpha [reflection factor]
		control.beta,    // In:  double beta [contraction and reduction factor]
		control.gamma,   // In:  double gamma [extension factor]
		control.trace,   // In:  int trace [if positive, print progress info]
		&out.fncount,    // Out: int *fncount [number of times the objective function was called]
		control.maxit    // In:  int maxit [maximum number of iterations]
	);

	out.par = Rcpp::NumericVector(xout, xout + n);
	return out;
}

}

#endif
