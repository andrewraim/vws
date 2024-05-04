// [[Rcpp::depends(vws)]]
#include "vws.h"

// [[Rcpp::export]]
int hello(double x)
{
	printf("Hello world %g\n", x);

	vws::RejectionControl control;
	printf("control.max_rejects_action = %d\n", control.get_max_rejects_action());

	return 0;
}
