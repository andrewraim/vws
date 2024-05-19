CustomConstRegion = R6::R6Class(

classname = "CustomConstRegion",
inherit = UnivariateConstRegion,

public = list(

initialize = function(a, b, w, g) {
	super$initialize(a, b, w, g)
},

optimize = function(maximize = TRUE, log = TRUE) {
	a = private$a
	b = private$b
	w = self$w

	y_star = exp(mu - sigma2)

	if (maximize) {
		if (y_star > b) {
			out = w(b, log = TRUE)
		} else if (y_star < a) {
			out = w(a, log = TRUE)
		} else {
			out = w(y_star, log = TRUE)
		}
	} else {
		out = min(w(a, log = TRUE), w(b, log = TRUE))
	}

	if (log) { return(out) } else { return(exp(out)) }
}

) # Close public
) # Close class
