library(vws)

helper = RUnivariateHelper$new()

helper = RUnivariateHelper$new(
	d = function(x, log) { dnorm(x = x, log = log) },
	p = function(x, lower, log) { pnorm(x = x, lower = lower, log = log) },
	q = function(x, lower, log) { qnorm(x = x, lower = lower, log = log) },
	s = function(x) { TRUE }
)
helper$d(10, TRUE)
