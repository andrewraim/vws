library(vws)

x = 0

h = new(RUnivariateHelper,
	pdf = function(x, log) { dnorm(x = x, log = log) },
	cdf = function(x, lower, log) { pnorm(q = x, lower.tail = lower, log.p = log) },
	quantile = function(x, lower, log) { qnorm(p = x, lower.tail = lower, log.p = log) },
	supp = function(x) { TRUE }
)

h$pdf(0, TRUE)
h$cdf(0.95, TRUE, FALSE)
h$quantile(0.95, TRUE, FALSE)
h$supp(0)

my_test_function(h)
my_test_function(x)
