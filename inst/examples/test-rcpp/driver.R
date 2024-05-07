library(vws)
library(tidyverse)
Rcpp::sourceCpp("sampler.cpp", rebuild = TRUE)

# ----- Generate data -----
z = 1
mu = 5
sigma2 = 0.5
lambda2 = 100

N = 30
n = 100000

set.seed(1234)
y_true = rlnorm(1, mu, sqrt(sigma2))
z = rnorm(1, y_true, sqrt(lambda2))

# ----- Sampler in R -----

w = function(y, log = TRUE) {
	out = -log(y) - (log(y) - mu)^2 / (2*sigma2) + log(y > 0)
	out[y == 0] = -Inf
	if (log) { return(out) } else { return(exp(out)) }
}

helper = normal_univariate_helper(mean = z, sd = sqrt(lambda2))

support = UnivariateConstRegion$new(a = 0, b = Inf, w = w, g = helper)
regions = list(support)
h_init = FMMProposal$new(regions)

adapt_out = adapt(h_init, N = 9, report = 10)
h = adapt_out$h
print(h, n = N - 1)

# bdd = tail(exp(adapt_out$log_bdd_hist), 1)

ctrl = rejection_control(report = 5000, extra_outputs = TRUE, max_rejects = 10000)
out1 = rejection(h, n = n, control = ctrl)
y = unlist(out1$draws)
hist(y)

sum(out1$rejects) / (sum(out1$rejects) + n)


# ----- Sampler in C++ -----
out2 = r_lognormal_normal(n = n, z = z, mu = mu, sigma2 = sigma2,
	lambda2 = lambda2, max_rejects = 10000, report_period = 10000)
hist(out2$draws)

sum(out2$rejects) / (sum(out2$rejects) + n)

ggplot() +
	geom_density(data = data.frame(x = y), aes(x = x)) +
	geom_density(data = data.frame(x = out2$draws), aes(x = x), col = "blue", lty = 2) +
	xlab("x") +
	ylab("Density") +
	theme_linedraw()
