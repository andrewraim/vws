library(vws)
library(tidyverse)

source("functions.R")
Rcpp::sourceCpp("sampler.cpp")

set.seed(1234)

# ----- Generate data -----
z = 1
mu = 5
sigma2 = 0.5
lambda2 = 100

y_true = rlnorm(1, mu, sqrt(sigma2))
z = rnorm(1, y_true, sqrt(lambda2))

# ----- Compute the (normalized) target density -----
d_target = get_target_density(mu, sigma2, z, lambda2)

# ----- Sampler in C++ -----
n = 10000
out = r_lognormal_normal(n = n, z = z, mu = mu, sigma2 = sigma2,
	lambda2 = lambda2, N = 10, max_rejects = n * 10, report_period = n / 10)
vws::printf("Empirical rejection rate: %g%%",
	100 * sum(out$rejects) / (sum(out$rejects) + n))


# ----- Compare distribution of the draws to the target density -----
gg = ggplot() +
	geom_function(fun = d_target, args = list(log = FALSE)) +
	geom_density(data = data.frame(x = out$draws), aes(x = x), col = "blue", lty = 2) +
	xlab("x") +
	ylab("Density") +
	theme_linedraw()
print(gg)
