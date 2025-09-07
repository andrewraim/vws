library(vws)

source(system.file("examples/lnorm-norm-const-v1/functions.R", package = "vws"))
source(system.file("examples/lnorm-norm-const-v2/CustomConstRegion.R", package = "vws"))

# ----- Generate data -----
mu = 5
sigma2 = 0.5
lambda2 = 100

set.seed(1234)
y_true = rlnorm(1, mu, sqrt(sigma2))
z = rnorm(1, y_true, sqrt(lambda2))
print(y_true)
print(z)

# ----- Construct proposal -----
w = function(y, log = TRUE) {
        out = -log(y) - (log(y) - mu)^2 / (2*sigma2) + log(y > 0)
        out[y == 0] = -Inf
        if (log) { return(out) } else { return(exp(out)) }
}

helper = normal_helper(mean = z, sd = sqrt(lambda2))

support = CustomConstRegion$new(a = 0, b = Inf, w = w, g = helper)
regions = list(support)

# ----- Refine proposal -----
h_init = FMMProposal$new(regions)
refine_out = refine(h_init, N = 30)
h = refine_out$h

# ----- Rejection sampling -----
ctrl = rejection_control(report = 5000, extra_outputs = TRUE)
out = rejection(h, n = 10000, control = ctrl)
y = unlist(out$draws)

cat("Percent of proposed draws which were rejected:",
	sum(out$rejects) / (length(y) + sum(out$rejects)) * 100)

# ----- Plot draws vs. target density -----
gg = data.frame(y = y) %>%
	ggplot() +
	geom_histogram(aes(x = y, y = after_stat(density)), col = "black",
		fill = NA, bins = 25) +
	geom_function(fun = d_target, args = list(mu = mu, sigma2 = sigma2, z = z,
		lambda2 = lambda2, log = FALSE), lty = 2) +
	xlab("y") +
	ylab("Density") +
	theme_minimal()
print(gg)

# ----- Plot log f(x) vs. log h(x) -----
log_h = function(y) { sapply(y, h$d, log = TRUE) }
log_f = function(y) { d_target(y, mu, sigma2, z, lambda2, log = TRUE) }

xlim = range(y)
ggplot() +
	geom_function(fun = log_f, col = "orange", lty = 1) +
	geom_function(fun = log_h, col = "black", lty = 2) +
	scale_x_continuous(limits = xlim) +
	xlab("y") +
	ylab("Density") +
	theme_minimal()
