source("../common/plots.R")
source("bessel.R")
Rcpp::sourceCpp("bessel-v1.cpp")
Rcpp::sourceCpp("bessel-v2.cpp")
Rcpp::sourceCpp("bessel-v3.cpp")

n = 20000
lambda = 10
nu = 2
N = 10
tol = 0.01
max_rejects = 50000
report = 10000

xseq = seq(0, 15)
fseq = d_bessel(xseq, lambda, nu)
lwseq = w(xseq, log = T)

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = r_bessel_v1(n, lambda, nu, N, tol, max_rejects*50, report, x = xseq)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))

plot_bounds(out$lbdd)

data.frame(x = xseq, h = out$hx) %>%
	ggplot() +
	geom_bar(aes(x, h), stat = "identity", fill = NA, col = "black") +
	geom_point(aes(x, fseq), pch = 2) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

ggplot(out$df_weight) +
	geom_segment(aes(x = lo, xend = hi, y = lwmax), col = "blue") +
	geom_point(aes(x = hi, y = lwmax), col = "blue") +
	geom_point(data = data.frame(x = xseq, w = lwseq), aes(x, w), pch = 2) +
	coord_cartesian(xlim = c(NA, 15), ylim = c(min(lwseq), NA)) +
	xlab("x") +
	ylab("Weight") +
	theme_minimal()

# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer.
out = r_bessel_v2(n, lambda, nu, N, tol, max_rejects, report, x = xseq)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))

plot_bounds(out$lbdd)

data.frame(x = xseq, h = out$hx) %>%
	ggplot() +
	geom_bar(aes(x, h), stat = "identity", fill = NA, col = "black") +
	geom_point(aes(x, fseq), pch = 2) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

ggplot(out$df_weight) +
	geom_segment(aes(x = lo, xend = hi, y = lwmax), col = "blue") +
	geom_point(aes(x = hi, y = lwmax), col = "blue") +
	geom_point(data = data.frame(x = xseq, w = lwseq), aes(x, w), pch = 2) +
	coord_cartesian(xlim = c(NA, 15), ylim = c(min(lwseq), NA)) +
	xlab("x") +
	ylab("Weight") +
	theme_minimal()

# ----- Version 3 -----
# Use custom optimization routine to compute constants in majorizer
out = r_bessel_v3(n, lambda, nu, lo = -0.1, hi = 1e5, N, tol, max_rejects,
	report, x = xseq)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))

plot_bounds(out$lbdd)

data.frame(x = xseq, h = out$hx) %>%
	ggplot() +
	geom_bar(aes(x, h), stat = "identity", fill = NA, col = "black") +
	geom_point(aes(x, fseq), pch = 2) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

out$df_weight %>%
	mutate(w_lo = beta0_max + beta1_max * lo) %>%
	mutate(w_hi = beta0_max + beta1_max * hi) %>%
	ggplot() +
	geom_segment(aes(x = lo, xend = hi, y = w_lo, yend = w_hi), col = "blue") +
	geom_point(aes(x = hi, y = w_hi), col = "blue") +
	geom_point(data = data.frame(x = xseq, y = lwseq), aes(x, y), pch = 2) +
	coord_cartesian(xlim = c(NA, 15), ylim = c(min(lwseq), NA)) +
	xlab("x") +
	ylab("Weight") +
	theme_minimal()
