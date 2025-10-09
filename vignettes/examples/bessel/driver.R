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
wseq = w(xseq, log = F)

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = r_bessel_v1(n, lambda, nu, N, tol, max_rejects*50, report, x = xseq)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))

plot_bounds(out$lbdd)

data.frame(x = xseq, h = out$hx) %>%
	ggplot() +
	geom_point(aes(x, h), pch = 1) +
	geom_point(aes(x, fseq), pch = 4) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

ggplot(out$df_weight) +
	geom_segment(aes(x = lo, xend = hi, y = exp(wmax)), col = "blue") +
	# geom_segment(aes(x = lo, xend = hi, y = exp(wmin)), col = "red") +
	coord_cartesian(xlim = c(NA, 15)) +
	# geom_function(fun = w, args = list(log = FALSE), xlim = c(0,15)) +
	geom_point(data = data.frame(x = xseq, w = wseq), aes(x, w)) +
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
	geom_point(aes(x, h), pch = 1) +
	geom_point(aes(x, fseq), pch = 4) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

ggplot(out$df_weight) +
	geom_segment(aes(x = lo, xend = hi, y = exp(wmax)), col = "blue") +
	# geom_segment(aes(x = lo, xend = hi, y = exp(wmin)), col = "red") +
	coord_cartesian(xlim = c(NA, 15)) +
	# geom_function(fun = w, args = list(log = FALSE), xlim = c(0,15)) +
	geom_point(data = data.frame(x = xseq, w = wseq), aes(x, w)) +
	ylab("Weight") +
	theme_minimal()

# ----- Version 3 -----
# Use custom optimization routine to compute constants in majorizer
out = r_bessel_v3(n, lambda, nu, N, lo = -0.1, hi = 100, tol, max_rejects,
	report, x = xseq)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))

plot_bounds(out$lbdd)

data.frame(x = xseq, h = out$hx) %>%
	ggplot() +
	geom_point(aes(x, h), pch = 1) +
	geom_point(aes(x, fseq), pch = 4) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

out$df_weight %>%
	mutate(w_lo = exp(beta0_max + beta1_max * lo)) %>%
	mutate(w_hi = exp(beta0_max + beta1_max * hi)) %>%
	ggplot() +
	geom_segment(aes(x = lo, xend = hi, y = w_lo, yend = w_hi), col = "blue") +
	geom_point(data = data.frame(x = xseq, w = wseq), aes(x, w)) +
	coord_cartesian(xlim = c(NA, 15)) +
	xlab("x") +
	ylab("Weight") +
	theme_minimal()




# ----- Debugging -----
xseq = 0:30

fseq1 = dpois(xseq, lambda = 10)
fseq2 = d_truncpois(xseq, lambda = 10, a = -1, b = 12)

cbind(xseq, fseq2)
sum(fseq2)

plot(xseq, fseq2)
points(xseq, fseq1, pch = 4)

incgamma(a = 0, x = 0, lower = F, log = T)

w_major = function(x, log = TRUE) {
	out = 0.693149 - 0.922785 * x
	if (log) { return(out) } else { return(exp(out)) }
}

w_minor = function(x, log = TRUE) {
	out = 0.693135 + 12.8155 * x
	if (log) { return(out) } else { return(exp(out)) }
}

curve(w(x), xlim = c(0, 100))
curve(w_major(x), add = T, col = "blue")
curve(w_minor(x), add = T, col = "red")



