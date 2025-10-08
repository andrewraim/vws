source("../common/plots.R")
source("bessel.R")
Rcpp::sourceCpp("bessel-v1.cpp")
Rcpp::sourceCpp("bessel-v2.cpp")
Rcpp::sourceCpp("bessel-v3.cpp")

n = 20000
lambda = 5
nu = 2
N = 20
tol = 0.01
max_rejects = 50000
report = 10000

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = r_bessel_v1(n, lambda, nu, N, tol, max_rejects, report)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, lambda, nu)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))
plot_bounds(out$lbdd)

# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer.
out = r_bessel_v2(n, lambda, nu, N, tol, max_rejects, report)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, lambda, nu)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))
plot_bounds(out$lbdd)

# ----- Version 3 -----
# Use custom optimization routine to compute constants in majorizer
out = r_bessel_v3(n, lambda, nu, N, lo = -0.5, hi = 1e8, tol, max_rejects, report)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, lambda, nu)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))
plot_bounds(out$lbdd)

# Debugging
xseq = 0:30

fseq1 = dpois(xseq, lambda = 10)
fseq2 = d_truncpois(xseq, lambda = 10, a = -1, b = 12)

cbind(xseq, fseq2)
sum(fseq2)

plot(xseq, fseq2)
points(xseq, fseq1, pch = 4)

incgamma(a = 0, x = 0, lower = F, log = T)

w = function(x, log = TRUE) {
	out = -lgamma(x + nu + 1)
	if (log) { return(out) } else { return(exp(out)) }
}

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



