source("../common/plots.R")
Rcpp::sourceCpp("vmf-v1.cpp")
Rcpp::sourceCpp("vmf-v1-xptr.cpp")
Rcpp::sourceCpp("vmf-v1-module.cpp")
Rcpp::sourceCpp("vmf-v2.cpp")
Rcpp::sourceCpp("vmf-v3.cpp")

n = 20000
kappa = 5
d = 4
N = 50
tol = 0.10
max_rejects = 10000
report = 1000

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = r_vmf_pre_v1(n, kappa, d, N, tol, max_rejects, report)

plot_density(out$draws) +
	geom_function(fun = d_target, args = list(kappa = kappa, d = d), lty = 2)
plot_bounds(out$lbdd)

# ----- Version 1 with Xptr -----
tryCatch({
	g = test_vector()
	refine(g, N, tol)
}, error = function(e) {
	cat("Caught an expected error:\n")
	print(e)
})

h = vmf_pre_v1_xptr(kappa, d)
lbdd = refine(h, N, tol)
out = draw(h, n, max_rejects, report)
rm(h)

plot_density(out$draws) +
	geom_function(fun = d_target, args = list(kappa = kappa, d = d), lty = 2)
plot_bounds(lbdd)

# ----- Version 1 with Rcpp Module -----
h = new(RcppProposal, kappa, d)
lbdd = h$refine(N-1, tol)
out = h$draw(n, max_rejects, report)

plot_density(out$draws) +
	geom_function(fun = d_target, args = list(kappa = kappa, d = d), lty = 2)
plot_bounds(lbdd)

# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer
out = r_vmf_pre_v2(n, kappa, d, N, tol, max_rejects, report)

plot_density(out$draws) +
	geom_function(fun = d_target, args = list(kappa = kappa, d = d), lty = 2)
plot_bounds(out$lbdd)

# ----- Version 3 -----
# Use a linear majorizer
out = r_vmf_pre_v3(n, kappa, d, N, tol, max_rejects, report)

plot_density(out$draws) +
	geom_function(fun = d_target, args = list(kappa = kappa, d = d), lty = 2)
plot_bounds(out$lbdd)
