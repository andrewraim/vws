source("../common/plots.R")
Rcpp::sourceCpp("vmf-v1.cpp")
Rcpp::sourceCpp("vmf-v2.cpp")

n = 20000
kappa = 5
d = 4
N = 50
tol = 0.25
max_rejects = 10000
report = 1000

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = r_vmf_pre_v1(n, kappa, d, N, tol, max_rejects, report)

plot_density(out$draws)
plot_bounds(out$lbdd)


# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer
out2 = r_vmf_pre_v2(n, kappa, d, N, tol, max_rejects, report)

plot_density(out2$draws)
plot_bounds(out2$lbdd)
