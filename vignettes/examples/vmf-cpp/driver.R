source("../common/plots.R")
Rcpp::sourceCpp("sample.cpp")
Rcpp::sourceCpp("sample2.cpp")

n = 20000
kappa = 5
d = 4
N = 50
tol = 0.25
max_rejects = 10000
report = 1000

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = sample(n, kappa, d, N, tol, max_rejects, report)

plot_density(out$draws)
plot_bounds(out$lbdd)


# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer
out2 = sample2(n, kappa, d, N, tol, max_rejects, report)

plot_density(out2$draws)
plot_bounds(out2$lbdd)
