source("../common/plots.R")
source("../common/bessel.R")

n = 20000
a = 5
nu = 2
N = 20
max_rejects = 50000
report = 10000

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
Rcpp::sourceCpp("sample.cpp")
out = sample(n, a, nu, N, max_rejects, report)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, a, nu)

plot_pmf(out$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))
plot_bounds(out$lbdd)

# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer
Rcpp::sourceCpp("sample2.cpp")
out2 = sample2(n, a, nu, N, max_rejects, report)

xseq = seq(0, max(out2$draws))
fseq = d_bessel(xseq, a, nu)

plot_pmf(out2$draws) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y))
plot_bounds(out2$lbdd)
