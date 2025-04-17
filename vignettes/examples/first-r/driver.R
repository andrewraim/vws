source("sample.R")

out = sample(n = 10000, kappa = 5, d = 4, N = 50)
hist(out$draws)
