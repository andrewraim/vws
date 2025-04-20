library(tidyverse)
Rcpp::sourceCpp("sample.cpp")

out = sample(n = 10000, kappa = 5, d = 4, N = 50, tol = 0.25)

data.frame(x = out$draws) %>%
	ggplot() +
	geom_density(aes(x)) +
	xlab("x") +
	ylab("Empirical Density") +
	theme_light()

data.frame(lbdd = out$lbdd) %>%
	mutate(step = row_number()) %>%
	ggplot() +
	geom_line(aes(step, lbdd)) +
	theme_light()

