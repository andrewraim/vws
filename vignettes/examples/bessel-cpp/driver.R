library(tidyverse)
source("../bessel-r/bessel.R")

Rcpp::sourceCpp("sample.cpp")
out = sample(n = 20000, a = 5, nu = 2, N = 20, max_rejects = 100000)
#out = sample(n = 1, a = 5, nu = 2, N = 25)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, a = 5, nu = 2)

data.frame(x = out$draws) %>%
	ggplot() +
	geom_bar(aes(x, y = after_stat(count / sum(count))), fill = "white", col = "black") +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y)) +
	xlab("x") +
	ylab("Empirical Density") +
	theme_light()

data.frame(lbdd = out$lbdd) %>%
	mutate(step = row_number()) %>%
	ggplot() +
	geom_line(aes(step, lbdd)) +
	theme_light()
