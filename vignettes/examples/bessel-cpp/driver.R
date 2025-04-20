library(tidyverse)
source("../bessel-r/bessel.R")

Rcpp::sourceCpp("sample.cpp")
# out = sample(n = 10000, a = 5, nu = 2, N = 50, max_rejects = 100000)
out = sample(n = 0, a = 5, nu = 2, N = 8)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, a = 5, nu = 2)

data.frame(x = out$draws) %>%
	ggplot() +
	geom_bar(aes(x, y = ..count.. / sum(..count..))) +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y)) +
	xlab("x") +
	ylab("Proportion") +
	theme_light()

data.frame(lbdd = out$lbdd) %>%
	mutate(step = row_number()) %>%
	ggplot() +
	geom_line(aes(step, lbdd)) +
	theme_light()
