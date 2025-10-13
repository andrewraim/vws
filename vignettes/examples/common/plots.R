library(tidyverse)

plot_density = function(x) {
	data.frame(x = x) %>%
		ggplot() +
		geom_density(aes(x)) +
		xlab("x") +
		ylab("Empirical Density") +
		theme_light()
}

plot_pmf = function(x) {
	data.frame(x = x) %>%
		ggplot() +
		geom_bar(aes(x, y = after_stat(count / sum(count))), fill = "white", col = "black") +
		xlab("x") +
		ylab("Empirical Probability") +
		theme_light()
}

plot_bounds = function(x) {
	data.frame(x = x) %>%
		mutate(step = row_number()) %>%
		ggplot() +
		geom_line(aes(step, x)) +
		xlab("Step") +
		ylab("Log of Rejection Bound") +
		theme_light()
}
