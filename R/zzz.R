# To enable Quarto vignettes, from David C. Norris' post at:
# <https://github.com/quarto-dev/quarto-cli/discussions/2307#discussioncomment-3571271>
.onLoad <- function(libname, pkgname)
{
	tools::vignetteEngine(name = "quarto",
		package = pkgname,
		pattern = "[.]qmd$",
		weave = function(file, ..., encoding = "UTF-8") {
			packageStartupMessage("Custom Quarto vignette engine")
			## TODO: What is Quarto's expectation about encoding?
			## NB: output_format = "all" below might make a better default
			quarto::quarto_render(file, ..., output_format = "pdf")
		},
			tangle = tools::vignetteEngine("knitr::rmarkdown")$tangle,
			aspell = tools::vignetteEngine("knitr::rmarkdown")$aspell
		)
}
