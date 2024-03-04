The vignette `vws.qmd` is a Quarto document. It can be compiled in Rstudio by
clicking the `Render` button of using the console as follows.

```r
quarto::quarto_render("vws.qmd")
```

It can be time-consuming to run the embedded code and completely render the
document. Evaluation of the embedded code can be enabled and disabled by
setting the following YAML header entry to `true` or `false`, respectively.

```yaml
execute:
  eval: false
```
