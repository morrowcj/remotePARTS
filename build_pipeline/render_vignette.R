
rmarkdown::render(
  input = "vignettes/Alaska.Rmd",
  output_format = NULL,
  output_file = "Alaska.html",
  output_dir = "docs/"
)
# TODO use devtools::build_rmd() instead?
file.copy(from = "docs/Alaska.html", to = "vignettes/", overwrite = TRUE)
