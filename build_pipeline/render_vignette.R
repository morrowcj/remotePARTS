rmarkdown::render(input = "vignettes/Alaska.Rmd",
                  output_format = "html_document",
                  output_file = "Alaska.html",
                  output_dir = "docs/")
file.copy(from = "docs/Alaska.html", to = "vignettes/")
