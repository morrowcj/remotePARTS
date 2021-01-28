# build.r
## needs to be executed from remotePARTS/

Rcpp::compileAttributes()
devtools::document()
devtools::build()
if(file.exists("doc/Alaska.html")){
  # rmarkdown::render("vignettes/Alaska.Rmd", output_dir = "docs/")
  file.copy(from = "doc/Alaska.html", to = "docs/Alaska.html",
            overwrite = TRUE)
}
