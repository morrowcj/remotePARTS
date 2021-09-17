# build.r
## needs to be executed from remotePARTS/

Rcpp::compileAttributes()
devtools::document()
devtools::build(vignettes = FALSE)
devtools::build_vignettes()
if(file.exists("doc/Alaska.html")){
  # rmarkdown::render("vignettes/Alaska.Rmd", output_dir = "docs/")
  file.copy(from = "doc/Alaska.html", to = "docs/Alaska.html",
            overwrite = TRUE)
}
if(file.exists("doc/GLS_optimization.html")){
  file.copy(from = "doc/GLS_optimization.html",
            to = "docs/GLS_optimization.html", overwrite = TRUE)
}
