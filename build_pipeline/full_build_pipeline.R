source("build_pipeline/document_package.R")
devtools::check(build_args = c("--resave-data"))
source("build_pipeline/build.R")
source("build_pipeline/render_vignette.R")
