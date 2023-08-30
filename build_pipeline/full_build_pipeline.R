source("build_pipeline/document_package.R")
devtools::check()
source("build_pipeline/build.R")
source("build_pipeline/render_vignette.R")
