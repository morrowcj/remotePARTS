source("build_pipeline/document_package.R")
devtools::check()
source("build.R")
source("render_vignette.R")
