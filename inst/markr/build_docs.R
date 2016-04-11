#library(staticdocs)
library(markr)
library(alakazam)

# Directories
pkg_path <- "/home/jason/workspace/igpipeline/alakazam"
doc_path <- "/home/jason/workspace/igpipeline/prototype-docs/docs"

# As a single step
build_mkdocs(pkg_path, doc_path=doc_path, yaml=F)

# As individual steps
#dir.create(doc_path, recursive=T)
#pkg <- staticdocs::as.sd_package(pkg=pkg_path, site_path=doc_path)
#pkg$topics <- build_md_topics(pkg)
#pkg$vignettes <- build_md_vignettes(pkg)
#pkg$index <- build_md_index(pkg)
#build_yaml(pkg)
