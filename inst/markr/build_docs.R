library(markr)
library(alakazam)

# Directories
pkg_path <- "/home/jason/workspace/igpipeline/alakazam"
doc_path <- "/home/jason/workspace/igpipeline/docs"

# Build
build_mkdocs(pkg_path, doc_path=doc_path, yaml=F)