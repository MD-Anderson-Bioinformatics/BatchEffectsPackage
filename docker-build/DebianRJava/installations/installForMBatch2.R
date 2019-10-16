message("starting")

install.packages("BiocManager", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")
BiocManager::install(c("limma", "RBGL", "graph", "Biobase"))

install.packages(c("oompaBase", "ClassDiscovery", "PreProcess"), dependencies=TRUE, repos=c("http://cran.r-project.org", "http://silicovore.com/OOMPA/"))

message("done")
