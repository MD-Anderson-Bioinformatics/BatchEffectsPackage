message("starting")

source("http://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("RBGL")
biocLite("graph")
biocLite("Biobase")

install.packages(c("oompaBase", "ClassDiscovery", "PreProcess"), dependencies=TRUE, repos=c("http://cran.revolutionanalytics.com", "http://silicovore.com/OOMPA/"))

message("done")
