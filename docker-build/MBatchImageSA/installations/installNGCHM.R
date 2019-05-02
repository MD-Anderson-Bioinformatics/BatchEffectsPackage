message("starting")

# do not use dependencies=TRUE, it loads way, way too much
install.packages("devtools", repos = "http://cloud.r-project.org")
library(devtools)
devtools::install_github("bmbroom/tsvio")
devtools::install_github("bmbroom/NGCHMR")

message("done")
