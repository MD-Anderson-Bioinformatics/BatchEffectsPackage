message("starting")

# do not use dependencies=TRUE, it loads way, way too much
install.packages("devtools", repos = "http://cran.r-project.org")
library(devtools)
install_github("js229/Vennerable")


message("done")
