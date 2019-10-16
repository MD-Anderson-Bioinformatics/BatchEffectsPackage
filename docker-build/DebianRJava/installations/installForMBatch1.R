message("starting")

install.packages("mclust", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")
install.packages("squash", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")
install.packages("rJava", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")

install.packages("httr", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")

# Used by NGCHM
install.packages("dunn.test", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")
install.packages("jsonlite", dependencies=TRUE, type="source", repos = "http://cran.r-project.org")

message("done")
