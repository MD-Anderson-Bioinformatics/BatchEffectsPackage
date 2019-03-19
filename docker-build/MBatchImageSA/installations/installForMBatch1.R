message("starting")

install.packages("mclust", dependencies=TRUE, type="source", repos = "http://cloud.r-project.org")
install.packages("squash", dependencies=TRUE, type="source", repos = "http://cloud.r-project.org")
install.packages("rJava", dependencies=TRUE, type="source", repos = "http://cloud.r-project.org")

install.packages("httr", dependencies=TRUE, type="source", repos = "http://cloud.r-project.org")

# Used by NGCHM
install.packages("dunn.test", dependencies=TRUE, type="source", repos = "http://cloud.r-project.org")
install.packages("jsonlite", dependencies=TRUE, type="source", repos = "http://cloud.r-project.org")

message("done")
