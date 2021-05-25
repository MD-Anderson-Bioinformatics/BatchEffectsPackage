\name{readAsDataFrame}
\alias{readAsDataFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{readAsDataFrame}
\description{
This function takes a cross-platform tab-delimited matrix file with column
and optional row labels and reads it in as a data.frame.
The data.frame data is handled as strings.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
readAsDataFrame(theFile, thePar = "-Xmx2000m", theUnknownString="Unknown")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{The cross-platform tab-delimited data.frame file with
  columns and optional row labels to be read. If row labels are included,
  there is an initial tab on the column label row.}
  \item{thePar}{Parameters to pass to the Java VM --
  by default this requests 8GB of memory.}
  \item{theUnknownString}{String to replace empty strings with. Defaults to 'Unknown'.}
}
\value{
Returns the data.frame read from file.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function (theFile, thePar = "-Xmx2000m")
#{
#    myJavaJars <- getJarsFromDir(dirname(system.file("ReadRJava", "ReadRJava.jar", package="MBatch")))
#    logDebug("readAsDataFrame - thePar ", thePar)
#    logDebug("readAsDataFrame - theFile ", theFile)
#    logDebug("readAsDataFrame - Calling .jinit ", myJavaJars)
#    .jinit(classpath = myJavaJars, force.init = TRUE, parameters = updateJavaParameters(thePar))
#    logDebug("readAsDataFrame - .jinit complete")
#    logDebug("readAsDataFrame before java")
#    objJavaFile <- .jcall("edu/mda/bcb/readrjava/ReadRJava",
#        returnSig = "Ledu/mda/bcb/readrjava/JavaFile;",
#        method = "loadStringData", .jnew("java/lang/String",
#            theFile), TRUE, FALSE, TRUE)
#    logDebug("readAsDataFrame after java")
#    myData <- .jcall(objJavaFile, "[S", "getmStringData")
#    myCols <- .jcall(objJavaFile, "[S", "getmColumns")
#    logDebug("readAsDataFrame - length(myData) ", length(myData))
#    logDebug("readAsDataFrame - length(myCols) ", length(myCols))
#    return(data.frame(matrixWithIssues(myData, ncol = length(myCols),
#        byrow = TRUE, dimnames = list(NULL, myCols))))
#  }
}
