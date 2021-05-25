\name{writeAsDataframe}
\alias{writeAsDataframe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writeAsDataframe}
\description{
This function takes a data.frame and writes it out as a cross-platform
tab-delimited matrix file with column and row labels.
The data.frame data is handled as strings.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
writeAsDataframe(theFile, theDataframe, thePar = "-Xmx2000m",
theIncludeRowNamesFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{Full path to the output file.}
  \item{theDataframe}{The data.frame to write out--requires column names.}
  \item{thePar}{Parameters to pass to the Java VM --
  by default this requests 8GB of memory.}
  \item{theIncludeRowNamesFlag}{By default this is false,
  which skips the row names.}
}
\value{
Returns TRUE if the save succeeded.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function (theFile, theDataframe, thePar = "-Xmx2000m",
#theIncludeRowNamesFlag = FALSE)
#{
#    myJavaJars <- getJarsFromDir(dirname(system.file("ReadRJava", "ReadRJava.jar", package="MBatch")))
#    logDebug("writeAsDataframe - thePar ", thePar)
#    logDebug("writeAsDataframe - theFile ", theFile)
#    myCols <- as.vector(colnames(theDataframe))
#    myRows <- NULL
#    myData <- as.vector(t(theDataframe))
#    if (TRUE == theIncludeRowNamesFlag) {
#        myRows <- as.vector(rownames(theDataframe))
#    }
#    logDebug("writeAsDataframe - length(myData) ", length(myData))
#    logDebug("writeAsDataframe - length(myCols) ", length(myCols))
#    logDebug("writeAsDataframe - length(myRows) ", length(myRows))
#    logDebug("writeAsDataframe - Calling .jinit ", myJavaJars)
#    .jinit(classpath = myJavaJars, force.init = TRUE, parameters = updateJavaParameters(thePar))
#    logDebug("writeAsDataframe - .jinit complete")
#    logDebug("writeAsDataframe before java")
#    if (is.null(myRows)) {
#        success <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig = "Z",
#            method = "writeStringData_Column", .jnew("java/lang/String",
#                theFile), .jcastToArray(myCols), .jcastToArray(myData))
#    }
#    else {
#        success <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig = "Z",
#            method = "writeStringData_All", .jnew("java/lang/String",
#                theFile), .jcastToArray(myCols), .jcastToArray(myRows),
#            .jcastToArray(myData))
#    }
#    logDebug("writeAsDataframe after java")
#    logDebug("writeAsDataframe success=", success)
#    return(success)
#  }
}
