\name{writeAsMatrix}
\alias{writeAsMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writeAsMatrix}
\description{
This function takes a matrix and writes it out as a cross-platform tab-delimited matrix file with column and row labels.
The matrix data should be doubles.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
writeAsMatrix(theFile, theMatrix, thePar = "-Xmx2000m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{Full path to the output file.}
  \item{theMatrix}{The matrix to write out--requires column names and double (numeric) data.}
  \item{thePar}{Parameters to pass to the Java VM -- by default this requests 8GB of memory.}
}
\value{
Returns TRUE if the save succeeded.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function (theFile, theMatrix, thePar = "-Xmx2000m")
#{
#    myClass1 <- system.file("ReadRJava", "ReadRJava.jar", package = "MBatch")
#    myJavaJars <- file.path(myClass1, fsep = .Platform$path.sep)
#    logDebug("writeAsMatrix - thePar ", thePar)
#    logDebug("writeAsMatrix - theFile ", theFile)
#    myCols <- as.vector(colnames(theMatrix))
#    myRows <- rownames(theMatrix)
#    myData <- as.numeric(as.vector(t(theMatrix)))
#    if (!is.null(myRows)) {
#        myRows <- as.vector(myRows)
#    }
#    logDebug("writeAsMatrix - length(myData) ", length(myData))
#    logDebug("writeAsMatrix - length(myCols) ", length(myCols))
#    logDebug("writeAsMatrix - length(myRows) ", length(myRows))
#    logDebug("writeAsMatrix - Calling .jinit ", myJavaJars)
#    .jinit(classpath = myJavaJars, force.init = TRUE, parameters = thePar)
#    logDebug("writeAsMatrix - .jinit complete")
#    logDebug("writeAsMatrix before java")
#    if (is.null(myRows)) {
#        success <- .jcall("org/mda/readrjava/ReadRJava", returnSig = "Z",
#            method = "writeDoubleData_Column", .jnew("java/lang/String",
#                theFile), .jarray(myCols), .jarray(myData))
#    }
#    else {
#        success <- .jcall("org/mda/readrjava/ReadRJava", returnSig = "Z",
#            method = "writeDoubleData_All", .jnew("java/lang/String",
#                theFile), .jarray(myCols), .jarray(myRows), .jarray(myData))
#    }
#    logDebug("writeAsMatrix after java")
#    logDebug("writeAsMatrix success=", success)
#    return(success)
#  }
}
