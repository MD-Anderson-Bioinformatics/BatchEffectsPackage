\name{readAsMatrix}
\alias{readAsMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{readAsMatrix}
\description{
This function takes a cross-platform tab-delimited matrix file with
column and optional row labels and reads it in as a matrix.
The data is handled as double (numeric) values.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
readAsMatrix(theFile, thePar = "-Xmx2000m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{The cross-platform tab-delimited matrix file with columns
  and optional row labels to be read. If row labels are included, there
  is an initial tab on the column label row.}
  \item{thePar}{Parameters to pass to the Java VM --
  by default this requests 8GB of memory.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The matrix read in.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function (theFile, thePar = "-Xmx2000m")
#{
#    myClass1 <- system.file("ReadRJava", "ReadRJava.jar", package = "MBatch")
#    myJavaJars <- file.path(myClass1, fsep = .Platform$path.sep)
#    logDebug("readAsMatrix - thePar ", thePar)
#    logDebug("readAsMatrix - theFile ", theFile)
#    logDebug("readAsMatrix - Calling .jinit ", myJavaJars)
#    .jinit(classpath = myJavaJars, force.init = TRUE, parameters = thePar)
#    logDebug("readAsMatrix - .jinit complete")
#    logDebug("readAsMatrix before java")
#    objJavaFile <- .jcall("org/mda/readrjava/ReadRJava",
#        returnSig = "Lorg/mda/readrjava/JavaFile;",
#        method = "loadDoubleData", .jnew("java/lang/String",
#            theFile), TRUE, TRUE, TRUE)
#    logDebug("readAsMatrix after java")
#    myData <- .jcall(objJavaFile, "[D", "getmDoubleData")
#    myCols <- .jcall(objJavaFile, "[S", "getmColumns")
#    myRows <- .jcall(objJavaFile, "[S", "getmRows")
#    myCols <- myCols[2:length(myCols)]
#    logDebug("readAsMatrix - length(myData) ", length(myData))
#    logDebug("readAsMatrix - length(myCols) ", length(myCols))
#    logDebug("readAsMatrix - length(myRows) ", length(myRows))
#    return(matrixWithIssues(myData, nrow = length(myRows),
#        ncol = length(myCols),
#        byrow = TRUE, dimnames = list(myRows, myCols)))
#  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
