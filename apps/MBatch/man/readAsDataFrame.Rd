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
readAsDataFrame(theFile, thePar = "-Xmx2000m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{The cross-platform tab-delimited data.frame file with
  columns and optional row labels to be read. If row labels are included,
  there is an initial tab on the column label row.}
  \item{thePar}{Parameters to pass to the Java VM --
  by default this requests 8GB of memory.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the data.frame read from file.
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
#    logDebug("readAsDataFrame - thePar ", thePar)
#    logDebug("readAsDataFrame - theFile ", theFile)
#    logDebug("readAsDataFrame - Calling .jinit ", myJavaJars)
#    .jinit(classpath = myJavaJars, force.init = TRUE, parameters = thePar)
#    logDebug("readAsDataFrame - .jinit complete")
#    logDebug("readAsDataFrame before java")
#    objJavaFile <- .jcall("org/mda/readrjava/ReadRJava",
#        returnSig = "Lorg/mda/readrjava/JavaFile;",
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
