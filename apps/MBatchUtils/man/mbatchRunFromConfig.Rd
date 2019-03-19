\name{mbatchRunFromConfig}
\alias{mbatchRunFromConfig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchRunFromConfig}
\description{Run MBatch from a configuration file. This is for reproduciblity or for use with Batch Effects Interface}
\usage{
mbatchRunFromConfig(theConfigFile, theOutputDir, theNaStrings,
theShaidyMapGen, theShaidyMapGenJava)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theConfigFile}{Full path to the configuration file. Documentation will be in the MBatch release documentation.}
  \item{theOutputDir}{Full path to directory to create output, and the function attempts to create the directory.}
  \item{theNaStrings}{This is the argument to mbatchLoadFiles, and should generally be "NA".}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar. Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
}
\details{
Does not return any meaningful value, but creates an empty file MBATCH_SUCCESS.txt in the output directory, if successful.
}
\value{
Does not return any meaningful value, but creates an empty file MBATCH_SUCCESS.txt in the output directory, if successful.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
