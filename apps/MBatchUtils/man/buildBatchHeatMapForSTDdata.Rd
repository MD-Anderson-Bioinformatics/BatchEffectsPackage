\name{buildBatchHeatMapForSTDdata}
\alias{buildBatchHeatMapForSTDdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{buildBatchHeatMapForSTDdata}
\description{Build an NGCHM for the given Standardized Data dataset.}
\usage{
buildBatchHeatMapForSTDdata(theFile, theBatchType, theOutputDir,
theShaidyMapGen,
theShaidyMapGenJava = "/usr/bin/java")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{Full path to the matrix_data.tsv file from Standardized Data.
  There should be a batches.tsv in the same directory.}
  \item{theBatchType}{Batch Type for which to create the NGCHM.
  (Should be a column in batches.tsv.)}
  \item{theOutputDir}{Output directory for NGCHM.}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar.
  Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
}
\details{
Build an NGCHM named <theBatchType>_matrix.ngchm
}
\value{Does not return a meaningful value.}
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
