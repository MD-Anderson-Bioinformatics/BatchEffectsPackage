\name{mbatchIncludeExcludeData}
\alias{mbatchIncludeExcludeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchIncludeExcludeData}
\description{Makes sure that the appropriate list(s) of included/excluded
samples/genes is enforced and returned the BEA Data object.}
\usage{
mbatchIncludeExcludeData(theBeaData, theIncludeSamples = NULL,
theIncludeGenes = NULL, theExcludeSamples = NULL, theExcludeGenes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBeaData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theIncludeSamples}{A list of sample ids.
  Any sample ids not in this list are removed from the data.}
  \item{theIncludeGenes}{A list of gene ids.
  Any gene ids not in this list are removed from the data.}
  \item{theExcludeSamples}{A list of sample ids.
  Any sample ids in this list are removed from the data.}
  \item{theExcludeGenes}{A list of gene ids.
  Any gene ids in this list are removed from the data.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an instance of BEA_DATA \code{\link{BEA_DATA-class}}
with the updated data and batch information.
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
