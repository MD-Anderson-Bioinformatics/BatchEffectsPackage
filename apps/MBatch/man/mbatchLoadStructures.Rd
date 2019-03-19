\name{mbatchLoadStructures}
\alias{mbatchLoadStructures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchLoadStructures}
\description{
Take the given structures and build a BEA_DATA object to contain them, after doing sorting and making sure data is consistent.
}
\usage{
mbatchLoadStructures(theGeneMatrix, theBatchDataframe, theCovariatedataframe=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theGeneMatrix}{Object of class \code{"matrix"} A matrix where the colnames are sample ids and the rownames are gene equivalents. All names should be strings, not factors. }
  \item{theBatchDataframe}{Object of class \code{"data.frame"} A data.frame where the column "names" are batch types. The first batch "type" is "Sample". All names and values should be strings, not factors or numeric. }
  \item{theCovariatedataframe}{Object of class \code{"data.frame"} NULL or a data.frame where the column "names" are covariate types. The first covariate "type" is "Sample". All names and values should be strings, not factors or numeric. }
}
\details{
}
\value{
An instance of BEA_DATA \code{\link{BEA_DATA-class}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tod Casasent \email{tdcasasent@mdanderson.org}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
