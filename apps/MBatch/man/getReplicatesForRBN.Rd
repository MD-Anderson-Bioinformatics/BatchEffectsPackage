\name{getReplicatesForRBN}
\alias{getReplicatesForRBN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getReplicatesForRBN}
\description{Function to find exact replicates in sample (column) names.}
\usage{
getReplicatesForRBN(matrix1, matrix2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix1}{First matrix, with sample ids in colnames() and
  feature names (like genes) in rownames().}
  \item{matrix2}{Second matrix, with sample ids in colnames() and
  feature names (like genes) in rownames().}
}
\details{
Find colnames that match between matrixes.
}
\value{
A vector of matching feature ids.
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
