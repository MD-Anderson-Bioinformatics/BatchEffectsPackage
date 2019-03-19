\name{mbatchTrimData}
\alias{mbatchTrimData}
\title{mbatchTrimData}
\description{
Based on IQR, shorten the number of columns (genes) to keep the data below the given size.
}
\usage{
mbatchTrimData(theMatrix, theMaxSize = 1.5e+07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theMatrix}{A matrix with row names being sample identifier strings, column names being gene strings, and the data being numeric.}
  \item{theMaxSize}{An integer giving the maximum number of cells to be in the matrix.}
}
\details{
This function can be used to reduce the size of a data set to ensure that it will fit within memory.
}
\value{
A matrix such that the row names (samples) are unchanged, while the columns (genes) are shortened such that row count times column count is less than theMaxSize.
Genes are selected for inclusion based on IQR, which higher ones selected first. This is useful for batch effects, but may not be desirable for other computations.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
