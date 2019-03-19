\name{mbatchFilterData}
\alias{mbatchFilterData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchFilterData}
\description{
Filter data based on IQR, SD, MAD and batches to remove and keep.
}
\usage{
mbatchFilterData(theBeaData,
theBatchTypeAndValuePairsToRemove=list(c("*", "unknown"), c("*", "Unknown")),
theBatchTypeAndValuePairsToKeep=list(list("Type", c("01", "03", "05"))),
theBatchTypesToRemove=NULL, theMinIqr = 0, theMinSd = 0, theMinMad = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBeaData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types) and
  the value to remove. list() indicates none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types) and
  a vector of the the value(s) to keep. list() indicates none while
  NULL will cause an error. }
  \item{theBatchTypesToRemove}{Object of class \code{"vector"}
  A vector of strings giving batch types to remove. NULL means none,
  c("Type") means remove Type.}
  \item{theMinIqr}{Object of class \code{"numeric"} Minimum acceptable
  Interquartile Range (MinIqr) value for gene values. (See \code{\link{IQR}}
  from the R stats package.) Zero (0) means use all samples. }
  \item{theMinSd}{Object of class \code{"numeric"} Minimum acceptable
  Standard Deviation (MinSd) value for gene values. (See \code{\link{sd}}
  from the R stats package.) Zero (0) means use all samples. }
  \item{theMinMad}{Object of class \code{"numeric"} Minimum acceptable
  Median Absolute Deviation (MinMad) value for gene values.
  (See \code{\link{mad}} from the R stats package.).
  Zero (0) means use all samples.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
