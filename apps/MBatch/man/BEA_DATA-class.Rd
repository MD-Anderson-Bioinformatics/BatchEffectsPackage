\name{BEA_DATA-class}
\Rdversion{1.1}
\docType{class}
\alias{BEA_DATA-class}
\alias{initialize,BEA_DATA-method}

\title{Class \code{"BEA_DATA"}}
\description{
An S4 class, BEA_DATA provides slots for the MBatch data matrix, batch data.frame and covariate data.frame.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BEA_DATA", theData, theBatches, theCovariates)}.
See the Slots section for a description of the objects.
If you have no covariate data, use an empty data.frame created with \code{data.frame()}
}
\section{Slots}{
  \describe{
    \item{\code{mData}:}{Object of class \code{"matrix"} A matrix where the colnames are sample ids and the rownames are gene equivalents. All names should be strings, not factors. }
    \item{\code{mBatches}:}{Object of class \code{"data.frame"} A data.frame where the column "names" are batch types. The first batch "type" is "Sample". All names and values should be strings, not factors or numeric. }
    \item{\code{mCovariates}:}{Object of class \code{"data.frame"} A data.frame where the column "names" are covariate types. The first covariate "type" is "Sample". All names and values should be strings, not factors or numeric. }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "BEA_DATA")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
Tod Casasent \email{tdcasasent@mdanderson.org}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
#showClass("BEA_DATA")
}
\keyword{classes}
