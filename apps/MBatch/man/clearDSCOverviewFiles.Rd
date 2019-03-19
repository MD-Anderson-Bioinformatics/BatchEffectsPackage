\name{clearDSCOverviewFiles}
\alias{clearDSCOverviewFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delete All DSC.RData Files
}
\description{
Deletes all *DSC.RData files from the directory structure passed in to it.

After running mbatchAssess, there are *DSC.RData files created -- these can either be use to create a DSC Overview file with buildDSCOverviewFile, or they can be deleted using clearDSCOverviewFiles. This is handled outside of these functions, so all data sets can be processed by multiple calls to mbatchAssess in a parallel or distributed environment, and then the DSC Overview created from these files at the end of the processing.
}
\usage{
clearDSCOverviewFiles(theStartDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theStartDir}{String indicating top of directory structure to search for *DSC.RData objects.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tod Casasent \email{tdcasasent@mdanderson.org}
Nianxiang Zhang \email{nzhang@mdanderson.org}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{buildDSCOverviewFile}}
}
\examples{
\dontrun{
example in doRunBEA_Structures
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
