\name{setLogging}
\alias{setLogging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
setLogging
}
\description{
Used to set up a shared logging object to be used by logging functions and
other routines in the package. See the Logging clas for details.
}
\usage{
setLogging(theLogger)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theLogger}{The logger object to share.}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
See Also as \code{\link{Logging-class}}
}
\examples{
\dontrun{
#	## see examples in doRunBEA_Files and doRunBEA_Structures
#	setLogging(new("Logging", theFile=theLogFile,
#theLevelNamesToLog=c('TIMING', 'ERROR', 'WARN', 'INFO')))
#
#	###setLogging(new("Logging", theFile=theLogFile))
#
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
