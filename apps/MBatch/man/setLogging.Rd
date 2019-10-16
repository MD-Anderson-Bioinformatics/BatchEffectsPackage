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
\author{
Tod Casasent \email{tdcasasent@mdanderson.org}
Nianxiang Zhang \email{nzhang@mdanderson.org}
}
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
