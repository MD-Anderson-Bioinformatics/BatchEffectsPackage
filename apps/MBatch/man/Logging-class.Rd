\name{Logging-class}
\Rdversion{1.1}
\docType{class}
\alias{Logging-class}

\title{Class \code{"Logging"}}
\description{
An S4 object detailing logging settings.
}
\section{Objects from the Class}{
See also \code{\link{new}}
Objects can be created by calls of the form \code{
	setLogging(new("Logging",
	     theFile=system.file("logging", "mbatch.log", package="MBatch")))
}
}
\section{Slots}{
  \describe{
    \item{\code{mFile}:}{Object of class \code{"character"}
    The file to write log data to.
    To log only to the console, set the file equal to an empty string, "" }
    \item{\code{mLevelNamesToLog}:}{Object of class \code{"vector"}
    Vector of string of types of logging to perform. Should be limited to the
    following logging types:
    'DEBUG', 'TIMING', 'INFO', 'WARN', 'PERCENT', 'ERROR' }
    \item{\code{mLevelNames}:}{Object of class \code{"vector"}
    Vector of the possible logging types for which there are corresponding
    logging functions. Should be left as the
    following a vector of the following strings:
    'DEBUG', 'TIMING', 'INFO', 'WARN', 'PERCENT', 'ERROR' }
    \item{\code{mSeparator}:}{Object of class \code{"character"}
    Separator character(s) to use
    between logging items.  Defaults to single space. }
    \item{\code{mConsole}:}{Object of class \code{"logical"}
    Determines if logged values should always be sent to console.
    Set to true if you wish to have the log messages
    sent to the console when a log file has been specified.
    If the log file is set to "", then logging
    messages will default to the console. }
  }
}
\section{Methods}{
No methods defined with class "Logging" in the signature.
}
\references{
%%  ~~put references to the literature/web site here~~
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
See also \code{\link{new}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
\dontrun{
#		###Set new logging object logging all 6 types
#   ###of messages to the console (default values).
#		logger<-new("Logging")
#
#		###Set up new logging object to write only error
#   ###and warning messages to the file
#   ###test.log and not to the console.
#		logger<-new("Logging", theFile="test.log",
#          mLevelNamesToLog=c("WARN","ERROR"),
#          theConsole=FALSE)
}
}
\keyword{classes}
