\name{getGlobalMBatchEnv}
\alias{getGlobalMBatchEnv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getGlobalMBatchEnv
}
\description{
Return Python environment name used by MBatch.
}
\usage{
getGlobalMBatchEnv()
}
%- maybe also 'usage' for other objects documented here.
\details{
Return Python environment name used by MBatch.
}
\value{
String Python environment name used by MBatch.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{setGlobalMBatchEnv}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function ()
{
    return(get("MBATCH_PYTHON_ENV", envir = GLOBAL_MBATCH_ENV))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
