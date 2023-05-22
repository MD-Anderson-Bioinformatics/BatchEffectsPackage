\name{setGlobalMBatchEnv}
\alias{setGlobalMBatchEnv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{setGlobalMBatchEnv}
\description{
Set Python Conda environment name.
Default expected environment is gendev.
}
\usage{
setGlobalMBatchEnv(theValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theValue}{String giving conda environment name for PyMBatch package.}
}
\details{
Set Python Conda environment name.
Default expected environment is gendev.
}
\value{
Returns results of call to getGlobalMBatchEnv.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{getGlobalMBatchEnv}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theValue)
{
    assign("MBATCH_PYTHON_ENV", theValue, GLOBAL_MBATCH_ENV)
  }
}
