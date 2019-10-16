\name{mutationBatchExtract}
\alias{mutationBatchExtract}
\title{Extract Mutation Data from DCC or GDC Files}
\description{This function extracts the mutations by type and by sample and
writes them into the output directory separated into disease sub-directories.
The output can be processed by a separate call to mutationBatchAssess.}
\usage{
mutationBatchExtract(theMafDir, theTypeCountDir, theGDCflag,
thePar = "-Xmx8000m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theMafDir}{Standardized Data timestamp-version directory, such as
  GDC directory, which uses the path to mutations files. Alternatively, the
  DCC directory and its files such as path to mutations files can be used.}
  \item{theTypeCountDir}{The output directory for collated mutation type data.
  Collated data will be broken down into disease directories, and within the
  disease directory into files separated by the categories: platform, level,
  reference genome, and type of mutation. Not all categories exist within all
  data sets.}
  \item{theGDCflag}{TRUE handle as GDC data. FALSE handle as DCC data.}
  \item{thePar}{String arguments for rJava. Default -Xmx8000m}
}
\value{No useful value returned.}
\examples{
 \dontrun{
mutationBatchExtract(
theMafDir="/2016_01_20_1306/",
theTypeCountDir="DCCCountTsv",
theGDCflag=FALSE)
mutationBatchExtract(
theMafDir="/current/",
theTypeCountDir="GDCCountTsv",
theGDCflag=TRUE)
 }
}
