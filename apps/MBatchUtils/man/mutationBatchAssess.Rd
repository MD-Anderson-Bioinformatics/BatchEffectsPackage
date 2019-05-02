\name{mutationBatchAssess}
\alias{mutationBatchAssess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mutationBatchAssess}
\description{This function processes the directories and files created by mutationBatchExtract
and does batch assessment processing on mutation data. Batch comparisons are both
between platforms, within platforms, and between disease types.}
\usage{
mutationBatchAssess(theTypeCountDir, theOutputDir, theJavaArgs = c("-Xms24000m", "-Djava.awt.headless=true"),
theThreads = 5, thePvalueCutoff = 1e-05, theZScoreCutoff = 1.96,
thePCAflag = FALSE, theBatchTypes = c("BatchId", "PlateId", "ShipDate", "TSS"),
theMutationTypes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theTypeCountDir}{full path to top level of directory containing disease sub-directories and files created by mutationBatchExtract}
  \item{theOutputDir}{full path to output directory for this function}
  \item{theJavaArgs}{Default: "-Xms24000m -Djava.awt.headless=true" - argument to Java for Boxplot and DSC computations.}
  \item{theThreads}{Default: 5 - argument to Java for # cores (threads) to use for DSC computations}
  \item{thePvalueCutoff}{Default: .00001 - P-value to use for cutoff for "significant" difference in Kruskal-Wallis and Dunn's tests}
  \item{theZScoreCutoff}{Default 1.96 - cut off for Dunn's test for "significant" confidence interval}
  \item{thePCAflag}{Default: FALSE - TRUE means to PCA with DSC (which is slow), FALSE means don't do PCA}
  \item{theBatchTypes}{Default: c("BatchId", "PlateId", "ShipDate", "TSS") - Batch types to check (correspond to columns for batch files)}
  \item{theMutationTypes}{Default: NULL. If not NULL, a vector of mutation types to process, such as "Total"}
}
\details{
This function processes the directories and files created by mutationBatchExtract
and does batch assessment processing on mutation data. Batch comparisons are both
between platforms, within platforms, and between disease types.
Output is divided into <BatchType>_<MutationType> subdirectories. Datasets without Kruskal-Wallis-called significant results
will ony have a PNG (FullMutCounts_<BatchType>_<MutationType>_Diagram.PNG) giving the Kruskal-Wallis/Dunn's Test results.
For datasets with significant calls, the aforementioned PNG will exist, as will files using the patterns
FullMutCounts_<BatchType>_<MutationType>_<DiseaseType>_MutDots_Diagram.PNG
		which contains all samples sorted by batch with mutations for all platforms plotted
NarrowBoxplot_<BatchType>_<MutationType>_<DiseaseType>_<Log10|ZScore>_<Platform_GenomeReference_MutationType>_Diagram.PNG
		which contains mutations for the given platform plotted by batches (with more than the average number of samples)
		and plots are provided both in terms of Log10 and ZScores.
WideBoxplot_<BatchType>_<MutationType>_<DiseaseType>_<Log10|ZScore>_Diagram.PNG
		which contains mutations for all platforms plotted by batches and plots are provided both in terms of Log10 and ZScores.
callReference.tsv
		is a tab delimited file with a header of MutationType, BatchTuype, MutationFile, and called Batches.
    Batches are comma delimited within parenthesis.

For Narrow and Wide Boxplots, the vertical axis is the number of mutations per sample, with the appropriate (Z-Score or Log10(+1)) modification/normalization. For MutDots plots, the vertical axis is the number of mutations per sample with a log10(+1) modification/normalization. Plotted points are features.
}
\value{Returns no useful value}
\references{}
\author{}
\note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
