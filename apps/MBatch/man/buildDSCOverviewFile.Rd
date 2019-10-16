\name{buildDSCOverviewFile}
\alias{buildDSCOverviewFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create DSC Text Output from DSC.RData Files
}
\description{
MBatch can create a human readable file containing all the DSC values from a run. These are built based on *DSC.RData files created during a call to mbatchAssess. The file uses the data's directory structure to create a table listing the existing data sets and the DSC results for each data set. This requires the data directory structure to match the download directories from the DCC for TCGA data, plus a top level "run date" directory. Additional directories are added to this to create places for the output.

For example, this is a valid directory structure for use with this output:

2011_11_30_0904\\GBM\\Expression-Genes\\BI__HT_HG-U133A\\Level_3\\original\\PCA\\PlateId\\ManyToMany\\PCA-Plus

This portion comes from the original data's directory structure: GBM\\Expression-Genes\\BI__HT_HG-U133A\\Level_3
Note that this matches the DCC Data Access Matrix directory structure for downloaded data.

The output file contains the following columns: "run-date", "disease-type", "data-type", "platform", "data-level", "correction-type", "PCA", "batch-type", "Overall-DSC", and "Overall-DSC-pvalue".

So, for the GBM output listed above, you end up with these values in the output, which makes it easier to track down interesting datasets using DSC values:

"run-date"	2011_11_30_0904			<- this is the date of the data run
"disease-type"	GBM					<- GBM disease data
"data-type"	Expression-Genes		<- type of data
"platform"	BI__HT_HG-U133A			<- platform
"data-level"	Level_3				<- level 3 DCC data
"correction-type"	original		<- this is "original" or uncorrected data. Assessments of corrected data will list the correction algorithm.
"PCA"	PCA							<- PCA assessment algorithm
"batch-type"	PlateId				<- assessment of the PlateId batch type

After running mbatchAssess, there are *DSC.RData files created -- these can either be use to create a DSC Overview file with buildDSCOverviewFile, or they can be deleted using clearDSCOverviewFiles. This is handled outside of these functions, so all data sets can be processed by multiple calls to mbatchAssess in a parallel or distributed environment, and then the DSC Overview created from these files at the end of the processing.
}
\usage{
buildDSCOverviewFile(theStartDir, theOutputFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theStartDir}{String indicating top of directory structure to search for *DSC.RData objects.}
  \item{theOutputFile}{String containing only the filename of output text file. Will be placed in theStartDir.}
}
\author{
Tod Casasent \email{tdcasasent@mdanderson.org}
Nianxiang Zhang \email{nzhang@mdanderson.org}
}
\seealso{
See also \code{\link{clearDSCOverviewFiles}}
}
\examples{
\dontrun{
example in doRunBEA_Files
}}
