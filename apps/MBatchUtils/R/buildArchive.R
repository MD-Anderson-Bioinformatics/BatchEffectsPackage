#MBatchUtils Copyright ? 2018 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(MBatch)

#############################################################################
## internal
#############################################################################

buildZip <- function(theSourceDir, theZip)
{
  on.exit(setwd(getwd()), add=TRUE)
  setwd(theSourceDir)
  myFiles <- list.files(".", all.files=TRUE, recursive=TRUE, include.dirs=TRUE)
  print(myFiles)
  zip(theZip, myFiles)

}

buildInternalIndex <- function(theDataRunDir, theZipPath, theDataSetName, theDataSetLabel, theDefaultLink,
                               theLevelLabels, theTooltipFile)
{
  ########
  message("buildInternalIndex::theDataRunDir=", theDataRunDir)
  message("buildInternalIndex::theZipPath=", theZipPath)
  message("buildInternalIndex::theDataSetName=", theDataSetName)
  message("buildInternalIndex::theDataSetLabel=", theDataSetLabel)
  message("buildInternalIndex::theDefaultLink=", theDefaultLink)
  message("buildInternalIndex::theLevelLabels=", theLevelLabels)
  message("buildInternalIndex::theTooltipFile=", theTooltipFile)
  ########
  .jcall("edu/mda/bioinfo/bevindex/BEVIndex", returnSig = "V",
         method='addIndexAndZipDataSet',
         .jnew("java/lang/String",theDataRunDir),
         .jnew("java/lang/String",theZipPath),
         .jarray(as.vector(as.character(theDefaultLink))),
         .jnew("java/lang/String",theDataSetName),
         .jnew("java/lang/String",theDataSetLabel),
         .jarray(as.vector(as.character(theLevelLabels))),
         .jnew("java/lang/String",theTooltipFile))
}

buildExternalIndex <- function(theDataRunDir, dataRunDirs, theDefaultLink, theDataSetName, theDataSetLabel,
                               theExternalIndexPath, theLevelLabels, theTooltipFile)
{
  ########
  message("buildExternalIndex::theDataRunDir=", theDataRunDir)
  message("buildExternalIndex::dataRunDirs=", dataRunDirs)
  message("buildExternalIndex::theDefaultLink=", theDefaultLink)
  message("buildExternalIndex::theDataSetName=", theDataSetName)
  message("buildExternalIndex::theDataSetLabel=", theDataSetLabel)
  message("buildExternalIndex::theExternalIndexPath=", theExternalIndexPath)
  message("buildExternalIndex::theLevelLabels=", theLevelLabels)
  message("buildExternalIndex::theTooltipFile=", theTooltipFile)
  ########
  .jcall("edu/mda/bioinfo/bevindex/BEVIndex", returnSig = "V",
         method='addExternalIndexToDataRun',
         .jnew("java/lang/String",theDataRunDir),
         .jarray(as.vector(as.character(dataRunDirs))),
         .jarray(as.vector(as.character(theDefaultLink))),
         .jarray(as.vector(as.character(theDataSetName))),
         .jarray(as.vector(as.character(theDataSetLabel))),
         .jarray(as.vector(as.character(theExternalIndexPath))),
         .jarray(as.vector(as.character(theLevelLabels))),
         .jnew("java/lang/String",theTooltipFile))
}

getLinks <- function(theSourceDir, theDefaultLink)
{
  files <- list.files(theSourceDir, theDefaultLink, recursive=TRUE, include.dirs=TRUE)[1]
  strsplit(files, .Platform$file.sep, fixed=TRUE)
}

#############################################################################
## exported
#############################################################################

# the SourceDir should contain a single directory, the data for which to build an archive
buildSingleArchive <- function(theSourceDir, theArchiveDir, theDataRunDir, theExternalIndexPath,
                               theDataSetName, theDataSetLabel,
                         theDefaultLink = "PCA",
                         theLevelLabels = c("Run Options", "Algorithm", "Diagram Type", "Sub-Type"),
                         theTooltipFile=system.file("BEVIndex", "tooltip_gdc.tsv", package="MBatchUtils"))
{
  ########
  # zip file path
  zipFile <- file.path(theArchiveDir, "ResultSet.zip")
  # data set names
  dataSetDirs <- basename(list.dirs(theSourceDir, recursive = FALSE)[1])
  # get link path to theDefaultLink
  defaultLinkVector <- getLinks(theSourceDir, theDefaultLink)
  ########
  # build the zip file
  buildZip(theSourceDir, zipFile)
  ########
  message("buildSingleArchive::theSourceDir=", theSourceDir)
  message("buildSingleArchive::theArchiveDir=", theArchiveDir)
  message("buildSingleArchive::theDataRunDir=", theDataRunDir)
  message("buildSingleArchive::theExternalIndexPath=", theExternalIndexPath)
  message("buildSingleArchive::theDataSetName=", theDataSetName)
  message("buildSingleArchive::theDataSetLabel=", theDataSetLabel)
  message("buildSingleArchive::zipFile=", zipFile)
  message("buildSingleArchive::dataSetDirs=", dataSetDirs)
  message("buildSingleArchive::defaultLinkVector=", defaultLinkVector)
  ########
  # jinit for BEVIndex
  myJavaJars <- file.path(system.file("BEVIndex", "BEVIndex.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "commons-codec-1.9.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "commons-io-2.5.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "commons-lang3-3.4.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "gson-2.7.jar", package="MBatchUtils"),
                          fsep=.Platform$path.sep)
  .jinit(classpath=myJavaJars, force.init = TRUE, parameters="-Xms4800m")
  ########
  # build the internal index
  buildInternalIndex(theDataRunDir, zipFile, theDataSetName, theDataSetLabel, defaultLinkVector, theLevelLabels, theTooltipFile)
  ########
  # build the external index
  buildExternalIndex(theDataRunDir, dataSetDirs, defaultLinkVector, theDataSetName, theDataSetLabel,
                     theExternalIndexPath, theLevelLabels, theTooltipFile)
}

#############################################################################
## exported
#############################################################################

# the SourceDir should contain a single directory, the data for which to build an archive
buildRunArchives <- function(theRunDir, theExternalIndexPath, theFinalPath,
                               theDataSetName, theDataSetLabel,
                               theDefaultLink = "PCA",
                               theLevelLabels = c("Data Release", "Program", "Disease", "Workflow", "Algorithm", "Diagram Type", "Sub-Type"),
                               theTooltipFile=system.file("BEVIndex", "tooltip_gdc.tsv", package="MBatchUtils"))
{

}
