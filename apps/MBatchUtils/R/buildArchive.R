# MBatchUtils Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

library(MBatch)

#############################################################################
## exported
#############################################################################

# the SourceDir should contain a single directory, the data for which to build an archive
buildSingleArchive <- function(mbatchID, originalDataJsonFile, mbatchResultsDir, zipDir)
{
  # runBEVIndex(mbatchID, versionStamp, versionType, originalDataJsonFile, mbatchResultsDir, zipDir)
  ########
  message("buildSingleArchive::mbatchID=", mbatchID)
  message("buildSingleArchive::originalDataJsonFile=", originalDataJsonFile)
  message("buildSingleArchive::mbatchResultsDir=", mbatchResultsDir)
  message("buildSingleArchive::zipDir=", zipDir)
  ########
  # jinit for BEVIndex
  myJavaJars <- file.path(system.file("BEVIndex", "BEVIndex.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "commons-codec-1.9.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "commons-io-2.5.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "commons-lang3-3.4.jar", package="MBatchUtils"),
                          system.file("BEVIndex", "gson-2.7.jar", package="MBatchUtils"),
                          fsep=.Platform$path.sep)
  .jinit(classpath=myJavaJars, force.init = TRUE, parameters=c("-Xms4800m", "-Djava.awt.headless=true"))
  ########
  .jcall("edu/mda/bioinfo/bevindex/BEVIndex", returnSig = "V",
         method='runBEVIndex',
         .jnew("java/lang/String",mbatchID),
         .jnew("java/lang/String",originalDataJsonFile),
         .jnew("java/lang/String",mbatchResultsDir),
         .jnew("java/lang/String",zipDir))
}

