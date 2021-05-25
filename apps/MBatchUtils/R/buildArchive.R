# MBatchUtils Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

####################################################################
###
####################################################################

getJarsFromDir <- function(theDir)
{
  allJars <- dir(theDir, "*.jar", full.names=TRUE)
  paste(allJars, collapse=.Platform$path.sep)
}

####################################################################
###
####################################################################

#############################################################################
## exported
#############################################################################

# the SourceDir should contain a single directory, the data for which to build an archive
buildSingleArchive <- function(theMbatchID, theResultDir, theDataDir, theZipDir)
{
  ########
  message("buildSingleArchive::theMbatchID=", theMbatchID)
  message("buildSingleArchive::theResultDir=", theResultDir)
  message("buildSingleArchive::theDataDir=", theDataDir)
  message("buildSingleArchive::theZipDir=", theZipDir)
  ########
  # jinit for BEVIndex
  myJavaJars <- getJarsFromDir(dirname(system.file("BEVIndex", "BEVIndex.jar", package="MBatchUtils")))
  .jinit(classpath=myJavaJars, force.init = TRUE, parameters=updateJavaParameters(c("-Xms4800m", "-Djava.awt.headless=true")))
  ########
  # String theMbatchID, String theResultDir, String theDataDir, String theZipDir
  .jcall("edu/mda/bcb/bevindex/BEVIndex", returnSig = "V",
         method='runBEVIndex',
         .jnew("java/lang/String",theMbatchID),
         .jnew("java/lang/String",theResultDir),
         .jnew("java/lang/String",theDataDir),
         .jnew("java/lang/String",theZipDir))
}

