# MBatchUtils Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
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

buildSingleArchive <- function(theResultDir, theDataDir, theZipDir)
{
  ########
  message("buildSingleArchive::theResultDir=", theResultDir)
  message("buildSingleArchive::theDataDir=", theDataDir)
  message("buildSingleArchive::theZipDir=", theZipDir)
  ########
  # :param the_results_dir: directory with MBatch results
  # :param the_data_dir: directory with actual data
  # :param the_zip_dir: directory in which to place ZIP file
  # :param the_info_dir: full path to directory containing TEST_<version> labels
  # :param the_new_data: new data object for latest analysis
  # :param the_std_list: dictionary of data objects
  # :return: full pathname for ZIP file
  message("buildSingleArchive - import(mbatch.index.index_api)")
  calc <- import("mbatch.index.index_api")
  # the_results_dir: str,
  # the_data_dir: str,
  # the_zip_dir: str,
  # the_info_dir: str,
  # the_update_only_flag: bool,
  # the_new_data: typing.Optional[StandardizedData],
  # the_std_list: typing.Optional[List[StandardizedData]])
  zipFile <- calc$create_index_archive(theResultDir, theDataDir, theZipDir,
                                       file.path(theResultDir, "info"), FALSE, NULL, NULL)
  message("buildSingleArchive - after Python")
  zipFile
}

