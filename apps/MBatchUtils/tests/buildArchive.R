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

require(MBatchUtils)

if (!is.null(getTestOutputDir()))
{
  ########################################################
  ########################################################
  # writes to input directory, so copy files to output
  sourceDir=cleanFilePath(getTestInputDir(), "configout")
  outDir=cleanFilePath(getTestOutputDir(), "archiveout")
  theOutputDirResults=cleanFilePath(outDir, "ZIP-RESULTS")
  theOutputDirData=cleanFilePath(outDir, "ZIP-DATA")
  unlink(outDir, recursive = TRUE)
  dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
  dir.create(theOutputDirResults, recursive=TRUE, showWarnings=FALSE)
  dir.create(theOutputDirData, recursive=TRUE, showWarnings=FALSE)
  file.copy(cleanFilePath(sourceDir, "ZIP-RESULTS"), outDir, recursive=TRUE)
  file.copy(cleanFilePath(sourceDir, "ZIP-DATA"), outDir, recursive=TRUE)
  buildSingleArchive(theOutputDirResults, theOutputDirData, outDir)
  TRUE
} else {
  message("No test data. Skip test.")
  TRUE
}
