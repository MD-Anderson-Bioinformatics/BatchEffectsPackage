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

library(MBatchUtils)

if (!is.null(getTestOutputDir()))
{
  ########################################################
  ########################################################
  # writes to input directory, so copy files to output
  sourceDir=file.path(getTestInputDir(), "configout")
  outDir=file.path(getTestOutputDir(), "archiveout", "2018-07-11-1200")
  unlink(outDir)
  dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
  #file.copy(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"), file.path(outDir, "MBatchConfig.tsv"))
  #buildSingleArchive mbatchID, originalDataJsonFile, mbatchResultsDir, zipDir
  buildSingleArchive("tmp-test-id", NULL, sourceDir, outDir)
  TRUE
} else {
  message("No test data. Skip test.")
  TRUE
}
