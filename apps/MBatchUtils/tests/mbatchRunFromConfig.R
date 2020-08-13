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
  #warnLevel<-getOption("warn")
  #on.exit(options(warn=warnLevel))
  # warnings are errors
  #options(warn=3)
  # if there is a warning, show the calls leading up to it
  #options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  #options(showErrorCalls=TRUE)
  ########################################################
  ########################################################
  # writes to input directory, so copy files to output
  outDir=file.path(getTestOutputDir(), "configout", "2018-07-11-1200")
  print(outDir)
  unlink(outDir, recursive=TRUE)
  dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
  print(file.exists(outDir))
  print(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"))
  print(file.path(outDir, "MBatchConfig.tsv"))
  file.copy(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"), file.path(outDir, "MBatchConfig.tsv"))
  print(file.exists(file.path(outDir, "MBatchConfig.tsv")))
  file.copy(file.path(getTestInputDir(), "config", "matrix_data.tsv"), file.path(outDir, "matrix_data.tsv"))
  file.copy(file.path(getTestInputDir(), "config", "batches.tsv"), file.path(outDir, "batches.tsv"))
  configFile=file.path(outDir, "MBatchConfig.tsv")
  ########################################################
  ########################################################
  mbatchRunFromConfig(theConfigFile=configFile, theOutputDir=outDir, theNaStrings="NA")
  file.exists(file.path(outDir, "MBATCH_SUCCESS.txt"))
} else {
  message("No test data. Skip test.")
  TRUE
}
