# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

require(MBatch)

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()

theBatchFile=cleanFilePath(inputDir, "brca_rnaseq2_batches.tsv")
theBatchFile2=cleanFilePath(inputDir, "brca_agi4502_batches.tsv")
theOutputDir=cleanFilePath(cleanFilePath(outputDir, "ebnplus"), "EBNPlus_CombineBatches")
theCompareFile=cleanFilePath(compareDir, "EBNPlus_CombineBatches.tsv")
theBatchId1="RNASeqV2"
theBatchId2="Agilent4502"

if (!is.null(inputDir))
{
  unlink(theOutputDir, recursive=TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  dataBatches <- EBNPlus_CombineBatches(readAsGenericDataframe(theBatchFile), readAsGenericDataframe(theBatchFile2), theBatchId1, theBatchId2)
  writeAsGenericDataframe(cleanFilePath(theOutputDir, "BatchData.tsv"), dataBatches)
  compareDF <- readAsGenericDataframe(theCompareFile)
  print(all(dataBatches==compareDF, na.rm=TRUE))
  (all(dataBatches==compareDF, na.rm=TRUE))
} else {
  message("No test data. Skip test.")
  TRUE
}
