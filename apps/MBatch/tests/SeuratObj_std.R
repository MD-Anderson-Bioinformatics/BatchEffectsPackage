# MBatch Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
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

compMatrixFile=cleanFilePath(compareDir, "seurat_matrix.tsv")
compBatchesFile=cleanFilePath(compareDir, "seurat_batches.tsv")

if (!is.null(compareDir))
{
  warnLevel<-getOption("warn")
  on.exit(options(warn=warnLevel))
  # warnings are errors
  options(warn=3)
  # if there is a warning, show the calls leading up to it
  options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  options(showErrorCalls=TRUE)
  # convert Standardized Data files to Seurat object
  theAssayToUse <- "RNA"
  theDataSlotToUse <- "counts"
  # theSeuratObj, theMatrixFilePath, theBatchFilePath, theAssayToUse, theDataSlotToUse
  seuratObj <- updateOrBuildSeuratObjectFromStdDataFiles(NULL,
                                                         compMatrixFile, compBatchesFile,
                                                         theAssayToUse, theDataSlotToUse)
  print(seuratObj)
  TRUE
} else {
  message("No test data. Skip test.")
  TRUE
}