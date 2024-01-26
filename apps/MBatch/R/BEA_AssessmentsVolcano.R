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

volcano_openAndWriteIssuesLogFile<-function(theOutputDir)
{
  unlink(theOutputDir, recursive = TRUE, force = TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Unable to calculate volcano plot results\n", file=myFile, append=TRUE)
}

calcAndWriteVolcano <- function(theData, theTitle, theOutputDir, theLogFrameFlag,
                                theDataVersion, theTestVersion)
{
  logWarn("VOLCANO / FOLD CHANGE CODE IS NOT TESTED - DO NOT USE OR RELY ON THIS OUTPUT")
  logDebug("calcAndWriteVolcano - theOutputDir=", theOutputDir)
  logDebug("calcAndWriteVolcano - dim(theData)[1]=", dim(theData)[1])
  logDebug("calcAndWriteVolcano - dim(theData)[2]=", dim(theData)[2])
  the_sample_id_col <- colnames(theData@mBatches)[1]
  batch_type_list <- as.vector(unlist(colnames(theData@mBatches)[c(2:length(colnames(theData@mBatches)))]))
  for(batch_type in batch_type_list)
  {
    tryCatch({
      logDebug("calcAndWriteVolcano - start trycatch calcAndWriteVolcano")
      checkIfTestError()
      batch_df <- theData@mBatches[c("Sample", batch_type)]
      logDebug("getDSCwithExcerpt - import(mbatch.volcano.api)")
      calc <- import("mbatch.volcano.api")
      # the_title: str, the_sample_id_col: str
      # the_matrix: pandas.DataFrame, the_features: List[str], the_samples: List[str],
      # the_batches: pandas.DataFrame,
      # the_batch_types: List[str], the_output_dir
      # matrix and batches need to be data.frame for Python/Reticulate
      calc$volcano_calc_plot_from_r(as.character(theTitle), as.character(the_sample_id_col),
                                    r_to_py(theData@mData), as.vector(unlist(rownames(theData@mData))), as.vector(unlist(colnames(theData@mData))),
                                    r_to_py(batch_df),
                                    np_array(batch_type), as.character(theOutputDir),
                                    theLogFrameFlag,
                                    as.character(theDataVersion),
                                    as.character(theTestVersion))
      logDebug("calcAndWriteVolcano - end trycatch calcAndWriteVolcano")
    },error=function(myError)
    {
      logError("calcAndWriteVolcano(): error in calcAndWriteVolcano, ERROR= ", myError)
      batch_type_outdir <- cleanFilePath(theOutputDir, batch_type)
      batch_type_outdir <- addVersionsIfNeeded(batch_type_outdir, theDataVersion, theTestVersion)
      checkDirForCreation(batch_type_outdir)
      volcano_openAndWriteIssuesLogFile(batch_type_outdir)
    })
  }
}
