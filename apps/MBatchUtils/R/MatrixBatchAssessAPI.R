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

openAndWriteIssuesDiscrete<-function(theOutputDir, theMessage)
{
  filename <- cleanFilePath(theOutputDir, "error.log")
  message("openAndWriteIssuesDiscrete filename=", filename)
  myFile <- file(filename, "w+")
  on.exit(close(myFile))
  cat(theMessage, file=myFile, append=TRUE)
}

writeErrorsForAllBatchTypesDiscrete<-function(theMessage, theBatchTypes, theOutputDir,
                                              theTitle, theDataVersion, theTestVersion)
{
  unlink(theOutputDir, recursive = TRUE, force = TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  for(batchType in theBatchTypes)
  {
    ### compile data and information for display
    message("writeErrorsForAllBatchTypesDiscrete batchType=", batchType)
    myBatchDir <- file.path(theOutputDir, batchType)
    myBatchDir <- addVersionsIfNeeded(myBatchDir, theDataVersion, theTestVersion)
    message("writeErrorsForAllBatchTypesDiscrete myBatchDir=", myBatchDir)
    # do errors for each batch type
    dir.create(myBatchDir, showWarnings=FALSE, recursive=TRUE)
    mesg=paste("No Kruskal-Wallis/Dunn's Test -", theMessage, " - ", theTitle, sep="")
    openAndWriteIssuesDiscrete(myBatchDir, mesg)
  }
}

mutBatchSingle <- function(theMBatchData, theTitle,
                           theOutputDir, theDataVersion, theTestVersion,
                           theThreads=5, thePvalueCutoff=.00001, theZScoreCutoff=1.96,
                           theBatchTypes=c("BatchId", "PlateId", "ShipDate", "TSS"))
{
  message("mutBatchSingle")
  origTitle <- theTitle
  message("origTitle=", origTitle)
  message("theOutputDir=", theOutputDir)
  baseDir <- theOutputDir
  theOutputDir <- cleanFilePath(theOutputDir, "Discrete")
  message("updated theOutputDir=", theOutputDir)
  tryCatch({
    checkIfTestError()
    pValueVector <- loadAndRunForKruskal(theMBatchData, theBatchTypes, thePvalueCutoff, theZScoreCutoff)
    message("pValueVector=", paste(pValueVector, collapse=" "))
    print(pValueVector)
    ####################
    # add p-values to result list
    resultList <- c(list(pValueVector))
    # add title to title vector
    titleVector <- c(theTitle)
    for (myBatchType in theBatchTypes)
    {
      # for each batch type
      message("Do Kruskal output for: ", myBatchType)
      # extract and convert p-values into logged p-values (-log10) for the given batch type
      loggedPvalues <- collectLoggedPvalues(resultList, myBatchType)
      # if the values for a batch type has batch names (in addition to the p-value)
      # collect the batch names found as significant by a Dunn's test
      collectBatches <- collectLoggedBatches(resultList, myBatchType)
      # subdir named for batch type and count type
      outdir <- cleanFilePath(theOutputDir, myBatchType)
      outdir <- addVersionsIfNeeded(outdir, theDataVersion, theTestVersion)
      # make sure sub-dir exists
      dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
      # output for Kruval - name the file "FullMutCounts_<batch-type>_<mutation-type>_Diagram.PNG"
      outfile <- cleanFilePath(outdir, "KW_Dunns_Diagram.PNG")
      outfileTSV <- cleanFilePath(outdir, "KW_Dunns_Diagram.tsv")
      # add names from title vector to logged p-values
      names(loggedPvalues) <- titleVector
      # if Dunn's test returned batch names, put them into parens to put in diagram otherwise, put an empty string
      batchAddOns <- getBatchAddOn(collectBatches)
      # write the PNG for Kruskal's output
      mytitle <- makeKruskalPng(outfile, loggedPvalues, collectBatches, myBatchType, batchAddOns, thePvalueCutoff)
      datasetTitle <- gsub(baseDir, "", titleVector)
      if (!startsWith(datasetTitle, "/"))
      {
        datasetTitle <- paste("/", datasetTitle, sep="")
      }
      names(loggedPvalues) <- datasetTitle
      saveTitle <- paste(origTitle, " / ", myBatchType, " / ", "Kruskal-Wallis Dunns-Test", sep="")
      makeKruskalTsv(outfileTSV, loggedPvalues, collectBatches, myBatchType, batchAddOns, thePvalueCutoff, saveTitle, outfile)
    }
  },
  # don't catch warnings
  #warning=function(e)
  #{
  #  message("1 Unable to generate Kruskal-Wallis (discrete)")
  #  writeErrorsForAllBatchTypesDiscrete("Kruskal-Wallis/Dunn's Test Failed",
  #                                      theBatchTypes, theOutputDir,
  #                                      theTitle, theDataVersion, theTestVersion)
  #},
  error=function(e)
  {
    message("2 Unable to generate Kruskal-Wallis (discrete)")
    writeErrorsForAllBatchTypesDiscrete("Kruskal-Wallis/Dunn's Test Failed",
                                        theBatchTypes, theOutputDir,
                                        theTitle, theDataVersion, theTestVersion)
  })
}
