
mutBatchSingle <- function(theDataFile, theBatchFile, theTitle, theOutputDir,
                           theJavaArgs=c("-Xms24000m", "-Djava.awt.headless=true"), theThreads=5,
                           thePvalueCutoff=.00001, theZScoreCutoff=1.96, thePCAflag=FALSE,
                           theBatchTypes=c("BatchId", "PlateId", "ShipDate", "TSS"))
{
  message("mutBatchSingle")
  message("theDataFile=", theDataFile)
  message("theBatchFile=", theBatchFile)
  pValueVector <- loadAndRunForKruskal(theDataFile, theBatchFile, theBatchTypes, thePvalueCutoff, theZScoreCutoff)
  message("pValueVector=", paste(pValueVector, collapse=" "))
  print(pValueVector)
  message("theOutputDir=", theOutputDir)
  theOutputDir <- file.path(theOutputDir, "MutBatch")
  message("updated theOutputDir=", theOutputDir)
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
    outdir <- file.path(theOutputDir, myBatchType)
    # make sure sub-dir exists
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    # output for Kruval - name the file "FullMutCounts_<batch-type>_<mutation-type>_Diagram.PNG"
    outfile <- file.path(outdir, "KW_Dunns_Diagram.PNG")
    outfileTSV <- file.path(outdir, "KW_Dunns_Diagram.tsv")
    # add names from title vector to logged p-values
    names(loggedPvalues) <- titleVector
    # if Dunn's test returned batch names, put them into parens to put in diagram otherwise, put an empty string
    batchAddOns <- getBatchAddOn(collectBatches)
    # write the PNG for Kruskal's output
    makeKruskalPng(outfile, loggedPvalues, collectBatches, myBatchType, batchAddOns, thePvalueCutoff)
    makeKruskalTsv(outfileTSV, loggedPvalues, collectBatches, myBatchType, batchAddOns, thePvalueCutoff)
  }
}
