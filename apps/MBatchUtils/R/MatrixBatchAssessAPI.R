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

mutBatchSingle <- function(theDataFile, theBatchFile, theTitle, theOutputDir,
                           theJavaArgs=c("-Xms24000m", "-Djava.awt.headless=true"), theThreads=5,
                           thePvalueCutoff=.00001, theZScoreCutoff=1.96,
                           theBatchTypes=c("BatchId", "PlateId", "ShipDate", "TSS"))
{
  message("mutBatchSingle")
  message("theDataFile=", theDataFile)
  message("theBatchFile=", theBatchFile)
  pValueVector <- loadAndRunForKruskal(theDataFile, theBatchFile, theBatchTypes, thePvalueCutoff, theZScoreCutoff)
  message("pValueVector=", paste(pValueVector, collapse=" "))
  print(pValueVector)
  message("theOutputDir=", theOutputDir)
  theOutputDir <- cleanFilePath(theOutputDir, "Discrete")
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
    outdir <- cleanFilePath(theOutputDir, myBatchType)
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
    makeKruskalPng(outfile, loggedPvalues, collectBatches, myBatchType, batchAddOns, thePvalueCutoff)
    makeKruskalTsv(outfileTSV, loggedPvalues, collectBatches, myBatchType, batchAddOns, thePvalueCutoff)
  }
}
