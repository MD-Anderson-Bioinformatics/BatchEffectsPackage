## ---- echo=TRUE----------------------------------------------------------
{
  library(MBatch)
  
  # set the paths
  theDataFile1="/bea_testing/MATRIX_DATA/brca_rnaseq2_matrix_data.tsv"
  theDataFile2="/bea_testing/MATRIX_DATA/brca_agi4502_matrix_data.tsv"

  # trim genes to get just gene symbols from standardized data
  trimGenes <- function(theGenes)
  {
    foo <- as.vector(unlist(
      sapply(theGenes, function(theGene)
      {
        # keep the same if it starts with ?
        if (TRUE==grepl("^[?]+", theGene))
        {
          return(theGene)
        }
        else
        {
          # split on the | and take the first argument
          # this makes no change if no pipe
          return(strsplit(theGene, "|", fixed=TRUE)[[1]][1])
        }
      })
    ))
    foo
  }
  
  # remove duplicates from columns (samples)
  removeDuplicatesFromColumns <- function(theMatrix)
  {
    indexOfDuplicates <- which(duplicated(colnames(theMatrix)))
    if (length(indexOfDuplicates) > 0)
    {
      # minus sign uses inverse of indexes
      theMatrix <- theMatrix[ ,-indexOfDuplicates]
    }
    return(theMatrix)
  }
  
  # remove duplicates from rows (genes/probes)
  removeDuplicatesFromRows <- function(theMatrix)
  {
    indexOfDuplicates <- which(duplicated(rownames(theMatrix)))
    if (length(indexOfDuplicates) > 0)
    {
      # minus sign uses inverse of indexes
      theMatrix <- theMatrix[-indexOfDuplicates, ]
    }
    return(theMatrix)
  }

  if ((!dir.exists(theDataFile1))&&(!dir.exists(theDataFile2)))
  {
    warnLevel<-getOption("warn")
    on.exit(options(warn=warnLevel))
    # warnings are errors
    options(warn=3)
    # if there is a warning, show the calls leading up to it
    options(showWarnCalls=TRUE)
    # if there is an error, show the calls leading up to it
    options(showErrorCalls=TRUE)
    #
    # read the files in. This can be done however you want
    theDataMatrix1 <- readAsGenericMatrix(theDataFile1)
    theDataMatrix2 <- readAsGenericMatrix(theDataFile2)
    # this is the reduce genes to just gene symbols, handling those from standardized data
    rownames(theDataMatrix1) <- trimGenes(rownames(theDataMatrix1))
    rownames(theDataMatrix2) <- trimGenes(rownames(theDataMatrix2))
    # remove any duplicates (this is a requirement for EBNplus)
    theDataMatrix1 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(theDataMatrix1))
    theDataMatrix2 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(theDataMatrix2))
    print("Is this data acceptable?")
    EBNPlus_CheckData_Structures(theDataMatrix1, theDataMatrix2)
    print("If you see this, it is.")
  }
}

