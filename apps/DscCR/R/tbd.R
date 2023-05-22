#DscCR Copyright 2023 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

setClass("DSCC_RESULT", representation(
  mListOfDSCbyGene="vector",
  mListOfDWbyGene="vector",
  mListOfDBbyGene="vector",
  mDSC="numeric",
  mDB="numeric",
  mDW="numeric"))

buildList<-function(theValue, theLength)
{
  ### looks weird so that it works with NA and NAN
  return(as.vector(unlist(lapply(c(1:theLength), function(x) {return(theValue)}))))
}

createDsccResult <- function(theLength)
{
  results <- new("DSCC_RESULT")
  results@mListOfDSCbyGene <- buildList(0.0, theLength)
  results@mListOfDWbyGene <- buildList(0.0, theLength)
  results@mListOfDBbyGene <- buildList(0.0, theLength)
  results@mDSC <- 0.0
  results@mDB <- 0.0
  results@mDW <- 0.0
  return(results)
}

dsccDoDscPerms <- function (theValues, theDim, theBatchIds, thePerms, theThreads)
{
  ###message("thePerms = ", thePerms)
  ###message("theThreads = ", theThreads)
  uniqueIds <- sort(unique(sort(theBatchIds)))
  ###message("NaN = ", sum(is.nan(as.double(theValues))))
  resultsDSCbyGene <- buildList(0.0, (thePerms*length(theBatchIds)))
  resultsDBbyGene <- buildList(0.0, (thePerms*length(theBatchIds)))
  resultsDWbyGene <- buildList(0.0, (thePerms*length(theBatchIds)))
  resultsDSC <- buildList(0.0, thePerms)
  resultsDB <- buildList(0.0, thePerms)
  resultsDW <- buildList(0.0, thePerms)
  foo <- .C("doDscPerms", as.double(theValues), as.integer(theDim), as.character(theBatchIds), as.character(uniqueIds), as.integer(length(uniqueIds)),
            as.integer(thePerms), as.integer(theThreads),
            as.double(resultsDSCbyGene), as.double(resultsDBbyGene), as.double(resultsDWbyGene), as.double(resultsDSC), as.double(resultsDB), as.double(resultsDW));
  resultsDSCbyGene <- foo[[8]]
  resultsDBbyGene <- foo[[9]]
  resultsDWbyGene <- foo[[10]]
  resultsDSC <- foo[[11]]
  resultsDB <- foo[[12]]
  resultsDW <- foo[[13]]
  ###message("Perm Results arg resultsDSC =", paste(resultsDSC,collapse=","))
  ###message("Perm Results arg resultsDB =", paste(resultsDB,collapse=","))
  ###message("Perm Results arg resultsDW =", paste(resultsDW,collapse=","))
  results <- c(1:thePerms)
  results <- sapply(results, function(x)
  {
    arrayListStart <- (1+(x-1)*length(theBatchIds))
    arrayListStop <- (x*length(theBatchIds))
    dsccResultObj <- createDsccResult(length(theBatchIds))
    dsccResultObj@mListOfDSCbyGene <- resultsDSCbyGene[arrayListStart:arrayListStop]
    dsccResultObj@mListOfDBbyGene <- resultsDBbyGene[arrayListStart:arrayListStop]
    dsccResultObj@mListOfDWbyGene <- resultsDWbyGene[arrayListStart:arrayListStop]
    dsccResultObj@mDSC <- resultsDSC[x]
    dsccResultObj@mDB <- resultsDB[x]
    dsccResultObj@mDW <- resultsDW[x]
    ###message("Perm ", x, " results.mDSC =" , resultsDSC[x], " results.mDB =" , resultsDB[x], " results.mDW =" , resultsDW[x]);
    return(dsccResultObj)
  })
  return(results)
}

dsccGetDSCwithExcerpt <- function (theValues, theDim, theBatchIds)
{
  dsccResultObj <- createDsccResult(length(theBatchIds))
  uniqueIds <- sort(unique(sort(theBatchIds)))
  ###message("length(theBatchIds) = ", length(theBatchIds))
  ###message("length(uniqueIds) = ", length(uniqueIds))
  ###message("theDim[1] number of samples = ", theDim[[1]])
  ###message("theDim[2] number of batch ids = ", theDim[[2]])
  foo <- .C("getDSCwithExcerpt", as.double(theValues), as.integer(theDim), as.character(theBatchIds), as.character(uniqueIds), as.integer(length(uniqueIds)),
            as.double(dsccResultObj@mListOfDSCbyGene), as.double(dsccResultObj@mListOfDBbyGene), as.double(dsccResultObj@mListOfDWbyGene),
            as.double(dsccResultObj@mDSC), as.double(dsccResultObj@mDB), as.double(dsccResultObj@mDW)
  )
  dsccResultObj@mListOfDSCbyGene <- foo[[6]]
  dsccResultObj@mListOfDBbyGene <- foo[[7]]
  dsccResultObj@mListOfDWbyGene <- foo[[8]]
  dsccResultObj@mDSC <- foo[[9]]
  dsccResultObj@mDB <- foo[[10]]
  dsccResultObj@mDW <- foo[[11]]
  ###message("dsccGetDSCwithExcerpt Results arg dsccResultObj@mListOfDSCbyGene =", paste(dsccResultObj@mListOfDSCbyGene,collapse=","))
  ###message("dsccGetDSCwithExcerpt Results arg dsccResultObj@mListOfDBbyGene =", paste(dsccResultObj@mListOfDBbyGene,collapse=","))
  ###message("dsccGetDSCwithExcerpt Results arg dsccResultObj@mListOfDWbyGene =", paste(dsccResultObj@mListOfDWbyGene,collapse=","))
  ###message("dsccGetDSCwithExcerpt Results arg dsccResultObj@mDSC =", paste(dsccResultObj@mDSC,collapse=","))
  ###message("dsccGetDSCwithExcerpt Results arg dsccResultObj@mDB =", paste(dsccResultObj@mDB,collapse=","))
  ###message("dsccGetDSCwithExcerpt Results arg dsccResultObj@mDW =", paste(dsccResultObj@mDW,collapse=","))
  return(dsccResultObj)
}

dsccTest <- function(theDouble, theInt, theString, theDoubleList, theIntList, theStringList, theMatrix)
{
  .C("dsccTest", as.double(theDouble), as.integer(theInt), as.character(theString),
     as.double(theDoubleList), as.integer(length(theDoubleList)),
     as.integer(theIntList), as.integer(length(theIntList)),
     as.character(theStringList), as.integer(length(theStringList)),
     as.double(as.vector(theMatrix)), as.integer(dim(theMatrix))
  );
}
