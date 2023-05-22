# MBatch Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

BeaRBN <- function(invMatrix, varMatrix, invReplicates, varReplicates, invGroupID = "",
                  varGroupID = "", matchedReplicates = TRUE, theCombineOnlyFlag = FALSE)
{
  # check that sufficient replicates are given for both matrices, warn otherwise
  if(length(invReplicates) == 0)
  {
    logWarn("No replicates provided for Invariant Matrix. Data will not be adjusted with RBN.")
  }
  else if(length(invReplicates) < 30)
  {
    logWarn("Less than 30 replicates provided for Invariant Matrix. RBN's performance may deteriorate.")
  }

  if(length(varReplicates) == 0)
  {
    logWarn("No replicates provided for Variant Matrix. Data will not be adjusted with RBN.")
  }
  else if(length(varReplicates) < 30)
  {
    logWarn("Less than 30 replicates provided for Variant Matrix. RBN's performance may deteriorate.")
  }

  # check that all replicates are present in both matrices
  stopifnotWithLogging(paste("The following invariant replicate name(s) were not found in the invariant matrix: ",
	 paste(invReplicates[!(invReplicates %in% rownames(invMatrix))], collapse=", "), sep=""),
	 all(invReplicates %in% rownames(invMatrix)))
  stopifnotWithLogging(paste("The following variant replicate name(s) were not found in the variant matrix: ",
	 paste(varReplicates[!(varReplicates %in% rownames(varMatrix))], collapse=", "), sep=""),
	 all(varReplicates %in% rownames(varMatrix)))

  # for matched replicates, the number of replicates in each matrix should be indentical
  if(matchedReplicates)
  {
    stopifnotWithLogging("Number of replicates from each matrix should be the same when matchedReplicates = TRUE.",
	    length(varReplicates) == length(invReplicates))
  }

  # get common features between both matrices, warning if none exist
  icols = colnames(invMatrix)
  vcols = colnames(varMatrix)
  common = icols[icols %in% vcols]
  warnifnotWithLogging("No common features found in both matrices", length(common) > 0)
  logInfo(paste("Found", length(common), "common features in both matrices."))
  # append group ID to row names of each matrix and save them, verifying that new row names are unique
  if(invGroupID != "")
  {
    irn = paste(rownames(invMatrix), invGroupID, sep="-")
  }
  else
  {
    irn = rownames(invMatrix)
  }

  if(varGroupID != "")
  {
    vrn = paste(rownames(varMatrix), varGroupID, sep="-")
  }
  else
  {
    vrn = rownames(varMatrix)
  }

  stopifnotWithLogging("Identical row names generated across both matrices. Please provide unique Group IDs.", !any(irn %in% vrn))

  if(FALSE==theCombineOnlyFlag)
  {
    # perform RBN adjustment on variant matrix only
    for(cc in common)
    {
      vm = varMatrix[varReplicates, cc]
      im = invMatrix[invReplicates, cc]
      if(matchedReplicates)
      {
        # if one replicate has NA value, set other replicate to NA as well
        vm = vm + (0 * im)
        im = im + (0 * vm)
      }

      # compute m and b using just the replicates with non-NA values
      m = sd(vm, na.rm = TRUE) / sd(im, na.rm = TRUE)
      b = mean(vm, na.rm = TRUE) - m * mean(im, na.rm = TRUE)

      # adjust the entire variant matrix column with parameters b and m
      if(!is.na(m) && !is.na(b))
      {
        varMatrix[, cc] = (varMatrix[, cc] - b) / m
      }
    }
  }

  # assign new row names to the matrices so the names are unique
  rownames(invMatrix) = irn
  rownames(varMatrix) = vrn

  # convert matrices to data frames so they can be merged properly
  invdataframe = data.frame(invMatrix, check.names=FALSE, stringsAsFactors=FALSE)
  vardataframe = data.frame(varMatrix, check.names=FALSE, stringsAsFactors=FALSE)

  # get column names that are only in one dataframe and insert in the other dataframe with NA values
  inames = colnames(invdataframe)
  vnames = colnames(vardataframe)
  notinvar = inames[!(inames %in% vnames)]
  notininv = vnames[!(vnames %in% inames)]
  vardataframe[,notinvar] = NA
  invdataframe[,notininv] = NA

  # concatenate both dataframes and return the result

  result = rbind(invdataframe, vardataframe, stringsAsFactors=FALSE)
  return(t(as.matrix(result)))
}
