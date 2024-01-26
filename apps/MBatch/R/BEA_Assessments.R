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


####################################################################
####################################################################
#
# generateFinalImageTitle<-function(theTitle, tempCentroids, thePathSubdir)
# {
# 	### cut off anything before the first occurance of thePathSubdir
# 	index <- regexpr(thePathSubdir, tempCentroids)
# 	resultString <- substr(tempCentroids,index-1, nchar(tempCentroids))
# 	### get rid of '_Diagram.png',
# 	resultString <- gsub("_Diagram.png", "", resultString, fixed=TRUE)
# 	### convert file underscores to spaces
# 	###resultString <- gsub("_", " ", resultString, fixed=TRUE)
# 	### add theTitle to replace what was removed
# 	resultString <- paste(theTitle, resultString, sep=" ")
# 	### convert directory separators to spaces
# 	###resultString <- gsub(.Platform$file.sep, " ", resultString, fixed=TRUE)
# 	### convert directory separators to space/space
# 	resultString <- gsub(.Platform$file.sep, paste(" ",.Platform$file.sep," ", sep=""), resultString, fixed=TRUE)
# 	### wrap to lengths
# 	resultString <- breakIntoTitle(resultString, theOldChar=" ", theWidth=50)
# 	return(resultString)
# }

####################################################################
####################################################################
