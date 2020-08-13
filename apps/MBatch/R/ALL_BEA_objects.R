# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

library(methods)

###Samples
###setClass("foo", representation(a = "character", b = "numeric"))
###setClass("bar", representation(d = "numeric", c = "numeric"))
###setClass("baz", contains = c("foo", "bar"))
###setMethod("initialize", "xx",
###		function(.Object, b)
###		{
###			.Object@b <- b
###			.Object@a <- nchar(b)
###			.Object
###		})

setClass("BEA_DATA", representation(
	mData="matrix",
	mBatches="data.frame",
	mCovariates="data.frame"
))

setMethod("initialize", "BEA_DATA",
					function(.Object,
									 theData, theBatches, theCovariates
					)
					{
						.Object@mData <- theData
						.Object@mBatches <- theBatches
						.Object@mCovariates <- theCovariates
						.Object
					})

setClass("Corrections_EBNP", representation(
	mEBNP_DoCorrectionFlag="logical",
	mEBNP_Data2="BEA_DATA",
	mEBNP_TrimBarcodesFunction="function",
	mEBNP_TrimGenesFunction="function",
	mEBNP_RemoveRowDuplicatesFunction="function",
	mEBNP_RemoveColDuplicatesFunction="function",
	mEBNP_Data1BatchId="character",
	mEBNP_Data2BatchId="character",
	mEBNP_BatchWithZero="character",
	mEBNP_FixDataSet="numeric",
	mEBNP_CorrectForZero="logical",
	mEBNP_ValidationRatio="numeric"
))


setClass("Logging", representation(
  mFile="character",
  mLevelNamesToLog="vector",
  mLevelNames="vector",
  mSeparator="character",
  mConsole="logical"
))

setMethod("initialize", "Logging",
          function(.Object,
                   theFile="",
                   theLevelNamesToLog=c('DEBUG', 'TIMING', 'INFO', 'WARN', 'PERCENT', 'ERROR'),
                   theLevelNames=     c('DEBUG', 'TIMING', 'INFO', 'WARN', 'PERCENT', 'ERROR'),
                   theSeparator=" ",
                   theConsole=TRUE)
          {
            .Object@mFile <- theFile
            .Object@mLevelNamesToLog <- theLevelNamesToLog
            .Object@mLevelNames <- theLevelNames
            .Object@mSeparator <- theSeparator
            .Object@mConsole <- theConsole
            .Object
          })

