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

require(MBatchUtils)
require(tools)

if (!is.null(getTestOutputDir()))
{
  warnLevel<-getOption("warn")
  on.exit(options(warn=warnLevel))
  # warnings are errors
  options(warn=3)
  # if there is a warning, show the calls leading up to it
  options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  options(showErrorCalls=TRUE)
  ########################################################
  ########################################################
  baseTestDir=getTestInputDir()
  baseOutputDir=getTestOutputDir()
  dccMinMut=cleanFilePath(baseTestDir, "DCCMin_MUT")
  gdcMinMut=cleanFilePath(baseTestDir, "GDCMin_MUT")
  dccMinOut=cleanFilePath(baseOutputDir, "DCCMin_OUT")
  gdcMinOut=cleanFilePath(baseOutputDir, "GDCMin_OUT")

  dccCompare=cleanFilePath(baseTestDir, "DCCMin_OUT")
  gdcCompare=cleanFilePath(baseTestDir, "GDCMin_OUT")
  ########################################################
  ########################################################
  unlink(dccMinOut, recursive=TRUE)
  dir.create(dccMinOut, showWarnings=FALSE, recursive=TRUE)
  unlink(gdcMinOut, recursive=TRUE)
  dir.create(gdcMinOut, showWarnings=FALSE, recursive=TRUE)
  ########################################################
  ########################################################
  mutationBatchAssess(theTypeCountDir=dccMinMut, theOutputDir=dccMinOut, theThreads=5,
                      thePvalueCutoff=.00001, theZScoreCutoff=1.96, thePCAflag=TRUE,
                      theBatchTypes=c("BatchId"),
                      theMutationTypes=c("FrameShiftIns", "Total"))
  mutationBatchAssess(theTypeCountDir=gdcMinMut, theOutputDir=gdcMinOut, theThreads=5,
                      thePvalueCutoff=.00001, theZScoreCutoff=1.96, thePCAflag=TRUE,
                      theBatchTypes=c("PlateId"),
                      theMutationTypes=c("FrameShiftIns", "Total"))
  ########################################################
  ########################################################
  htmlMutationBatchEffects(dccMinOut)
  htmlMutationBatchEffects(gdcMinOut)
  ########################################################
  ########################################################
  outIndex <- cleanFilePath(cleanFilePath(dccMinOut, "BatchId_FrameShiftIns"), "index.html")
  cmpIndex <- cleanFilePath(cleanFilePath(gdcCompare, "BatchId_FrameShiftIns"), "index.html")
  outCallRef <- cleanFilePath(cleanFilePath(dccMinOut, "BatchId_FrameShiftIns"), "callReference.tsv")
  cmpCallRef <- cleanFilePath(cleanFilePath(gdcCompare, "BatchId_FrameShiftIns"), "callReference.tsv")
  (md5sum(outCallRef) == md5sum(cmpCallRef))&&(md5sum(outIndex) == md5sum(cmpIndex))
} else {
  message("No test data. Skip test.")
  TRUE
}
