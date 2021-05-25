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

require(MBatchUtils)

if (!is.null(getTestOutputDir()))
{
  #warnLevel<-getOption("warn")
  #on.exit(options(warn=warnLevel))
  # warnings are errors
  #options(warn=3)
  # if there is a warning, show the calls leading up to it
  #options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  #options(showErrorCalls=TRUE)
  ########################################################
  ########################################################
  baseTestDir=getTestInputDir()
  dccMinMaf=file.path(baseTestDir, "DCCMin_MAF")
  gdcMinMaf=file.path(baseTestDir, "GDCMin_MAF")
  dccMinMutCompare=file.path(baseTestDir, "DCCMin_MUT")
  gdcMinMutCompare=file.path(baseTestDir, "GDCMin_MUT")
  ########################################################
  baseOutputDir=getTestOutputDir()
  dccMinMut=file.path(baseOutputDir, "DCCMin_MUT")
  gdcMinMut=file.path(baseOutputDir, "GDCMin_MUT")
  ########################################################
  ########################################################
  unlink(dccMinMut, recursive=TRUE)
  dir.create(dccMinMut, showWarnings=FALSE, recursive=TRUE)
  unlink(gdcMinMut, recursive=TRUE)
  dir.create(gdcMinMut, showWarnings=FALSE, recursive=TRUE)
  ########################################################
  ########################################################
  print(paste("mutationBatchExtract(theMafDir=", dccMinMaf, ", theGDCflag=FALSE, theTypeCountDir=",dccMinMut))
  mutationBatchExtract(theMafDir=dccMinMaf, theGDCflag=FALSE, theTypeCountDir=dccMinMut)
  print(paste("mutationBatchExtract(theMafDir=", dccMinMaf, ", theGDCflag=TRUE, theTypeCountDir=",dccMinMut))
  mutationBatchExtract(theMafDir=gdcMinMaf, theGDCflag=TRUE, theTypeCountDir=gdcMinMut)
  ########################################################
  ########################################################
  print(paste("readAsGenericMatrix(", file.path(dccMinMut, "acc", "acc.illuminaga.hgsc_bcm_edu.Level_2_HG19_Total.tsv"),")"))
  correctedMatrix <- readAsGenericMatrix(file.path(dccMinMut, "acc", "acc.illuminaga.hgsc_bcm_edu.Level_2_HG19_Total.tsv"))
  print(paste("readAsGenericMatrix(", file.path(dccMinMutCompare, "acc", "acc.illuminaga.hgsc_bcm_edu.Level_2_HG19_Total.tsv"),")"))
  compareMatrix <- readAsGenericMatrix(file.path(dccMinMutCompare, "acc", "acc.illuminaga.hgsc_bcm_edu.Level_2_HG19_Total.tsv"))
  print("compared1")
  compared1 <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared1)
  print(file.path(gdcMinMut, "GDCMin_MAF", "GDCMin_MAF.TCGA_HG38_Total.tsv"))
  print(file.path(gdcMinMutCompare, "TCGA-ACC", "TCGA-ACC.VarScan2VariantAggregationandMasking_HG38_Total.tsv"))
  print("correctedMatrix 1")
  correctedMatrix <- readAsGenericMatrix(file.path(gdcMinMut, "GDCMin_MAF", "GDCMin_MAF.TCGA_HG38_Total.tsv"))
  print("compareMatrix 1")
  compareMatrix <- readAsGenericMatrix(file.path(gdcMinMutCompare, "TCGA-ACC", "TCGA-ACC.VarScan2VariantAggregationandMasking_HG38_Total.tsv"))
  print("compared2")
  compared2 <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared2)
  compared1&&compared2
} else {
  message("No test data. Skip test.")
  TRUE
}
