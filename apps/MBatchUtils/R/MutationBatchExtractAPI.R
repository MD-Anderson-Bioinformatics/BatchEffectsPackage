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

################################################################################
#### public functions
################################################################################

mutationBatchExtract <- function(theMafDir, theTypeCountDir, theGDCflag, thePar = "-Xmx8000m")
{
	# START HERE
	# theMafDir - Standardized Data timestamp directory
	#             GDC directory
	#                 path to mutations files
	#             DCC directory
	#                 path to mutations files
	# theTypeCountDir - output directory for collation mutation type data
	# theGDCflag - TRUE then handle as GDC data, FALSE handle as DCC data
	#
	# This function extracts the mutations by type and by sample and puts them into the output directory
	# the output is then processed by a separate call to mutationBatchAssess
	##############################################################################
  message(mbatchUtilVersion())
	message("mutationBatchExtract")
	message("find MAF files in ", theMafDir)
	# collect the mutation data. This file name (pattern) is for GDC-derived data and newer DCC-derived data
	fileVector <- list.files(theMafDir, pattern="mutations.tsv", full.names=TRUE, recursive=TRUE)
	message("mutations.tsv length=", length(fileVector))
	if (0==length(fileVector))
	{
	  message("MAF files")
	  # collect the mutation data. This file name (pattern) is for older DCC-derived data
		fileVector <- list.files(theMafDir, pattern=".*.maf", full.names=TRUE, recursive=TRUE)
		message("maf length=", length(fileVector))
		if (0==length(fileVector))
		{
		  message("mutation_matrix.tsv files")
		  # collect the mutation data. This file name (pattern) is for older DCC-derived data
		  fileVector <- list.files(theMafDir, pattern="mutation_matrix.tsv", full.names=TRUE, recursive=TRUE)
		  message("mutation_matrix.tsv length=", length(fileVector))
		}
	}
	# make sure the output directory exists
	dir.create(theTypeCountDir, showWarnings=FALSE, recursive=TRUE)
	##############################################################################
	# get all gene symbols and variant classifications
	# gene symbols used to populate individual files with same numbers of row
	# variable classifications used to create files, so all platforms can be compared (even if they had no mutations of that type)
	message("Collect Variant Classifications and Gene Symbols")
	variantClassifications <- c()
	geneSymbols <- c()
	for(myMaf in fileVector)
	{
		# for each file
		print(myMaf)
		# build the mutation data frame for this file
		df <- buildMutationDataframe(myMaf)
		# extract variant classifications from dataframe, and add to vector
		# TODO: Take Variant_Classification column identifier as argument
		variantClassifications <- unique(sort(unlist(c(variantClassifications, df$Variant_Classification))))
		# extract gene symbols and add to vector
		if(isTRUE(theGDCflag))
		{
			geneSymbols <- unique(sort(unlist(c(geneSymbols, df$Gene))))
		}
		else
		{
			if (is.element("Hugo_Symbol", colnames(df)))
			{
				geneSymbols <- unique(sort(unlist(c(geneSymbols, df$Hugo_Symbol))))
			}
			else
			{
				geneSymbols <- unique(sort(unlist(c(geneSymbols, df$Gene))))
			}
		}
	}
	print(variantClassifications)
	##############################################################################
	# copy the batch files to the destination for use by mutation batch assessment
	message("Copy Batch Files")
	for(myMaf in fileVector)
	{
		# for each mutation file
	  message("Copy Batch Files for each mutation file")
	  print(myMaf)
		# create the destination directory with disease sub-directory
	  message("create the destination directory with disease sub-directory")
	  dir.create(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), showWarnings=FALSE, recursive=TRUE)
		# copy batches file into disease sub-directory
	  message("copy batches file into disease sub-directory")
	  file.copy(cleanFilePath(dirname(myMaf), "batches.tsv"),
	            cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_batches.tsv", sep="")))
	}
	##############################################################################
	# pull out and collate the mutation type counts by gene for each file
	message("Collate Type Count Files")
	for(myMaf in fileVector)
	{
		# for each mutation file
		print(myMaf)
		# build the mutation data frame for this file
		df <- buildMutationDataframe(myMaf, theGDCflag)
		##########################################################################
		message("Save as Total Mutations by Gene")
		############################################################################
		# make sure disease sub-directory exists
		dir.create(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), showWarnings=FALSE, recursive=TRUE)
		# process GDC-derived and DCC-derived data differently
		# TODO: make ref build selection an argument (same as mutation type specific code above)
		if (isTRUE(theGDCflag))
		{
		  # process GDC data
		  message("Process GRCh38 Totals")
		  message("myMaf=",myMaf)
		  # name output file in disease subdirectory
		  # file name is <disease-type>.<platform>_HG38_Total.tsv
		  outputFile <- cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_HG38_Total.tsv", sep=""))
		  message("outputFile=", outputFile)
		  # extra data matrix totals from dataframe for gene symbols
		  # note call to DCC-specific matrix processing
		  # TODO: take arguments for matrix processing
		  message("calling getMatrixOfTotalsGDC")
		  dataMatrix <- getMatrixOfTotalsGDC(df, geneSymbols)
		  # write the matrix to the output file
		  message("calling writeAsMatrix")
		  writeAsGenericMatrix(outputFile, dataMatrix)
		  message("after writeAsMatrix")
		}
		else
		{
		  # process DCC data
		  refBuilds <- df$NCBI_Build
		  my37 <- (refBuilds=="37"|refBuilds=="hg19")
		  message("my37 length is ", length(my37))
		  if (sum(my37)>5)
		  {
		    message("Process GRCh37 Totals")
		    # name output file in disease subdirectory
		    # file name is <disease-type>.<platform>.<institution>.<level>_HG19_Total.tsv
		    outputFile <- cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_HG19_Total.tsv", sep=""))
		    # extract data matrix totals from dataframe for gene symbols
		    # note call to DCC-specific matrix processing
		    # TODO: take arguments for matrix processing
		    dataMatrix <- getMatrixOfTotalsDCC(df, geneSymbols, my37)
		    # write the matrix to the output file
		    writeAsGenericMatrix(outputFile, dataMatrix)
		  }
		  ##########################################################################
		  my36 <- (refBuilds=="36"|refBuilds=="hg18")
		  message("my36 length is ", length(my36))
		  if (sum(my36)>5)
		  {
		    message("Process GRCh36 Totals")
		    # name output file in disease subdirectory
		    # file name is <disease-type>.<platform>.<institution>.<level>_HG18_Total.tsv
		    outputFile <- cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_HG18_Total.tsv", sep=""))
		    # extract data matrix totals from dataframe for gene symbols
		    # note call to DCC-specific matrix processing
		    # TODO: take arguments for matrix processing
		    dataMatrix <- getMatrixOfTotalsDCC(df, geneSymbols, my36)
		    # write the matrix to the output file
		    writeAsGenericMatrix(outputFile, dataMatrix)
		  }
  	}
		############################################################################
		message("Break up by Variant Classifications by Gene")
		############################################################################
		# break up single mutation data file into one file per variant classification
		# determine which variant classificiation to process
		# GDC-derived data has only HG38, DCC-derived data may have HG18 and/or HG19
		refBuilds <- "HG38"
		if (!isTRUE(theGDCflag))
		{
			refBuilds <- df$NCBI_Build
		}
		# iterate through variant classifications, writing a single file for each mutation type
		for (callType in variantClassifications)
		{
			# make sure disease name sub-directory exists
			dir.create(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), showWarnings=FALSE, recursive=TRUE)
			# extract data and write file
			if (isTRUE(theGDCflag))
			{
				# for GDC data, for this mutation type
				message(callType)
				message("Process GRCh38")
				# name output file in disease subdirectory
				# file name is <disease-type>.<platform>_HG38_<mutation-type>.tsv
				outputFile <- cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_HG38_",  cleanString(callType), ".tsv", sep=""))
				# extract data matrix from dataframe for gene symbols and mutation type
				# note call to GDC-specific matrix processing
				# TODO: take arguments for matrix processing
				message("browser 1")
				#browser()
				dataMatrix <- getMatrixOfTypeGDC(df, geneSymbols, callType)
				# write the matrix to the output file
				writeAsGenericMatrix(outputFile, dataMatrix)
			}
			else
			{
				# for DCC data, for this mutation type, look for HG19 and HG18 data (may have either or both)
				##########################################################################
				# some DCC files use 37 and some use hg19 for reference build id
				# TODO: make ref build selection an argument
				my37 <- (refBuilds=="37"|refBuilds=="hg19")
				message("my37 length is ", length(my37))
				if (sum(my37)>5)
				{
					message("Process GRCh37")
					# name output file in disease subdirectory
					# file name is <disease-type>.<platform>.<institution>.<level>_HG19_<mutation-type>.tsv
					outputFile <- cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_HG19_",  cleanString(callType), ".tsv", sep=""))
					# extract data matrix from dataframe for gene symbols and mutation type
					# note call to DCC-specific matrix processing
					# TODO: take arguments for matrix processing
					message("browser 2")
					#browser()
					dataMatrix <- getMatrixOfTypeDCC(df, geneSymbols, my37, callType)
					# write the matrix to the output file
					writeAsGenericMatrix(outputFile, dataMatrix)
				}
				##########################################################################
				my36 <- (refBuilds=="36"|refBuilds=="hg18")
				message("my36 length is ", length(my36))
				if (sum(my36)>5)
				{
					message("Process GRCh36")
					# name output file in disease subdirectory
					# file name is <disease-type>.<platform>.<institution>.<level>_HG18_<mutation-type>.tsv
					outputFile <- cleanFilePath(cleanFilePath(theTypeCountDir, getProjectDiseaseIdentifier(myMaf, theGDCflag)), paste(convertPathToName(myMaf, theGDCflag), "_HG18_",  cleanString(callType), ".tsv", sep=""))
					# extract data matrix from dataframe for gene symbols and mutation type
					# note call to DCC-specific matrix processing
					# TODO: take arguments for matrix processing
					message("browser 3")
					#browser()
					dataMatrix <- getMatrixOfTypeDCC(df, geneSymbols, my36, callType)
					# write the matrix to the output file
					writeAsGenericMatrix(outputFile, dataMatrix)
				}
				##########################################################################
			}
		}
	}
}

################################################################################
#### internal functions
################################################################################

getMatrixOfTotalsGDC <- function(theDF, theGenes)
{
	# create a matrix of samples (columns) by genes (rows) with mutation totals
	# theDF - dataframe loaded from mutation data file (mitochondrial and other non 1-22-X-Y data remove, chromosomes may be numbers or chr#)
	# theGenes - vector of gene symbols used for building matrix
	#
	# Note, this uses dataframe column ids specific to GDC-derived MAF/mutation data files
	# TODO: combine/simplify getMatrix type functions into single function
	#
	#Get sorted list of sample ids (barcodes)
	samples <- unique(sort(theDF$Tumor_Sample_Barcode))
	# build matrix of zeros for samples (columns) and genes (rows)
	myM <- matrix(data=0, nrow=length(theGenes), ncol=length(samples), dimnames=list(theGenes, samples))
	# populate the matrix by counting how many times each barcode is in MAF file
	# MAF files have one line per mutation call, which is why this works
	for (index in 1:length(theDF$Tumor_Sample_Barcode))
	{
		# pull out gene called
	  hugoSymbol <- NULL
	  # pull out gene called
	  if (is.element("Hugo_Symbol", colnames(theDF)))
	  {
	    hugoSymbol <- theDF$Hugo_Symbol[index]
	  }
	  else
	  {
	    hugoSymbol <- theDF$Gene[index]
	  }
		# pull out sample (barcode) tested
		sampleId <- theDF$Tumor_Sample_Barcode[index]
		# increment value in matrix
		if (isTRUE(sampleId %in% samples ))
		{
		  if (isTRUE(hugoSymbol %in% theGenes ))
		  {
		    myM[hugoSymbol, sampleId] <- myM[hugoSymbol, sampleId] + 1
		  }
		  else
		  {
		    message("getMatrixOfTotalsGDC: hugoSymbol not found in matrix = ", hugoSymbol)
		  }
		}
		else
		{
		  message("getMatrixOfTotalsGDC: sampleId not found in matrix = ", sampleId)
		}
	}
	myM
}

getMatrixOfTotalsDCC <- function(theDF, theGenes, theUseSelection)
{
	# create a matrix of samples (columns) by genes (rows) with mutation totals
	# theDF - dataframe loaded from mutation data file (mitochondrial and other non 1-22-X-Y data remove, chromosomes may be numbers or chr#)
	# theGenes - vector of gene symbols used for building matrix
	# theUseSelection - a vector of TRUE/FALSE indicating whether to count each row, used to separate HG18 and HG19 calls
	#
	# Note, this uses dataframe column ids specific to DCC-derived MAF/mutation data files
	# TODO: combine/simplify getMatrix type functions into single function
	#
	#Get sorted list of sample ids (barcodes)
	samples <- unique(sort(theDF$Tumor_Sample_Barcode[theUseSelection]))
	# build matrix of zeros for samples (columns) and genes (rows)
	myM <- matrix(data=0, nrow=length(theGenes), ncol=length(samples), dimnames=list(theGenes, samples))
	# populate the matrix by counting how many times each barcode is in MAF file
	# MAF files have one line per mutation call, which is why this works
	for (index in 1:length(theUseSelection))
	{
		# check if this row is in the desired reference build
		if (isTRUE(theUseSelection[index]))
		{
			hugoSymbol <- NULL
			# pull out gene called
			if (is.element("Hugo_Symbol", colnames(theDF)))
			{
				hugoSymbol <- theDF$Hugo_Symbol[index]
			}
			else
			{
				hugoSymbol <- theDF$Gene[index]
			}
			# pull out sample (barcode) tested
			sampleId <- theDF$Tumor_Sample_Barcode[index]
			# increment value in matrix
			myM[hugoSymbol, sampleId] <- myM[hugoSymbol, sampleId] + 1
		}
	}
	myM
}

getMatrixOfTypeGDC <- function(theDF, theGenes, theClassificationType)
{
	# create a matrix of samples (columns) by genes (rows) with mutation totals
	# theDF - dataframe loaded from mutation data file (mitochondrial and other non 1-22-X-Y data remove, chromosomes may be numbers or chr#)
	# theGenes - vector of gene symbols used for building matrix
	# theClassificationType - a string giving the type of mutation to count for this dataframe
	#
	# Note, this uses dataframe column ids specific to GDC-derived MAF/mutation data files
	# TODO: combine/simplify getMatrix type functions into single function
	#
	#Get sorted list of sample ids (barcodes)
	samples <- unique(sort(theDF$Tumor_Sample_Barcode))
	# build matrix of zeros for samples (columns) and genes (rows)
	myM <- matrix(data=0, nrow=length(theGenes), ncol=length(samples), dimnames=list(theGenes, samples))
	# populate the matrix by counting how many times each barcode is in MAF file with the passed-in mutation type
	# MAF files have one line per mutation call, which is why this works
	for (index in 1:length(theDF$Tumor_Sample_Barcode))
	{
	  hugoSymbol <- NULL
	  # pull out gene called
	  if (is.element("Hugo_Symbol", colnames(theDF)))
	  {
	    hugoSymbol <- theDF$Hugo_Symbol[index]
	  }
	  else
	  {
	    hugoSymbol <- theDF$Gene[index]
	  }
		# pull out sample (barcode) tested
		sampleId <- theDF$Tumor_Sample_Barcode[index]
		# pull out variant type called
		classType <- theDF$Variant_Classification[index]
		# check if variant called is the type we are counting
		if (theClassificationType==classType)
		{
			# increment value in matrix, if of correct type
			myM[hugoSymbol, sampleId] <- myM[hugoSymbol, sampleId] + 1
		}
	}
	myM
}

getMatrixOfTypeDCC <- function(theDF, theGenes, theUseSelection, theClassificationType)
{
	# create a matrix of samples (columns) by genes (rows) with mutation totals
	# theDF - dataframe loaded from mutation data file (mitochondrial and other non 1-22-X-Y data remove, chromosomes may be numbers or chr#)
	# theGenes - vector of gene symbols used for building matrix
	# theUseSelection - a vector of TRUE/FALSE indicating whether to count each row, used to separate HG18 and HG19 calls
	# theClassificationType - a string giving the type of mutation to count for this dataframe
	#
	# Note, this uses dataframe column ids specific to DCC-derived MAF/mutation data files
	# TODO: combine/simplify getMatrix type functions into single function
	#
	#Get sorted list of sample ids (barcodes)
	samples <- unique(sort(theDF$Tumor_Sample_Barcode[theUseSelection]))
	# build matrix of zeros for samples (columns) and genes (rows)
	myM <- matrix(data=0, nrow=length(theGenes), ncol=length(samples), dimnames=list(theGenes, samples))
	# populate the matrix by counting how many times each barcode is in MAF file with the passed-in mutation type
	# MAF files have one line per mutation call, which is why this works
	for (index in 1:length(theUseSelection))
	{
		# check if this row is in the desired reference build
		if (isTRUE(theUseSelection[index]))
		{
			hugoSymbol <- NULL
			# pull out gene called
			if (is.element("Hugo_Symbol", colnames(theDF)))
			{
				hugoSymbol <- theDF$Hugo_Symbol[index]
			}
			else
			{
				hugoSymbol <- theDF$Gene[index]
			}
			# pull out sample (barcode) tested
			sampleId <- theDF$Tumor_Sample_Barcode[index]
			# pull out variant type called
			classType <- theDF$Variant_Classification[index]
			# check if variant called is the type we are counting
			if (theClassificationType==classType)
			{
				# increment value in matrix, if of correct type
				myM[hugoSymbol, sampleId] <- myM[hugoSymbol, sampleId] + 1
			}
		}
	}
	myM
}

cleanString <- function(theString)
{
	# used to reduce platform names and mutation type names into something safe for a file name
	# remove everything but alphanumeric
	gsub("[^[:alnum:]=\\.]", "", theString)
}

convertPathToName <- function(thePath, theGDCflag)
{
  message("convertPathToName")
	# Pull the "unique" data out of path to form file names for GDC and DCC data
	# thePath - path to data from which to extract name
	# theGDCflag - TRUE means use GDC format, FALSE means use DCC format
	#
	# TODO: use arguments to create or pass in
	# for GDC <disease-type>.<platform>.<institution>.<level>_batches.tsv
	# for DCC <disease-type>.<platform>_batches.tsv and
	name <- NULL
	if (isTRUE(theGDCflag))
	{
		# for GDC <disease-type>.<platform>.<institution>.<level>_batches.tsv
		name <- paste(
		        getProjectDiseaseIdentifier(thePath, theGDCflag), # TCGA-ACC
		        basename(dirname(dirname(dirname(dirname(thePath))))), # MuSEVariant...
						sep=".")
	}
	else
	{
		# for DCC <disease-type>.<platform>_batches.tsv and
		name <- paste(basename(dirname(dirname(dirname(dirname(thePath))))), # kirp
									paste(head(strsplit(x=basename(dirname(dirname(thePath))), split="_", fixed=TRUE)[[1]], n=1), sep="_"), # mixed_dnaseq_curated
									gsub(x=tail(strsplit(x=basename(dirname(dirname(thePath))), split="_", fixed=TRUE)[[1]], n=1), pattern=".", replacement="_", fixed=TRUE), # hgsc_bcm_edu
									basename(dirname(thePath)), # Level_2
									sep=".")
	}
	return(name)
}

getProjectDiseaseIdentifier <- function(theFile, theGDCflag)
{
  message("getProjectDiseaseIdentifier")
	# Pull the disease type out of path to form sub-directory
	# theFile - file and path to file. Use path to get disease type
	# theGDCflag - TRUE means use GDC format, FALSE means use DCC format
	#
	# TODO: use arguments to create or pass in
	pdId <- NULL
	if (isTRUE(theGDCflag))
	{
		# GDC disease types
		pdId <- basename(dirname(dirname(dirname(dirname(dirname(dirname(theFile)))))))
	}
	else
	{
		# DCC disease types
		pdId <- basename(dirname(dirname(dirname(dirname(theFile)))))
	}
	return(pdId)
}

isChromsomeP <- function(theList)
{
	# theList - the lsit of chromosome values, without the "chr" that some files have
	# Return: return TRUE (is a 1-22, X, Y Chromosome or FALSE is not
	# this is used to filter out mitochondrial DNA
	as.vector(unlist(lapply(theList, function(theChromosome)
	{
		result <- FALSE
		if (("1"==theChromosome)||("2"==theChromosome)||("3"==theChromosome)||("4"==theChromosome)||
				("5"==theChromosome)||("6"==theChromosome)||("7"==theChromosome)||("8"==theChromosome)||
				("9"==theChromosome)||("10"==theChromosome)||("11"==theChromosome)||("12"==theChromosome)||
				("13"==theChromosome)||("14"==theChromosome)||("15"==theChromosome)||("16"==theChromosome)||
				("17"==theChromosome)||("18"==theChromosome)||("19"==theChromosome)||("20"==theChromosome)||
				("21"==theChromosome)||("22"==theChromosome)||("Y"==theChromosome)||("X"==theChromosome))
		{
			result <- TRUE
		}
		result
	})))
}

readDataframeScore <- function(file)
{
	# read everything as strings, don't modify column names, and don't convert to factors
	read.delim(file, sep="\t", as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE, colClasses=c("character"))
}

buildMutationDataframe <- function(theMafFile, theGDCflag=FALSE)
{
	# theMafFile - the full path and file name for MAF file to read
	# theGDCflag - TRUE means use GDC processing, FALSE means DCC
	# used internally to build a mutation dataframe from the MAF file data
	message("read MAF ", theMafFile)
	# read as all-character dataframe
	myDF <- readDataframeScore(theMafFile)
	# remove mitochondrial and other DNA
	message("remove mitochondrial and other non 1-22-X-Y data")
	keep <- NULL
	if (isTRUE(theGDCflag))
	{
		# remove "chr" from chromosome name vector, and pass to isChromosome
		keep <- isChromsomeP(gsub(myDF$Position_Chromosome, pattern="chr", replacement="", fixed=TRUE))
	}
	else
	{
		# pass to isChromosome
		if (is.element("Chromosome", colnames(df)))
		{
			keep <- isChromsomeP(myDF$Chromosome)
		}
		else
		{
			keep <- isChromsomeP(myDF$Position_Chromosome)
		}
	}
	# keep only 1-21, X, and Y chromosome data
	if (sum(keep)>0)
	{
		myDF <- myDF[keep,]
	}
	myDF
}
