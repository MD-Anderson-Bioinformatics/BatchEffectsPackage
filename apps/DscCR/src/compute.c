#include <string.h>
#include <math.h>

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>

int countOccurances(char ** theArray, int theArraySize, char * theString)
{
	int count = 0;
	int x = 0;
	for (x=0;x<theArraySize;x++)
	{
		if (0==strcmp(theArray[x], theString))
		{
			count = count + 1;
		}
	}
	return count;
}

void getBatchOccurances(gsl_vector * theBatchOccurances, gsl_vector * theRowVector, 
		char ** theBatchIds, int theNumberOfBatchIds, char * theBatchName)
{
	int count = 0;
	int x = 0;
	for (x = 0; x<theNumberOfBatchIds; x++)
	{
		if (0==strcmp(theBatchIds[x], theBatchName))
		{
			gsl_vector_set(theBatchOccurances, count, gsl_vector_get(theRowVector, x));
			count = count + 1;
		}
	}
}

void removeNaN(gsl_vector * theBatchOccurances, gsl_vector * theBatchVector, int theNumSamplesInBatch)
{
	int x = 0;
	int count = 0;
	for (x = 0; x<theNumSamplesInBatch; x++)
	{
		double val = gsl_vector_get(theBatchVector, x);
		if (!isnan(val))
		{
			gsl_vector_set(theBatchOccurances, count, val);
			count = count + 1;
		}
	}
}

void displayVector(gsl_vector * theVector, int theLimit)
{
	int x = 0;
	int count = 0;
	Rprintf("length = %d, ", theVector->size);
	if (0==theLimit)
	{
		theLimit = theVector->size;
	}
	for (x = 0; x<theLimit; x++)
	{
		Rprintf("%.60f, ", gsl_vector_get(theVector, x));
	}
	Rprintf("\n");
}

double epsilonZeroCheck(double theValue)
{
	// if theValue is less than or equal 1x10^-7, return 0, otherwise return theValue
	if (isnan(theValue))
	{
		return(theValue);
	}
	else if (theValue<=0.0000001)
	{
		return(0);
	}
	else
	{
		return(theValue);
	}
}

void copyArrayToWithEpsilonCheck(double * theToArray, double * theFromArray, int theLength)
{
	int x = 0;
	for(x=0;x<theLength;x++)
	{
		theToArray[x] = epsilonZeroCheck(theFromArray[x]);
	}
}

void copyArrayTo(double * theToArray, double * theFromArray, int theLength)
{
	int x = 0;
	for(x=0;x<theLength;x++)
	{
		theToArray[x] = theFromArray[x];
	}
}

void compute(gsl_matrix * theMatrix, gsl_vector * theRowVector,
		gsl_vector ** theBatchVectorList, gsl_vector ** theValuesForBatch, 
		int theNumberOfSamples, 
		char ** theBatchIds, int theNumberOfBatchIds, 
		char ** theSortedUniqueBatchIds, int theNumberOfUniqueBatchIds,
		double * retListOfGeneDSC, double * retListOfGeneDB, double * retListOfGeneDW, 
		double * retDSC, double * retDB, double * retDW)
{
	//Rprintf("compute\n");
//JAVA// 	PcaDsc myObj = new PcaDsc();
	// number of samples
	//Rprintf("size1 %d\n", theMatrix->size1);
	//Rprintf("size2 %d\n", theMatrix->size2);
	//Rprintf("tda %d\n", theMatrix->tda);
	//Rprintf("theNumberOfSamples %d\n", theNumberOfSamples);
	//Rprintf("theNumberOfBatchIds %d\n", theNumberOfBatchIds);
	//Rprintf("theNumberOfUniqueBatchIds %d\n", theNumberOfUniqueBatchIds);
	int numberOfSamples = theNumberOfSamples;
	//Rprintf("mean of across all the samples in the dataset\n");
	double meanListAcrossSamples[numberOfSamples];
	//Rprintf("samples %d\n", numberOfSamples);
	//displayVector(rowVector);
	//Rprintf("Mean across samples\n");
	int x = 0;
	int y = 0;
	for(x=0;x<numberOfSamples;x++)
	{
		//displayMatrixRow(theMatrix, x, numberOfSamples);
		//Rprintf("copy row to vector %d\n", x);
		gsl_matrix_get_row(theRowVector, theMatrix, x);
		//Rprintf("calculate mean\n");
		meanListAcrossSamples[x] = gsl_stats_mean(theRowVector->data,theRowVector->stride,theRowVector->size);
		//Rprintf("%.60f, ", meanListAcrossSamples[x]);
	}
	//Rprintf("\n");
	double dwValue = 0.0;
	double dbValue = 0.0;
	// this removes duplicates in batch ids (it also sorts)
	char ** batchIdList = theSortedUniqueBatchIds;
	//
	double gDSC[theNumberOfBatchIds];
	double gDW[theNumberOfBatchIds];
	double gDB[theNumberOfBatchIds];
	// for each gene
	//Rprintf("for each gene %d\n", theNumberOfBatchIds);
	for(x=0;x<theNumberOfBatchIds;x++)
	{
		//Rprintf("loop 1\n");
		double geneDw = 0.0;
		double geneDb = 0.0;
		//Rprintf("theNumberOfUniqueBatchIds %d\n", theNumberOfUniqueBatchIds);
		for (y=0;y<theNumberOfUniqueBatchIds;y++)
		{
			//Rprintf("loop 2\n");
			char * batchName = batchIdList[y];
			//Rprintf("loop 3\n");
			gsl_matrix_get_row(theRowVector, theMatrix, x);
			//Rprintf("loop 4\n");
			int numSamplesInBatch = countOccurances(theBatchIds, theNumberOfBatchIds, batchName);
			//Rprintf("loop 5\n");
			getBatchOccurances(theBatchVectorList[y], theRowVector, theBatchIds, theNumberOfBatchIds, batchName);
			//Rprintf("batchVector\n");
			//displayVector(batchVector);
			//Rprintf("loop 6\n");
			removeNaN(theValuesForBatch[y], theBatchVectorList[y], numSamplesInBatch);
			//Rprintf("valuesForBatch\n");
			//displayVector(valuesForBatch);
			//Rprintf("loop 7\n");
			double batchMean = gsl_stats_mean(theValuesForBatch[y]->data,theValuesForBatch[y]->stride,theValuesForBatch[y]->size);
			//Rprintf("loop 8\n");
			double sampleVariance = gsl_stats_variance(theValuesForBatch[y]->data,theValuesForBatch[y]->stride,theValuesForBatch[y]->size);
			//Rprintf("loop 9\n");
			geneDw = geneDw + ( ((double)(numSamplesInBatch-1)) / ((double)(numberOfSamples)) * sampleVariance );
			//Rprintf("loop 10\n");
			geneDb = geneDb + ( (((double)(numSamplesInBatch)) / (double)(numberOfSamples)) * 
					( batchMean - meanListAcrossSamples[x] ) * 
					( batchMean - meanListAcrossSamples[x] ) );
			//Rprintf("b=" + batchName + " nsb=" + numSamplesInBatch + " ns=" + numberOfSamples + " bm=" + batchMean + " sv=" + sampleVariance + " geneDw=" + geneDw + " geneDb=" + geneDb + " m=" + meanListAcrossSamples[x]);
			//Rprintf("loop 11\n");
			//Rprintf("b=%s nsb=%d ns=%d bm=%.60f sv=%.60f geneDw=%.60f geneDb=%.60f m=%.60f\n", 
			//		batchName, numSamplesInBatch, numberOfSamples, batchMean, sampleVariance, 
			//		geneDw, geneDb, meanListAcrossSamples[x]);
		}
		//Rprintf("after y loop\n");
		dwValue = dwValue + geneDw;
		dbValue = dbValue + geneDb;
		// take square root
		geneDw = sqrt(geneDw);
		geneDb = sqrt(geneDb);
		//
		double geneDSC = NAN;
		//if ((geneDw!=BigDecimal.NaN)&&(1==geneDw.compareTo(BigDecimal.ZERO)))
		if ((!isnan(geneDw)) &&(geneDw>0))
		{
			geneDSC = geneDb/geneDw;
		}
		gDSC[x] = geneDSC;
		gDW[x] = geneDw;
		gDB[x] = geneDb;
	}			
	//Rprintf("after x loop\n");
	// take square root
	dwValue = sqrt(dwValue);
	dbValue = sqrt(dbValue);
	// avoid divide by zero error
	// message("avoid divide by zero error")
	// message("dwValue ", dwValue)
	double DSC = NAN;
	if ((!isnan(dwValue))&&(dwValue > 0))
	{
		DSC = dbValue/dwValue;
	}
	if (theNumberOfBatchIds<2)
	{
		//Rprintf("set arrays to 0.0\n");
		for(x=0;x<theNumberOfBatchIds;x++)
		{
			gDSC[x] = 0.0; //nrow(mMatrix)
			gDB[x] = 0.0; //nrow(mMatrix))
		}
		DSC = 0.0;
		dbValue = 0.0;
	}
	copyArrayToWithEpsilonCheck(retListOfGeneDSC, gDSC, theNumberOfBatchIds);
	copyArrayTo(retListOfGeneDB, gDW, theNumberOfBatchIds);
	copyArrayTo(retListOfGeneDW, gDB, theNumberOfBatchIds);
	*retDSC = epsilonZeroCheck(DSC);
	*retDB = dbValue;
	*retDW = dwValue;
	//Rprintf("RESULTS results.mDSC =%.30f\n" , *retDSC);
	//Rprintf("RESULTS results.mDB =%.30f\n" , *retDB);
	//Rprintf("RESULTS results.mDW =%.30f\n" , *retDW);
}

