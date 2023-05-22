//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <math.h>
//#include <pthread.h>
//#include <errno.h>
//#include <sys/types.h>
//#include <time.h>
//#include <sys/time.h>

//#include <unistd.h>
/* #include <sys/syscall.h> */
//#include <sched.h>
/* TDC replaced with below #include <alloca.h> */
//#include <malloc.h>

//#define CSTACK_DEFNS 1;

/*--*/
#include <R.h>
#include <Rdefines.h>
//#include <Rinterface.h>
#include <Rmath.h>

#include <pthread.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>

#include "compute.h"
#include "perm.h"
#include "c_threads.h"

// called from R
void getDSCwithExcerpt(double * theValues, int * theDim, char ** theBatchIds, 
		char ** theSortedUniqueBatchIds, int * theNumberOfUniqueBatchIds,
		double * retListOfGeneDSC, double * retListOfGeneDB, double * retListOfGeneDW, 
		double * retDSC, double * retDB, double * retDW)
{
	//R_CStackLimit = (uintptr_t)-1;
	gsl_matrix * myMatrix = gsl_matrix_alloc(theDim[0], theDim[1]);
	copyMatrixFromArray(myMatrix, theValues, theDim);
	gsl_vector * rowVector = gsl_vector_calloc(theDim[1]);
	gsl_vector * batchVectorList[*theNumberOfUniqueBatchIds];
	gsl_vector * valuesForBatch[*theNumberOfUniqueBatchIds];
	int x = 0;
	for(x=0;x<*theNumberOfUniqueBatchIds;x++)
	{
		char * batchName = theSortedUniqueBatchIds[x];
		int numSamplesInBatch = countOccurances(theBatchIds, theDim[1], batchName);
		gsl_vector * temp = gsl_vector_calloc(numSamplesInBatch);
		batchVectorList[x] = temp;
		temp = gsl_vector_calloc(numSamplesInBatch);
		valuesForBatch[x] = temp;
	}
	compute(myMatrix, rowVector, batchVectorList, valuesForBatch,
			theDim[0], theBatchIds, theDim[1], theSortedUniqueBatchIds, *theNumberOfUniqueBatchIds,
			retListOfGeneDSC, retListOfGeneDB, retListOfGeneDW, retDSC, retDB, retDW);
	//Rprintf("getDSCwithExcerpt Results DSC =%.30f\n" , *retDSC);
	//Rprintf("getDSCwithExcerpt Results DB =%.30f\n" , *retDB);
	//Rprintf("getDSCwithExcerpt Results DW =%.30f\n" , *retDW);
	for(x=0;x<*theNumberOfUniqueBatchIds;x++)
	{
		gsl_vector_free(batchVectorList[x]);
		gsl_vector_free(valuesForBatch[x]);
	}
	gsl_vector_free(rowVector);
	gsl_matrix_free(myMatrix);
}

//PcaDsc 
//	public PcaDsc [] doDscPerms(double [] theValues, int [] theDim, String [] theBatchIds, int thePerms, int theThreads )
void doDscPerms(double * theValues, int * theDim, char ** theBatchIds, 
		char ** theSortedUniqueBatchIds, int * theNumberOfUniqueBatchIds,
		int * thePerms, int * theThreads,
		double * retListOfGeneDSC, double * retListOfGeneDB, double * retListOfGeneDW, 
		double * retDSC, double * retDB, double * retDW)
{
	//R_CStackLimit = (uintptr_t)-1;
	//Rprintf("thePerms %d\n", *thePerms);
	//Rprintf("theThreads %d\n", *theThreads);
	//ForkJoinPool pool = new ForkJoinPool(theThreads);
	//ArrayList<PcaDsc> myresults = pool.invoke(new Perm(matrix, theBatchIds, thePerms));
	//System.out.println("results contains " + myresults.size());
	//results = myresults.toArray(new PcaDsc[0]);
	// setup for threads
	NUM_PERMUTATIONS = *thePerms;
	NEXT_PERMUTATION = 0;
	ARG_Values = theValues;
	ARG_Dim = theDim;
	ARG_BatchIds = theBatchIds;
	ARG_SortedUniqueBatchIds = theSortedUniqueBatchIds;
	ARG_NumberOfUniqueBatchIds = theNumberOfUniqueBatchIds;
	ARG_retListOfGeneDSC = retListOfGeneDSC;
	ARG_retListOfGeneDB = retListOfGeneDB;
	ARG_retListOfGeneDW = retListOfGeneDW;
	ARG_retDSC = retDSC;
	ARG_retDB = retDB;
	ARG_retDW = retDW;
	pthread_t myThreads[*theThreads];
	int x = 0;
	for(x=0;x<*theThreads;x++)
	{
		pthread_create(&myThreads[x], NULL, threadPermCompute, NULL);
	}
	for(x=0;x<*theThreads;x++)
	{
		pthread_join(myThreads[x], NULL);
	}
}

void dsccTest(double * theDouble, int * theInt, char ** theString,
		double * theDoubleList, int * theDLLength, 
		int * theIntList, int * theILLength, 
		char ** theStringList, int * theSLLength,
		double * theMatrixValues, int * theMatrixDim)
{
	Rprintf("5\n");
	Rprintf("theDouble %f\n", *theDouble);
	Rprintf("theInt %d\n", *theInt);
	Rprintf("theString %s\n", *theString);
	int x = 0;
	int y = 0;
	Rprintf("double list ");
	for(x=0;x<*theDLLength;x++)
	{
		Rprintf("%f, ", theDoubleList[x]);
	}
	Rprintf("\n");
	Rprintf("int list ");
	for(x=0;x<*theILLength;x++)
	{
		Rprintf("%d, ", theIntList[x]);
	}
	Rprintf("\n");
	Rprintf("string list ");
	for(x=0;x<*theSLLength;x++)
	{
		Rprintf("%s, ", theStringList[x]);
	}
	Rprintf("\n");
	Rprintf("theMatrixValues list ");
	for(x=0;x<(theMatrixDim[0]*theMatrixDim[1]);x++)
	{
		Rprintf("%f, ", theMatrixValues[x]);
	}
	Rprintf("\n");
	Rprintf("matrix list\n");
	Rprintf("theMatrixDim[0] %d\n", theMatrixDim[0]);
	Rprintf("theMatrixDim[1] %d\n", theMatrixDim[1]);
	// create a matrix (two dimensional array) from the array of theValues
	gsl_matrix * myMatrix = gsl_matrix_alloc(theMatrixDim[0], theMatrixDim[1]);
	int index = 0;
	for(y=0;y<theMatrixDim[1];y++)
	{
		for(x =0;x<theMatrixDim[0];x++)
		{
			gsl_matrix_set(myMatrix, x, y, theMatrixValues[index]);
			index = index + 1;
		};
	};
	for(x=0;x<theMatrixDim[0];x++)
	{
		for(y=0;y<theMatrixDim[1];y++)
		{
			Rprintf("%f, ", gsl_matrix_get(myMatrix, x, y));
		}
		Rprintf("\n");
	}
	gsl_matrix_free(myMatrix);
}
