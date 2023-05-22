#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <pthread.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>

#include "compute.h"
#include "c_threads.h"

void copyMatrixFromArray(gsl_matrix * theMatrix, double * theMatrixValues, int * theMatrixDim)
{
	//Rprintf("createMatrix\n");
	// create a matrix (two dimensional array) from the array of theValues
	//gsl_matrix * myMatrix = gsl_matrix_alloc(theMatrixDim[0], theMatrixDim[1]);
	int index = 0;
	int x = 0;
	int y = 0;
	for(y=0;y<theMatrixDim[1];y++)
	{
		for(x =0;x<theMatrixDim[0];x++)
		{
			gsl_matrix_set(theMatrix, x, y, theMatrixValues[index]);
			index = index + 1;
		};
	};
	//return(myMatrix);
}

void doPermutation(gsl_matrix * theMatrix, gsl_vector * theRowVector,
		gsl_vector ** theBatchVectorList, gsl_vector ** theValuesForBatch, 
		double * theValues, int * theDim, char ** theBatchIds, 
		char ** theSortedUniqueBatchIds, int * theNumberOfUniqueBatchIds,
		double * retListOfGeneDSC, double * retListOfGeneDB, double * retListOfGeneDW, 
		double * retDSC, double * retDB, double * retDW, gsl_rng * theGlsRng)
{
	//Rprintf("doPermutation\n");
	copyMatrixFromArray(theMatrix, theValues, theDim);
	//Rprintf("vSize %d\n", vSize);
	// permutate the values in each row
	//Rprintf("range size %d\n", gsl_rng_size(r));
	int x = 0;
	//Rprintf("rowVector calloc\n");
	for(x=0;x<theDim[0];x++)
	{
		gsl_matrix_get_row(theRowVector, theMatrix, x);
		gsl_ran_shuffle(theGlsRng, theRowVector->data, theRowVector->size, sizeof(double));
		gsl_matrix_set_row(theMatrix, x, theRowVector);
	}
	//Rprintf("after loop\n");
	compute(theMatrix, theRowVector, theBatchVectorList, theValuesForBatch,
			theDim[0], theBatchIds, theDim[1], theSortedUniqueBatchIds, *theNumberOfUniqueBatchIds,
			retListOfGeneDSC, retListOfGeneDB, retListOfGeneDW, retDSC, retDB, retDW);
	//Rprintf("Perm Results DSC =%.30f\n" , *retDSC);
	//Rprintf("Perm Results DB =%.30f\n" , *retDB);
	//Rprintf("Perm Results DW =%.30f\n" , *retDW);
}
