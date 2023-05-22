#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <pthread.h>

#include "perm.h"

int NUM_PERMUTATIONS;
int NEXT_PERMUTATION;

double * ARG_Values;
int * ARG_Dim;
char ** ARG_BatchIds;
char ** ARG_SortedUniqueBatchIds;
int * ARG_NumberOfUniqueBatchIds;
double * ARG_retListOfGeneDSC;
double * ARG_retListOfGeneDB;
double * ARG_retListOfGeneDW;
double * ARG_retDSC;
double * ARG_retDB;
double * ARG_retDW;

pthread_mutex_t LOCK_FOR_PERM_COUNT = PTHREAD_MUTEX_INITIALIZER;

void *threadPermCompute(void * theArg)
{
	//Rprintf("threadPermCompute\n");
	//Rprintf("ARG_Dim[0] = %d\n", ARG_Dim[0]);
	//Rprintf("ARG_Dim[1] = %d\n", ARG_Dim[1]);
	//Rprintf("*ARG_NumberOfUniqueBatchIds = %d\n", *ARG_NumberOfUniqueBatchIds);
	int arrayIndex = 0;
	int arrayListIndex = 0;
	//
	const gsl_rng_type * T = gsl_rng_default;
	gsl_rng * glsRng = gsl_rng_alloc(T);
	//
	gsl_matrix * myMatrix = gsl_matrix_alloc(ARG_Dim[0], ARG_Dim[1]);
	copyMatrixFromArray(myMatrix, ARG_Values, ARG_Dim);
	gsl_vector * rowVector = gsl_vector_calloc(ARG_Dim[1]);
	gsl_vector * batchVectorList[*ARG_NumberOfUniqueBatchIds];
	gsl_vector * valuesForBatch[*ARG_NumberOfUniqueBatchIds];
	int x = 0;
	for(x=0;x<*ARG_NumberOfUniqueBatchIds;x++)
	{
		char * batchName = ARG_SortedUniqueBatchIds[x];
		int numSamplesInBatch = countOccurances(ARG_BatchIds, ARG_Dim[1], batchName);
		batchVectorList[x] = gsl_vector_calloc(numSamplesInBatch);
		valuesForBatch[x] = gsl_vector_calloc(numSamplesInBatch);
	}

	//Rprintf("mutex lock\n");
	pthread_mutex_lock( &LOCK_FOR_PERM_COUNT );
	if (NEXT_PERMUTATION>=NUM_PERMUTATIONS)
	{
		arrayIndex = -1;
	}
	else
	{
		arrayIndex = NEXT_PERMUTATION;
		NEXT_PERMUTATION = NEXT_PERMUTATION + 1;
	}
	pthread_mutex_unlock( &LOCK_FOR_PERM_COUNT );
	//Rprintf("mutex unlock\n");
	while(arrayIndex>=0)
	{
		gsl_rng_set(glsRng, arrayIndex+1);
		//Rprintf("loop arrayIndex=%d\n", arrayIndex);
		//
		arrayListIndex = arrayIndex * ARG_Dim[1];
		//Rprintf("loop arrayListIndex=%d\n", arrayListIndex);
		doPermutation(myMatrix, rowVector, batchVectorList, valuesForBatch,
				ARG_Values, ARG_Dim, ARG_BatchIds, ARG_SortedUniqueBatchIds, ARG_NumberOfUniqueBatchIds,
				&(ARG_retListOfGeneDSC[arrayListIndex]), &(ARG_retListOfGeneDB[arrayListIndex]), &(ARG_retListOfGeneDW[arrayListIndex]), 
				&(ARG_retDSC[arrayIndex]), &(ARG_retDB[arrayIndex]), &(ARG_retDW[arrayIndex]), glsRng);
		//Rprintf("DSCCPerm arrayIndex =%d\n" , arrayIndex);
		//Rprintf("DSCCPerm Results DSC =%.30f\n" , ARG_retDSC[arrayIndex]);
		//Rprintf("DSCCPerm Results DB =%.30f\n" , ARG_retDB[arrayIndex]);
		//Rprintf("DSCCPerm Results DW =%.30f\n" , ARG_retDW[arrayIndex]);
		////
		pthread_mutex_lock( &LOCK_FOR_PERM_COUNT );
		if (NEXT_PERMUTATION>=NUM_PERMUTATIONS)
		{
			arrayIndex = -1;
		}
		else
		{
			arrayIndex = NEXT_PERMUTATION;
			NEXT_PERMUTATION = NEXT_PERMUTATION + 1;
		}
		pthread_mutex_unlock( &LOCK_FOR_PERM_COUNT );
	}
	for(x=0;x<*ARG_NumberOfUniqueBatchIds;x++)
	{
		gsl_vector_free(batchVectorList[x]);
		gsl_vector_free(valuesForBatch[x]);
	}
	gsl_rng_free(glsRng);
	gsl_vector_free(rowVector);
	gsl_matrix_free(myMatrix);
}
