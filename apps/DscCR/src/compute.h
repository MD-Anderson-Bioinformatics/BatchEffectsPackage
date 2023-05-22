/* 
 * File:   compute.h
 * Author: tdcasasent
 *
 * Created on May 2, 2012, 9:38 AM
 */

#ifndef COMPUTE_H
#define	COMPUTE_H

#ifdef	__cplusplus
extern "C"
{
#endif

extern void displayVector(gsl_vector * theVector, int theLimit);
extern void compute(gsl_matrix * theMatrix, gsl_vector * theRowVector,
		gsl_vector ** theBatchVectorList, gsl_vector ** theValuesForBatch, 
		int theNumberOfSamples, char ** theBatchIds, int theNumberOfBatchIds, char ** theSortedUniqueBatchIds, int theNumberOfUniqueBatchIds,
		double * retListOfGeneDSC, double * retListOfGeneDB, double * retListOfGeneDW, 
		double * retDSC, double * retDB, double * retDW);


#ifdef	__cplusplus
}
#endif

#endif	/* COMPUTE_H */

