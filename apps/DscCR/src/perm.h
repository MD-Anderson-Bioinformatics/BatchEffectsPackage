/* 
 * File:   perm.h
 * Author: tdcasasent
 *
 * Created on May 3, 2012, 12:28 PM
 */

#ifndef PERM_H
#define	PERM_H

#ifdef	__cplusplus
extern "C"
{
#endif

extern void copyMatrixFromArray(gsl_matrix * theMatrix, double * theMatrixValues, int * theMatrixDim);

extern void doPermutation(gsl_matrix * theMatrix, gsl_vector * theRowVector,  
		gsl_vector ** theBatchVectorList, gsl_vector ** theValuesForBatch, 
		double * theValues, int * theDim, char ** theBatchIds, 
		char ** theSortedUniqueBatchIds, int * theNumberOfUniqueBatchIds,
		double * retListOfGeneDSC, double * retListOfGeneDB, double * retListOfGeneDW, 
		double * retDSC, double * retDB, double * retDW, gsl_rng * theGlsRng);


#ifdef	__cplusplus
}
#endif

#endif	/* PERM_H */

