/* 
 * File:   c_threads.h
 * Author: tdcasasent
 *
 * Created on May 8, 2012, 10:46 AM
 */

#ifndef C_THREADS_H
#define	C_THREADS_H

#ifdef	__cplusplus
extern "C"
{
#endif

extern int NUM_PERMUTATIONS;
extern int NEXT_PERMUTATION;
extern double * ARG_Values;
extern int * ARG_Dim;
extern char ** ARG_BatchIds;
extern char ** ARG_SortedUniqueBatchIds;
extern int * ARG_NumberOfUniqueBatchIds;
extern double * ARG_retListOfGeneDSC;
extern double * ARG_retListOfGeneDB;
extern double * ARG_retListOfGeneDW;
extern double * ARG_retDSC;
extern double * ARG_retDB;
extern double * ARG_retDW;

extern void *threadPermCompute(void * theArg);
		
#ifdef	__cplusplus
}
#endif

#endif	/* C_THREADS_H */

