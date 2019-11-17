#ifndef _tridiagonalSystems_h
#define _tridiagonalSystems_h

#include "memoryHandling.h"
#include "arrayOperations.h"



double cond(int n, Carray upper, Carray lower, Carray mid);
double errBack(int n, Carray upper, Carray lower, Carray mid);
/* Compute condition number u
 * **************************
 *
 * See BUENO, M. Isabel & DOPICO, FROILÁN M. for the meaning of condition
 * numbers in "Stability and sensitivity of tridiagonal LU factorization"
 *
 * A good condition number is < 1E10 (assure at least 6 digits precision).
 *
 * The Backward error should be of 10E-15
 *
 * **********************************************************************/





double complex checkLU(int n, Carray upper, Carray lower, Carray mid);
/* Check if the tridiagonal system has L.U. factorization
 * ******************************************************
 *
 * n is the size of the system
 * upper is upper diagonal and has size (n-1)
 * lower is lower diagonal and has size (n-1)
 * mid is the main diagonal of size (n)
 *
 * The process ends up with the determinant, otherwise it
 * gives zero(a source of breakdown to LU factorization)
 *
 * ******************************************************/





void triDiag(int, doublec, doublec, doublec, Carray, Carray);

#endif
