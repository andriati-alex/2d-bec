#ifndef _tridiagonalSystems_h
#define _tridiagonalSystems_h

#include "memoryHandling.h"
#include "arrayOperations.h"

void triDiag(int, doublec, doublec, doublec, Carray, Carray);
void LU(int, doublec, doublec, doublec, Carray, Carray);
void triDiagLU(int, Carray, Carray, doublec, Carray, Carray);

#endif
