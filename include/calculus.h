#ifndef _calculus_h
#define _calculus_h

#include <mkl.h>
#include <mkl_dfti.h>
#include <math.h>

#include "memoryHandling.h"
#include "arrayOperations.h"

double complex Csimps1D(int, Carray, double);

double Rsimps1D(int, Rarray, double);

double complex Csimps2D(int, int, Carray, double, double);

double Rsimps2D(int, int, Rarray, double, double);

void DfDx(int, int, Carray, double, Carray);

void DfDy(int, int, Carray, double, Carray);

void renormalize(int, int, Carray, double, double, double);





#endif
