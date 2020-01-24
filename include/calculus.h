#ifndef _calculus_h
#define _calculus_h

#include <math.h>

#include "memoryHandling.h"
#include "arrayOperations.h"

double complex Csimps1D(int, Carray, double);

double Rsimps1D(int, Rarray, double);

double complex Csimps2D(int, int, Carray, double, double);

double Rsimps2D(int, int, Rarray, double, double);

void DfDx(int, int, Carray, double, Carray);

void DfDx_real(int, int, Rarray, double, Rarray);

void DfDy(int, int, Carray, double, Carray);

void DfDy_real(int, int, Rarray, double, Rarray);

void renormalize(int, int, Carray, double, double, double);

void renormalizeReal(int, int, Rarray, double, double, double);





#endif
