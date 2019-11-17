#ifndef _imagtimeIntegrator_h
#define _imagtimeIntegrator_h

#include "tridiagonalSystems.h"
#include "arrayOperations.h"
#include "observables.h"
#include "inout.h"

void ExplicitY(int, int, int, doublec, double, double, double, Rarray, Carray,
     Carray);

void ExplicitX(int, int, int, doublec, double, double, double, Rarray, Carray,
     Carray);

int SplitStepPR(EqDataPkg, int, double, Carray);

#endif
