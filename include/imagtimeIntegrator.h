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

void ExplicitX_alongX(int nx, int ny, doublec dt, double hx, double b,
     double Ome, double yj, Carray in, Carray out);

int SplitStepPR(EqDataPkg, int, double, Carray);

int SplitStepDYakonov(EqDataPkg EQ, int N, double realDT, Carray S);

#endif
