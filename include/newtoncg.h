#ifndef _newtoncg_h
#define _newtoncg_h

#include <mkl.h>
#include <mkl_dfti.h>
#include "inout.h"
#include "calculus.h"
#include "observables.h"
#include "arrayOperations.h"

double maxNorm(int N, Rarray fr, Rarray fi);

void StationaryOp(EqDataPkg,double,Rarray,Rarray,Rarray,Rarray);

void linearizedOp(EqDataPkg,Rarray,Rarray,double,Rarray,Rarray,Rarray,Rarray);

int conjgrad(EqDataPkg,double,double,Rarray,Rarray,Rarray,Rarray,Rarray,Rarray);

void stationaryNewton(EqDataPkg,Carray,double);

#endif
