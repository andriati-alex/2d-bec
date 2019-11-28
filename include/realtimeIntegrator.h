#ifndef _realtimeIntegrator_h
#define _realimeIntegrator_h

#include "imagtimeIntegrator.h"

int RealSplitStepPR(EqDataPkg,int,double,Carray,int,char []);

int RealSplitStepDYakonov(EqDataPkg,int,double,Carray,int, char []);

#endif
