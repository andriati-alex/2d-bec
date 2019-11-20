#ifndef _observables_h
#define _observables_h

#include "calculus.h"




doublec Chem(int, int, double, double, double, double, double, Rarray,
        Rarray, Rarray, Carray);

doublec Energy(int, int, double, double, double, double, double, Rarray,
        Rarray, Rarray, Carray);



double MeanR(int,int,Carray,double,double,Rarray,Rarray);



double MaxResidue(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f, double mu);





#endif
