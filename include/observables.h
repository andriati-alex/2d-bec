#ifndef _observables_h
#define _observables_h

#include "calculus.h"




doublec Chem(int, int, double, double, double, double, double, Rarray,
        Rarray, Rarray, Carray);



doublec Energy(int, int, double, double, double, double, double, Rarray,
        Rarray, Rarray, Carray);

doublec angularMom(int, int, double, double, Rarray, Rarray, Carray);



double Kinect(int nx, int ny, double hx, double hy, double b, Carray f);



double Potential(int nx, int ny, double hx, double hy, Rarray V, Carray f);



double Interacting(int nx, int ny, double hx, double hy, double g, Carray f);



double Virial(int nx, int ny, double hx, double hy, double b, double g,
       Rarray V, char Vname [], Carray f);



double MeanR(int,int,Carray,double,double,Rarray,Rarray);



double MaxResidue(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f, double mu);



double AvgResidue(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f, double mu);





#endif
