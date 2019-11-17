#ifndef _dataStructures_h
#define _dataStructures_h

#include <complex.h>
#include <mkl.h>



/* DATATYPES OF THE PACKAGE
 * *********************************************************
 *
 * Description
 * ===========
 *
 * Header with all datatypes used on this project
 *
 *
 *
 * The Compressed-Column-Storage(CCS) of a Matrix
 * ==============================================
 *
 * This method consist in save memory when  storing a
 * sparse Matrix by shifting all the nonzero elements
 * to the left. Then the resulting matrix has only  m
 * columns,  m being  the maximum  number of  nonzero
 * elements in a same row, and the same number of rows
 * than the former one.  Thereafter  this  matrix  is
 * organized in a  vector  following  a  Column Major
 * ordering as well as its former column positions are
 * stored in a vector of integers of the same size.
 *
 *
 *
 * The Equation Package structure
 * ==============================================
 *
 * It is created a structure to hold all equation's coefficients
 * and external potential (usually trap),  for simplification on
 * the integrators API.
 *
 * *************************************************************/





#define PI 3.141592653589793



typedef double complex doublec;



typedef MKL_Complex16 * CMKLarray;



typedef double *  Rarray;
typedef double ** Rmatrix;



typedef double complex *  Carray;
typedef double complex ** Cmatrix;



/* Compressed Column storage of an n x n Real Matrix
 * ================================================= */

struct _CCScmat
{
    int  m;     // max number of non-zero elements in a same row
	int * col;  // Column index of elemetns.
	Carray vec; // column oriented vector.
};

typedef struct _CCScmat * CCScmat;



/* Compressed Column storage of an n x n Real Matrix
 * ================================================= */

struct _CCSrmat
{
    // same fields from complex case above
    int  m;
	int * col;
	Rarray vec;
};

typedef struct _CCSrmat * CCSrmat;



/* Structure to store all relevant information to call integrator
 * ============================================================== */

struct _EquationDataPkg
{

    char
        Vname[80]; // One-Body external potential name

    int 
        nx,
        ny;

    double
        hx,     // grid spacing in x-direction
        hy,     // grid spacing in y-direction
        b,      // factor of second order derivatives
        g,      // contact interaction strength
        Ome;    // rotation frequency

    double
        p[3];   // Parameters to generate external potential

    Rarray
        x,
        y,
        V;

};

typedef struct _EquationDataPkg * EqDataPkg;

#endif
