#include "memoryHandling.h"



/* ========================================================================
 
                               MEMORY ALLOCATION                            

   ======================================================================== */



Rarray rarrDef(int n)
{
    double * ptr;

    ptr = (double * ) malloc( n * sizeof(double) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for double\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Carray carrDef(int n)
{
    double complex * ptr;

    ptr = (double complex * ) malloc( n * sizeof(double complex) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



CMKLarray cmklDef(int n)
{
    MKL_Complex16 * ptr;

    ptr = (MKL_Complex16 *) malloc( n * sizeof(MKL_Complex16) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex(mkl)\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Rmatrix rmatDef(int m, int n)
{

/** Real matrix of m rows and n columns **/

    int i;

    double ** ptr;

    ptr = (double ** ) malloc( m * sizeof(double *) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for real matrix\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = rarrDef(n);

    return ptr;
}



Cmatrix cmatDef(int m, int n)
{

/** Complex matrix of m rows and n columns **/

    int i;

    double complex ** ptr;

    ptr = (double complex ** ) malloc( m * sizeof(double complex *) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex matrix\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = carrDef(n);

    return ptr;
}



CCScmat ccscmatDef(int n, int max_nonzeros)
{

/** Return empty CCS representation of matrix of  n  rows
  * and the maximum number of nonzeros elements in a same
  * row given in 'max_nonzeros'                       **/

    CCScmat
        M;

    M = (struct _CCScmat *) malloc(sizeof(struct _CCScmat));

    if (M == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for CCS matrix\n\n");
        exit(EXIT_FAILURE);
    }

    M->m = max_nonzeros;
    M->vec = carrDef(max_nonzeros * n);
    M->col = (int *) malloc( max_nonzeros * n * sizeof(int) );

    if (M->col == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integers\n\n");
        exit(EXIT_FAILURE);
    }

    return M;
}



CCSrmat ccsrmatDef(int n, int max_nonzeros)
{

    CCSrmat
        M;

    M = (struct _CCSrmat *) malloc(sizeof(struct _CCSrmat));

    if (M == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for CCS matrix\n\n");
        exit(EXIT_FAILURE);
    }

    M->m = max_nonzeros;
    M->vec = rarrDef(max_nonzeros * n);
    M->col = (int *) malloc( max_nonzeros * n * sizeof(int) );

    if (M->col == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integers\n\n");
        exit(EXIT_FAILURE);
    }

    return M;
}



EqDataPkg PackEqData(int nx,int ny,double xi,double xf,double yi, double yf,
          double b,double g,double Ome,char Vname[],double p[])
{

/** Return pointer to a basic data structure with all needed information
  * to solve the Gross-Pitaevskii equation                           **/

    EqDataPkg gp = (EqDataPkg) malloc(sizeof(struct _EquationDataPkg));

    if (gp == NULL)
    {
        printf("\n\n\nMEMORY ERROR : malloc fail for EqData structure\n\n");
        exit(EXIT_FAILURE);
    }

    gp->x = rarrDef(nx);
    gp->y = rarrDef(ny);
    gp->V = rarrDef(nx*ny);

    gp->nx = nx;
    gp->ny = ny;
    gp->hx = (xf - xi) / (nx - 1);
    gp->hy = (yf - yi) / (ny - 1);
    gp->g = g;
    gp->b = b;
    gp->Ome = Ome;

    rarrFillInc(nx, xi, gp->hx, gp->x);
    rarrFillInc(ny, yi, gp->hy, gp->y);

    gp->p[0] = p[0];
    gp->p[1] = p[1];
    gp->p[2] = p[2];
    strcpy(gp->Vname,Vname);

    GetPotential(Vname,nx,ny,gp->x,gp->y,gp->V,p);

    return gp;
}



/* ========================================================================
 
                               MEMORY RELEASE

   ======================================================================== */



void rmatFree(int m, Rmatrix M)
{

/** Release a real matrix of m rows **/

    int i;

    for (i = 0; i < m; i++) free(M[i]);

    free(M);
}



void cmatFree(int m, Cmatrix M)
{

/** Release a complex matrix of m rows **/

    int i;

    for (i = 0; i < m; i++) free(M[i]);

    free(M);
}



void ccscmatFree(CCScmat M)
{

/** Release Compressed-Column Storaged matrix **/

    free(M->col);
    free(M->vec);
    free(M);
}



void ccsrmatFree(CCSrmat M)
{

/** Release Compressed-Column Storaged matrix **/

    free(M->col);
    free(M->vec);
    free(M);
}



void ReleaseEqDataPkg(EqDataPkg gp)
{
    free(gp->V);
    free(gp->x);
    free(gp->y);
    free(gp);
}
