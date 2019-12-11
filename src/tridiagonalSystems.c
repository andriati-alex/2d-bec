#include "tridiagonalSystems.h"



/*          ***********************************************

                 CHECK SOLVABILITY of LU DECOMPOSITION

            ***********************************************          */



double cond(int n, Carray upper, Carray lower, Carray mid)
{

/** Indicator of accuracy for solution of tridiagonal systems **/

    unsigned int i;

    // condition numbers

    double
        maxCond,
        condK;

    // L . U decomposition

    double complex
        u,
        l;

    maxCond = 1;
    condK = 1;

    u = mid[0];

    for (i = 0;  i < n - 1; i++)
    {
        l = lower[i] / u;

        u = mid[i + 1] - l * upper[i];

        condK = 1 + cabs(upper[i] * l / u) * (2 + condK);

        if (maxCond < condK) maxCond = condK;
    }

    return maxCond;
}





double errBack(int n, Carray upper, Carray lower, Carray mid)
{

/** Indicator of accuracy for solution of tridiagonal systems **/

    unsigned int i;

    // Backward error

    double
        maxErrBack,
           ErrBack;

    // L . U decomposition

    double complex
        u,
        l;

    maxErrBack = 1;
    ErrBack = 1;

    u = mid[0];
    for (i = 0;  i < n - 1; i++)
    {
        l = lower[i] / u;
        u = mid[i+1] - l * upper[i];

        ErrBack = (cabs(l) * cabs(upper[i]) + cabs(u)) / cabs(mid[i+1]);

        if (maxErrBack < ErrBack) maxErrBack = ErrBack;
    }

    return maxErrBack * 1E-16;
}





double complex checkLU(int n, Carray upper, Carray lower, Carray mid)
{

/** Check using the three diagonals the solvability by LU decomposition
  * In positive case return the determinant  of the tridiagonal matrix,
  * and return 0 otherwise. **/

    unsigned int i;

    double complex
        det0,
        det1,
        detK;

    det0 = 1;
    det1 = mid[0];

    for (i = 1; i < n; i++)
    {
        detK = mid[i] * det1 - lower[i-1] * upper[i-1] * det0;

        if (cabs(detK) == 0) return 0;

        det0 = det1;
        det1 = detK;
    }

    return detK;
}



void triDiag(int n, doublec upper, doublec lower, doublec mid, Carray RHS,
     Carray ans)
{

    unsigned int
        i,
        k;

    // Intermediate steps - L . U decomposition

    Carray
        u,
        l,
        z;



    u = carrDef(n);
    l = carrDef(n - 1);
    z = carrDef(n);

    u[0] = mid;
    z[0] = RHS[0];

    for (i = 0;  i < n - 1; i++)
    {
        k = i + 1;
        l[i] = lower / u[i];
        u[k] = mid - l[i] * upper;
        z[k] = RHS[k] - l[i] * z[i];
    }

    ans[n-1] = z[n-1] / u[n-1];

    for (i = 2; i <= n; i++)
    {
        k = n - i;
        ans[k] = (z[k] - upper * ans[k+1]) / u[k];
    }

    // Free local allocated memory
    free(u);
    free(l);
    free(z);
}





void LU(int n, doublec upper, doublec lower, doublec mid, Carray l, Carray u)
{

    unsigned int
        i,
        k;

    u[0] = mid;

    for (i = 0;  i < n - 1; i++)
    {
        k = i + 1;
        l[i] = lower / u[i];
        u[k] = mid - l[i] * upper;
    }

}





void triDiagLU(int n, Carray l, Carray u, doublec upper, Carray RHS,
	       Carray ans)
{

    unsigned int
        i,
        k;

    Carray
        z;



    z = carrDef(n);

    z[0] = RHS[0];

    for (i = 0;  i < n - 1; i++)
    {
        k = i + 1;
        z[k] = RHS[k] - l[i] * z[i];
    }

    ans[n-1] = z[n-1] / u[n-1];

    for (i = 2; i <= n; i++)
    {
        k = n - i;
        ans[k] = (z[k] - upper * ans[k+1]) / u[k];
    }

    free(z);
}
