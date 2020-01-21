#include "tridiagonalSystems.h"



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



void triDiagLU(int n, Carray l, Carray u, doublec upper, Carray RHS, Carray ans)
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
