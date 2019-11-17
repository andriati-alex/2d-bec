#include "linearPotential.h"



void harmonic(int nx, int ny, Rarray x, Rarray y, Rarray V, double wx, double wy)
{
    int
        i,
        j;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            V[i + j*nx] = 0.5 * (wx*wx*x[i]*x[i] + wy*wy*y[j]*y[j]);
        }
    }
}



void GetPotential(char name [], int nx, int ny, Rarray x, Rarray y, Rarray V,
                  double p [])
{
    if (strcmp(name, "harmonic") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0)
        {
            printf("\n\nHarmonic Trap parameters must be positive\n");
            exit(EXIT_FAILURE);
        }

        harmonic(nx,ny,x,y,V,p[0],p[1]);
        return;
    }

    printf("\n\n\nERROR: Potential '%s' not implemented\n\n", name);
    exit(EXIT_FAILURE);
}
