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



void quartic(int nx, int ny, Rarray x, Rarray y, Rarray V, double wx, double wy)
{
    int
        i,
        j;

    double
        x4,
        y4;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x4 = x[i] * x[i] * x[i] * x[i];
            y4 = y[j] * y[j] * y[j] * y[j];
            V[i + j*nx] = 0.5 * (wx*wx*x4 + wy*wy*y4);
        }
    }
}



void GaussianRing(int nx, int ny, Rarray x, Rarray y, Rarray V,
     double wx, double wy, double height, double width)
{
    int
        i,
        j;

    double
        x2,
        y2;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x2 = x[i] * x[i];
            y2 = y[j] * y[j];
            V[i + j*nx] = 0.5 * (wx*wx*x2 + wy*wy*y2) + \
                          height * exp(-(x2 + y2) / width / width);
        }
    }
}



void GetPotential(char name [], int nx, int ny, Rarray x, Rarray y, Rarray V,
                  double p [])
{

    if (strcmp(name, "zero") == 0)
    {
        rarrFill(nx*ny,0,V);
        return;
    }

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

    if (strcmp(name, "quartic") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0)
        {
            printf("\n\nQuartic Trap parameters must be positive\n");
            exit(EXIT_FAILURE);
        }

        quartic(nx,ny,x,y,V,p[0],p[1]);
        return;
    }

    if (strcmp(name, "GaussianRing") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0 || p[2] <= 0)
        {
            printf("\n\nHarmonic parameters must be positive and");
            printf(" gaussian height must be positive.\n");
            exit(EXIT_FAILURE);
        }

        GaussianRing(nx,ny,x,y,V,p[0],p[1],p[2],p[3]);
        return;
    }

    printf("\n\n\nERROR: Potential '%s' not implemented\n\n", name);
    exit(EXIT_FAILURE);
}
