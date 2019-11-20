#include "observables.h"





doublec Chem(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f)
{

/** Gross-Pitaesvkii functional
  * ----------------------------
  *
  *  To compute the energy the value passed on g must be the contact
  *  interaction strength divided by 2, otherwise gives the chemical
  *  potential of the system. It divides by the norm of the function
  *  in the end to enforce the result obtained is the particle average
  *
  *  Parameters
  *
  *  (nx,ny) - number of grid points in each direction
  *  (hx,hy) - grid step in each direction
  *  b - number to multiply second order derivatives
  *  Ome - Angular rotation frequency
  *  g - see description above
  *  V - Potential values at grid points
  *  f - function
  *
**/

    int
        j,
        i;

    double
        norm;

    double complex
        E;

    Carray
        dfdx,
        dfdy,
        Integ;

    Rarray
        abs2f,
        abs2dfdx,
        abs2dfdy;



    dfdx = carrDef(nx*ny);
    dfdy = carrDef(nx*ny);
    Integ = carrDef(nx*ny);
    abs2f = rarrDef(nx*ny);
    abs2dfdx = rarrDef(nx*ny);
    abs2dfdy = rarrDef(nx*ny);



    DfDx(nx,ny,f,hx,dfdx);
    DfDy(nx,ny,f,hy,dfdy);

    carrAbs2(nx*ny,dfdx,abs2dfdx);
    carrAbs2(nx*ny,dfdy,abs2dfdy);

    carrAbs2(nx*ny,f,abs2f);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            // Sum potential part
            E = (V[i + j*nx] + g * abs2f[i + j*nx]) * abs2f[i + j*nx];
            // Sum kinetic part
            E = E - b * (abs2dfdx[i + j*nx] + abs2dfdy[i + j*nx]);
            // Sum rotation part
            E = E + I * Ome * conj(f[i + j*nx]) * (x[i]*dfdy[i + j*nx] \
                    - y[j]*dfdx[i + j*nx]);
            // Setup the integrand
            Integ[i + j*nx] = E;
        }
    }

    E = Csimps2D(nx,ny,Integ,hx,hy);
    norm = Rsimps2D(nx,ny,abs2f,hx,hy);

    // release memory
    free(dfdx); free(dfdy); free(abs2f); free(abs2dfdx); free(abs2dfdy);
    free(Integ);

    return E / norm;

}





doublec Energy(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f)
{
    return Chem(nx,ny,hx,hy,b,Ome,g/2,V,x,y,f);
}





double MeanR(int nx, int ny, Carray f, double hx, double hy,
       Rarray x, Rarray y)
{

/** Compute Mean Square value of normalized complex function/distribution **/

    int
        j,
        i;

    double
        r;

    Rarray
        Integ;

    Integ = rarrDef(nx*ny);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            r = sqrt(x[i]*x[i] + y[j]*y[j]);
            Integ[i + j*nx] = r * ( creal(f[i + j*nx])*creal(f[i + j*nx]) + \
                       cimag(f[i + j*nx])*cimag(f[i + j*nx]) );
        }
    }

    r = Rsimps2D(nx,ny,Integ,hx,hy);

    free(Integ);

    return r;

}





double MaxResidue(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f, double mu)
{

    int
        j,
        i;

    double
        sentinel,
        maxres;

    double complex
        Hf;

    Carray
        dfdx,
        dfdy,
        d2fdx2,
        d2fdy2;

    Rarray
        abs2f;



    dfdx = carrDef(nx*ny);
    dfdy = carrDef(nx*ny);
    d2fdx2 = carrDef(nx*ny);
    d2fdy2 = carrDef(nx*ny);
    abs2f = rarrDef(nx*ny);



    carrAbs2(nx*ny,f,abs2f);

    DfDx(nx,ny,f,hx,dfdx);
    DfDy(nx,ny,f,hy,dfdy);

    DfDx(nx,ny,dfdx,hx,d2fdx2);
    DfDy(nx,ny,dfdy,hy,d2fdy2);

    maxres = 0;

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            // Sum potential part
            Hf = (V[i + j*nx] + g * abs2f[i + j*nx]) * f[i + j*nx];
            // Sum kinetic part
            Hf = Hf + b * (d2fdx2[i + j*nx] + d2fdy2[i + j*nx]);
            // Sum rotation part
            Hf = Hf + I * Ome * (x[i]*dfdy[i + j*nx] - y[j]*dfdx[i + j*nx]);

            sentinel = cabs(Hf - mu * f[i + j*nx]);

            if (sentinel > maxres) maxres = sentinel;
        }
    }

    // release memory
    free(dfdx); free(dfdy); free(abs2f); free(d2fdx2); free(d2fdy2);

    return maxres;

}





double AvgResidue(int nx, int ny, double hx, double hy, double b, double Ome,
        double g, Rarray V, Rarray x, Rarray y, Carray f, double mu)
{

    int
        j,
        i;

    double
        res;

    double complex
        Hf;

    Carray
        dfdx,
        dfdy,
        d2fdx2,
        d2fdy2;

    Rarray
        integ,
        abs2f;



    dfdx = carrDef(nx*ny);
    dfdy = carrDef(nx*ny);
    d2fdx2 = carrDef(nx*ny);
    d2fdy2 = carrDef(nx*ny);
    abs2f = rarrDef(nx*ny);
    integ = rarrDef(nx*ny);



    carrAbs2(nx*ny,f,abs2f);

    DfDx(nx,ny,f,hx,dfdx);
    DfDy(nx,ny,f,hy,dfdy);

    DfDx(nx,ny,dfdx,hx,d2fdx2);
    DfDy(nx,ny,dfdy,hy,d2fdy2);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            // Sum potential part
            Hf = (V[i + j*nx] + g * abs2f[i + j*nx]) * f[i + j*nx];
            // Sum kinetic part
            Hf = Hf + b * (d2fdx2[i + j*nx] + d2fdy2[i + j*nx]);
            // Sum rotation part
            Hf = Hf + I * Ome * (x[i]*dfdy[i + j*nx] - y[j]*dfdx[i + j*nx]);

            integ[i + j*nx] = cabs(Hf - mu * f[i + j*nx]);
        }
    }

    res = Rsimps2D(nx,ny,integ,hx,hy);

    // release memory
    free(dfdx); free(dfdy); free(abs2f); free(d2fdx2); free(d2fdy2);
    free(integ);

    return res;

}
