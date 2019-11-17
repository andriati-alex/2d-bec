#include "calculus.h"





double complex Csimps1D(int n, Carray f, double h)
{

    int
        i;

    double complex
        sum;

    sum = 0;

    if (n < 3)
    {
        printf("\n\n\tERROR : less than 3 point to integrate by simps !\n\n");
        exit(EXIT_FAILURE);
    }

    if (n % 2 == 0)
    {

    //  Case the number of points is even then must integrate the last
    //  chunk using simpson's 3/8 rule to maintain accuracy

        for (i = 0; i < (n - 4); i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals
        sum = sum + (f[n-4] + 3 * (f[n-3] + f[n-2]) + f[n-1]) * 3 * h / 8;

    }

    else
    {

        for (i = 0; i < n - 2; i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals

    }

    return sum;

}





double Rsimps1D(int n, Rarray f, double h)
{

    int
        i;

    double
        sum;

    sum = 0;

    if (n % 2 == 0)
    {

    //  Case the number of points is even then must integrate the last
    //  chunk using simpson's 3/8 rule to maintain accuracy

        for (i = 0; i < (n - 4); i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals
        sum = sum + (f[n-4] + 3 * (f[n-3] + f[n-2]) + f[n-1]) * 3 * h / 8;

    }

    else
    {

        for (i = 0; i < n - 2; i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals

    }

    return sum;

}





double complex Csimps2D(int nx, int ny, Carray f, double hx, double hy)
{

/** INTEGRATION OF A FUNCTION OF 2 VARIABLES
  * ----------------------------------------
  *
  * Here the function must be stored in a vector, where given point a
  * point in the grid (xi,yj) then
  *
  *     funtion(xi,yj) = f[i + nx * j]
  *
  * where i = 0, 1, ..., nx-1 and j = 0, 1, ..., ny-1 and the grid points
  * are linearly spaced xi = x0 + i*hx as well as yj = y0 + j*hy      **/

    unsigned int
        j;

    double complex
        result;

    Carray
        fy;

    fy = carrDef(ny);

    for (j = 0; j < ny; j++)
    {
        // Integrate in x-direction and end up with a function of y
        fy[j] = Csimps1D(nx,&f[j*nx],hx);
    }

    result = Csimps1D(ny,fy,hy);

    free(fy);

    return result;

}





double Rsimps2D(int nx, int ny, Rarray f, double hx, double hy)
{

/** INTEGRATION OF A FUNCTION OF 2 VARIABLES
  * ----------------------------------------
  *
  * Here the function must be stored in a vector, where given point a
  * point in the grid (xi,yj) then
  *
  *     funtion(xi,yj) = f[i + nx * j]
  *
  * where i = 0, 1, ..., nx-1 and j = 0, 1, ..., ny-1 and the grid points
  * are linearly spaced xi = x0 + i*hx as well as yj = y0 + j*hy      **/

    unsigned int
        j;

    double
        result;

    Rarray
        fy;

    fy = rarrDef(ny);

    for (j = 0; j < ny; j++)
    {
        // Integrate in x-direction and end up with a function of y
        fy[j] = Rsimps1D(nx,&f[j*nx],hx);
    }

    result = Rsimps1D(ny,fy,hy);

    free(fy);

    return result;

}





void renormalize(int nx, int ny, Carray f, double hx, double hy, double norm)
{

    int
        i;

    double
        renorm;

    Rarray
        ModSquared;

    ModSquared = rarrDef(nx*ny);

    carrAbs2(nx*ny,f,ModSquared);

    renorm = norm * sqrt(1.0 / Rsimps2D(nx,ny,ModSquared,hx,hy));

    for (i = 0; i < nx*ny; i++) f[i] = f[i] * renorm;

    free(ModSquared);
}





void DfDx(int nx, int ny, Carray f, double hx, Carray dfdx)
{

/** Compute partial derivative in x-direction
  *
  * Output parameter : dfdx is a  function of the
  * grid points df/dx(xi,yi) = dfdx[i + j*nx] **/

    int
        i,
        j,
        s;

    double
        r;

    r = 1.0 / (12 * hx); // ratio for a fourth-order scheme

    for (j = 0; j < ny; j++)
    {
        s = j * nx; // stride to get the yj grid line

        // Do separately the computation near  the boundary
        // since involves evaluation where the  'f' is zero
        // The zeros are kept to track the rule in the loop
        // every time the loop yield an evaluation  out  of
        // the domain the function is assumed to give  zero

        dfdx[0+s]   = ( 0 - f[2+s] + 8 * (f[1+s] - 0) ) * r;

        dfdx[1+s]   = ( 0 - f[3+s] + 8 * (f[2+s] - f[0+s]) ) * r;

        dfdx[nx-2+s] = ( f[nx-4+s] - 0 + 8 * (f[nx-1+s] - f[nx-3+s]) ) * r;

        dfdx[nx-1+s] = ( f[nx-3+s] - 0 + 8 * (0 - f[nx-2+s]) ) * r;

        for (i = 2; i < nx - 2; i++)
        {
            dfdx[i+s] = ( f[i-2+s] - f[i+2+s] + 8*(f[i+1+s] - f[i-1+s]) )*r;
        }
    }

}





void DfDy(int nx, int ny, Carray f, double hy, Carray dfdy)
{

/** Compute partial derivative in x-direction
  *
  * Output parameter : dfdy is a  function of the
  * grid points df/dy(xi,yi) = dfdy[i + j*nx] **/

    int
        i,
        j,
        s,
        s1,
        s2,
        sm1,
        sm2;

    double
        r;

    r = 1.0 / (12 * hy); // ratio for a fourth-order scheme

    for (i = 0; i < nx; i++)
    {

        // Do separately the computation near  the boundary
        // since involves evaluation where the  'f' is zero
        // The zeros are kept to track the rule in the loop
        // every time the loop yield an evaluation  out  of
        // the domain the function is assumed to give  zero

        s = 0*nx;
        s1  = 1*nx;
        s2  = 2*nx;
        dfdy[i+s]   = ( 0 - f[i+s2] + 8 * (f[i+s1] - 0) ) * r;

        s = 1*nx;
        s1  = 2*nx;
        s2  = 3*nx;
        sm1 = 0*nx;
        dfdy[i+s]   = ( 0 - f[i+s2] + 8 * (f[i+s1] - f[i+sm1]) ) * r;

        s = (ny-2)*nx;
        s1  = (ny-1)*nx;
        sm1 = (ny-3)*nx;
        sm2 = (ny-4)*nx;
        dfdy[i+s] = ( f[i+sm2] - 0 + 8 * (f[i+s1] - f[i+sm1]) ) * r;

        s = (ny-1)*nx;
        sm1 = (ny-2)*nx;
        sm2 = (ny-3)*nx;
        dfdy[i+s] = ( f[i+sm2] - 0 + 8 * (0 - f[i+sm1]) ) * r;

        for (j = 2; j < ny - 2; j++)
        {
            s = j * nx; // stride to get the yj point
            s1  = (j + 1) * nx;
            s2  = (j + 2) * nx;
            sm1 = (j - 1) * nx;
            sm2 = (j - 2) * nx;

            dfdy[i+s] = ( f[i+sm2] - f[i+s2] + 8*(f[i+s1] - f[i+sm1]) )*r;
        }
    }

}
