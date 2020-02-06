#include "newtoncg.h"



void checkFFTstatus(MKL_LONG status)
{
    if (status != 0)
    {
        printf("\n\nERROR IN FFT ROUTINE : status returned %d\n",status);
        printf("Error message given : %s\n\n",DftiErrorMessage(status));
        exit(EXIT_FAILURE);
    }
}





void compressComplex(int nx, int ny, Rarray re, Rarray im, Carray f)
{

/** AUXILIAR ROUTINE
    ****************
    Transfer data from two arrays of real number corresponding to real
    and imaginary part of complex numbers of a function of 2 variables
    to a complex array

    Output parameters : 'f' with f[j*nx + i] = f(xi,yj)
    Input parameters  : 're' and 'im' matrices in row-major format
                        re[i + j*nx] = Re(f(xi,yj))
                        im[i + j*nx] = Im(f(xi,yj))
*********************************************************************/

    int
        i,
        j;

    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++) f[j*nx + i] = re[i+j*nx] + I * im[i+j*nx];
    }
}





void extractComplex(int nx, int ny, Rarray re, Rarray im, Carray f)
{

/** AUXILIAR ROUTINE
    ****************
    Invert the input and output parameters from compressComplexMat routine
    For more information see the description there
*************************************************************************/

    int
        i,
        j;

    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            re[i+j*nx] = creal(f[j*nx + i]);
            im[i+j*nx] = cimag(f[j*nx + i]);
        }
    }
}





double maxNorm(int N, Rarray fr, Rarray fi)
{

/** AUXILIAR ROUTINE
    ****************
    Given a function split in two arrays with its real and imaginary  parts
    compute the Maximum of complex absolute value of its elements. Works as
    norm to stop iterative methods
*************************************************************************/

    int
        i;

    double
        maxRes;

    maxRes = 0;

    for (i = 0; i < N; i++)
    {
        if (sqrt(fr[i]*fr[i] + fi[i]*fi[i]) > maxRes)
        {
            maxRes = sqrt(fr[i]*fr[i] + fi[i]*fi[i]);
        }
    }

    return maxRes;
}





void applyCond(int nx, int ny, double Lx, double Ly, double mu, double b,
               Carray f)
{
/** OPTIONAL ROUTINE
    ****************
    Routine to apply pre-conditioning using second order derivatives which
    appears in the equation as b*(d^2/dx^2 + d^2/dy^2).  Not much  helpful
    if the equation contain others terms like  trap potential and rotation

    Input Parameters  : nx & ny are grid dimension in each direction
                        Lx & Ly the domain extension
                        mu & b equation coefficients
                        f function values in the grid f[j*nx + i] = f(xi,yj)

    Output Parameters : f (overwritten)
**************************************************************************/

    int
        i,
        j;

    double
        freqx,
        freqy,
        bfy2,
        bfx2;

    MKL_LONG
        dim_sizes[2];

    MKL_LONG
        st;

    DFTI_DESCRIPTOR_HANDLE
        desc;



    dim_sizes[0] = ny;
    dim_sizes[1] = nx;

    st = DftiCreateDescriptor(&desc,DFTI_DOUBLE,DFTI_COMPLEX,2,dim_sizes);
    checkFFTstatus(st);
    st = DftiSetValue(desc,DFTI_FORWARD_SCALE,1.0/sqrt(nx*ny));
    checkFFTstatus(st);
    st = DftiSetValue(desc,DFTI_BACKWARD_SCALE,1.0/sqrt(nx*ny));
    checkFFTstatus(st);
    st = DftiCommitDescriptor(desc);
    checkFFTstatus(st);

    st = DftiComputeForward(desc, f);
    checkFFTstatus(st);

    for (j = 0; j < ny; j++)
    {
        if (j <= (ny - 1) / 2) { freqy = j / Ly;        }
        else                   { freqy = (j - ny) / Ly; }

        bfy2 = b * ( - 4 * PI * PI * freqy * freqy);

        for (i = 0; i < nx; i++)
        {
            if (i <= (nx - 1) / 2) { freqx = i / Lx;        }
            else                   { freqx = (i - nx) / Lx; }

            bfx2 = b * ( - 4 * PI * PI * freqx * freqx);

            f[j*nx + i] = f[j*nx + i] / (- mu + bfx2 + bfy2);
        }
    }

    st = DftiComputeBackward(desc, f);
    checkFFTstatus(st);
    st = DftiFreeDescriptor(&desc);
    checkFFTstatus(st);
}





void stationaryOp(EqDataPkg EQ, double mu, Rarray fr, Rarray fi, Rarray Fr,
                  Rarray Fi)
{

/** MAIN ROUTINE - Nonlinead operators that define the differential equation
    ************************************************************************
    Act with the nonlinear operator on a function given through its real and
    imaginary parts separately. If the given function is a  solution for the
    equation then it must give zero (in all grid points)

    Input Parameters  : Equation structure and 'lagrange multiplier' mu
                        fr and fi real and imag part of the input function

    Output Parameters : Fr and Fi real and imag part of the function after
                        the action of the operator
****************************************************************************/

    int
        i,
        j,
        nx,
        ny;

    double
        local,
        Ome,
        fr2,
        fi2,
        hx,
        hy,
        g,
        b;

    Rarray
        x,
        y,
        V,
        fr_dx,
        fi_dx,
        fr_dy,
        fi_dy,
        fr_dxdx,
        fi_dxdx,
        fr_dydy,
        fi_dydy;

    b = EQ->b;
    g = EQ->g;
    hx = EQ->hx;
    hy = EQ->hy;
    Ome = EQ->Ome;
    V = EQ->V;
    x = EQ->x;
    y = EQ->y;

    nx = EQ->nx;
    ny = EQ->ny;

    fr_dx = rarrDef(nx*ny);
    fi_dx = rarrDef(nx*ny);
    fr_dy = rarrDef(nx*ny);
    fi_dy = rarrDef(nx*ny);
    fr_dxdx = rarrDef(nx*ny);
    fi_dxdx = rarrDef(nx*ny);
    fr_dydy = rarrDef(nx*ny);
    fi_dydy = rarrDef(nx*ny);

    DfDx_real(nx,ny,fr,hx,fr_dx);
    DfDx_real(nx,ny,fi,hx,fi_dx);

    DfDy_real(nx,ny,fr,hy,fr_dy);
    DfDy_real(nx,ny,fi,hy,fi_dy);

    D2fDx2_real(nx,ny,fr,hx,fr_dxdx);
    D2fDx2_real(nx,ny,fi,hx,fi_dxdx);

    D2fDy2_real(nx,ny,fr,hy,fr_dydy);
    D2fDy2_real(nx,ny,fi,hy,fi_dydy);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fr2 = fr[i+j*nx]*fr[i+j*nx];
            fi2 = fi[i+j*nx]*fi[i+j*nx];

            local = - mu +  V[i+j*nx] + g * (fr2 + fi2);

            Fr[i + j*nx] = b * (fr_dxdx[i+j*nx] + fr_dydy[i+j*nx])  - \
                    Ome * (x[i]*fi_dy[i+j*nx] - y[j]*fi_dx[i+j*nx]) + \
                    local * fr[i+j*nx];

            Fi[i + j*nx] = b * (fi_dxdx[i+j*nx] + fi_dydy[i+j*nx])  + \
                    Ome * (x[i]*fr_dy[i+j*nx] - y[j]*fr_dx[i+j*nx]) +
                    local * fi[i+j*nx];
        }
    }

    free(fr_dx);
    free(fr_dy);
    free(fi_dx);
    free(fi_dy);
    free(fr_dxdx);
    free(fr_dydy);
    free(fi_dxdx);
    free(fi_dydy);
}





void linearizedOp(EqDataPkg EQ, Rarray phir, Rarray phii, double mu,
                  Rarray infr, Rarray infi, Rarray outfr, Rarray outfi)
{

/** LINEARIZED OPERATOR FROM NEWTON METHOD
    **************************************

    Work in block structure from real and imaginary part to apply the
    linearize operator according to the Newton current solution \phi.
    \phi change according to the Newton iteration,  and is passsed in
    real and imaginary parts separately.

    Input Parameters  : 'inf' the input function in real and imag. parts
                        'phi' the function around the operator was linearized
                        Equation structure and 'lagrange mult.' mu
****************************************************************************/

    int
        i,
        j,
        nx,
        ny;

    double
        b,
        Ome,
        phir2,
        phii2,
        hx,
        hy,
        g,
        Lrr,
        Lri,
        Lii,
        local;

    Rarray
        x,
        y,
        V,
        fr_dx,
        fi_dx,
        fr_dy,
        fi_dy,
        fr_dxdx,
        fi_dxdx,
        fr_dydy,
        fi_dydy;

    b = EQ->b;
    g = EQ->g;
    hx = EQ->hx;
    hy = EQ->hy;
    Ome = EQ->Ome;
    V = EQ->V;
    x = EQ->x;
    y = EQ->y;

    nx = EQ->nx;
    ny = EQ->ny;

    fr_dx = rarrDef(nx*ny);
    fi_dx = rarrDef(nx*ny);
    fr_dy = rarrDef(nx*ny);
    fi_dy = rarrDef(nx*ny);
    fr_dxdx = rarrDef(nx*ny);
    fi_dxdx = rarrDef(nx*ny);
    fr_dydy = rarrDef(nx*ny);
    fi_dydy = rarrDef(nx*ny);

    DfDx_real(nx,ny,infr,hx,fr_dx);
    DfDx_real(nx,ny,infi,hx,fi_dx);

    DfDy_real(nx,ny,infr,hy,fr_dy);
    DfDy_real(nx,ny,infi,hy,fi_dy);

    D2fDx2_real(nx,ny,infr,hx,fr_dxdx);
    D2fDx2_real(nx,ny,infi,hx,fi_dxdx);

    D2fDy2_real(nx,ny,infr,hy,fr_dydy);
    D2fDy2_real(nx,ny,infi,hy,fi_dydy);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            phir2 = phir[i+j*nx]*phir[i+j*nx];
            phii2 = phii[i+j*nx]*phii[i+j*nx];

            local = - mu +  V[i+j*nx] + g * (3*phir2 + phii2);

            Lrr = b * (fr_dxdx[i+j*nx] + fr_dydy[i+j*nx]) + \
                  local * infr[i+j*nx];

            Lri = - Ome * (x[i]*fi_dy[i+j*nx] - y[j]*fi_dx[i+j*nx]) + \
                  2 * g * phir[i+j*nx] * phii[i+j*nx] * infi[i+j*nx];

            outfr[i+j*nx] = Lrr + Lri;

            local = - mu +  V[i+j*nx] + g * (3*phii2 + phir2);

            Lii = b * (fi_dxdx[i+j*nx] + fi_dydy[i+j*nx]) + \
                  local * infi[i+j*nx];

            Lri = Ome * (x[i]*fr_dy[i+j*nx] - y[j]*fr_dx[i+j*nx]) + \
                  2 * g * phir[i+j*nx] * phii[i+j*nx] * infr[i+j*nx];

            outfi[i+j*nx] = Lii + Lri;
        }
    }

    free(fr_dx);
    free(fr_dy);
    free(fi_dx);
    free(fi_dy);
    free(fr_dxdx);
    free(fr_dydy);
    free(fi_dxdx);
    free(fi_dydy);
}





int conjgrad(EqDataPkg EQ, double mu, double tol, Rarray Fr, Rarray Fi,
             Rarray phir, Rarray phii, Rarray fr, Rarray fi)
{

/** CONJUGATE GRADIENT ITERATIVE METHOD FOR 2D-BEC
  * **********************************************
  *
  * Since the Newton method is dereived from a complex equation the system
  * can be splited in two parts, here we chose in real and imaginaty parts
  *
  * The differential operator that is discretized in the given gird assume
  * a 2x2 block form if organized in a matrix. The final vector space size
  * is 2 * nx * ny, where nx and ny are the number of points in each direc
  * tion in the discretized domain.
  *
  * Therefore any time the vector is updated, we call twice  the underlying
  * routine for the first nx*ny elements corresponding to real part and and
  * the later nx*ny related to imaginary part
***************************************************************************/

    int
        N,
        l,
        nx,
        ny,
        maxiter;

    double
        a,
        beta;

    Rarray
        lfr,
        lfi,
        realRes,
        imagRes,
        realDir,
        imagDir,
        prev_fr,
        prev_fi,
        prev_realRes,
        prev_imagRes;

    nx = EQ->nx;
    ny = EQ->ny;

    l = 0;
    N = nx * ny; // Length of the vector is 2 * N
    maxiter = 3 * N;

    lfr = rarrDef(N); // real part after aplly linearized operator
    lfi = rarrDef(N); // imag part after apply linearized operator
    realRes = rarrDef(N);
    imagRes = rarrDef(N);
    realDir = rarrDef(N);
    imagDir = rarrDef(N);
    prev_fr = rarrDef(N);
    prev_fi = rarrDef(N);
    prev_realRes = rarrDef(N);
    prev_imagRes = rarrDef(N);

    linearizedOp(EQ,phir,phii,mu,fr,fi,lfr,lfi);

    // compute residue
    rarrSub(nx*ny,Fr,lfr,realRes);
    rarrSub(nx*ny,Fi,lfi,imagRes);

    // Initialize direction
    rarrCopy(nx*ny,realRes,realDir);
    rarrCopy(nx*ny,imagRes,imagDir);
    
    while (rarrMod(nx*ny,realRes) + rarrMod(nx*ny,imagRes) > tol)
    {
        // lf hold the result of operator acting on direction
        linearizedOp(EQ,phir,phii,mu,realDir,imagDir,lfr,lfi);

        // scalar to update solution and residue
        a = ( rarrDot(N,realRes,realRes) + rarrDot(N,imagRes,imagRes) ) \
            / ( rarrDot(N,realDir,lfr) + rarrDot(N,imagDir,lfi) );

        // record residue from this iteration to safely update
        rarrCopy(N, realRes, prev_realRes);
        rarrCopy(N, imagRes, prev_imagRes);

        // update residue
        rarrUpdate(N, prev_realRes, (-1) * a, lfr, realRes);
        rarrUpdate(N, prev_imagRes, (-1) * a, lfi, imagRes);

        // record solution at this iteration to safely update
        rarrCopy(N, fr, prev_fr);
        rarrCopy(N, fi, prev_fi);

        // update solution
        rarrUpdate(N, prev_fr, a, realDir, fr);
        rarrUpdate(N, prev_fi, a, imagDir, fi);

        // scalar to update direction
        beta = (rarrMod2(N, realRes) + rarrMod2(N, imagRes)) \
               / (rarrMod2(N, prev_realRes) + rarrMod2(N, prev_imagRes));

        rarrUpdate(N,realRes,beta,realDir,realDir);
        rarrUpdate(N,imagRes,beta,imagDir,imagDir);

        l = l + 1; // Update iteration counter

        if (l == maxiter)
        {
            printf("\n\nWARNING : exit before achieve desired residual ");
            printf("value in Conjugate Gradient method due to max number ");
            printf("of iterations given =  %d\n\n", maxiter);
            break;
        }
    }

    // Free function local memory
    free(lfr);
    free(lfi);
    free(realRes);
    free(imagRes);
    free(realDir);
    free(imagDir);
    free(prev_realRes);
    free(prev_imagRes);
    free(prev_fr);
    free(prev_fi);

    return l;
}





int conjgradCond(EqDataPkg EQ, double mu, double tol, Rarray Fr, Rarray Fi,
             Rarray phir, Rarray phii, Rarray fr, Rarray fi)
{

/** PRE-CONDITION VERSION OF CONJUGATE GRADIENT METHOD **/

    int
        N,
        l,
        nx,
        ny,
        maxiter;

    double
        a,
        Lx,
        Ly,
        beta,
        condScalar;

    Rarray
        lfr,
        lfi,
        realRes,
        imagRes,
        realDir,
        imagDir,
        prev_fr,
        prev_fi,
        prev_realRes,
        prev_imagRes;

    Rarray
        realCond,
        imagCond;

    Carray
        fmat;

    nx = EQ->nx;
    ny = EQ->ny;
    Lx = EQ->x[nx-1] - EQ->x[0];
    Ly = EQ->y[ny-1] - EQ->y[0];

    l = 0;
    N = nx * ny; // Length of the vector is 2 * N
    maxiter = 3 * N;

    lfr = rarrDef(N); // real part after aplly linearized operator
    lfi = rarrDef(N); // imag part after apply linearized operator
    realRes = rarrDef(N);
    imagRes = rarrDef(N);
    realDir = rarrDef(N);
    imagDir = rarrDef(N);
    prev_fr = rarrDef(N);
    prev_fi = rarrDef(N);
    prev_realRes = rarrDef(N);
    prev_imagRes = rarrDef(N);

    // New variables needed for pre-conditioning
    realCond = rarrDef(N);
    imagCond = rarrDef(N);
    fmat = carrDef(ny*nx);

    linearizedOp(EQ,phir,phii,mu,fr,fi,lfr,lfi);

    // compute residue
    rarrSub(N,Fr,lfr,realRes);
    rarrSub(N,Fi,lfi,imagRes);

    compressComplex(nx,ny,realRes,imagRes,fmat);
    applyCond(nx,ny,Lx,Ly,mu,EQ->b,fmat);
    extractComplex(nx,ny,realCond,imagCond,fmat);

    // Initialize direction
    extractComplex(nx,ny,realDir,imagDir,fmat);

    // Scalar product with pre-conditioned residue
    condScalar = rarrDot(N,realRes,realCond) + rarrDot(N,imagRes,imagCond);

    while (rarrMod(nx*ny,realRes) + rarrMod(nx*ny,imagRes) > tol)
    {

        // lf hold the result of operator acting on direction
        linearizedOp(EQ,phir,phii,mu,realDir,imagDir,lfr,lfi);

        // scalar to update solution and residue
        a = condScalar / (rarrDot(N,realDir,lfr) + rarrDot(N,imagDir,lfi));

        // record residue from this iteration to safely update
        rarrCopy(N, realRes, prev_realRes);
        rarrCopy(N, imagRes, prev_imagRes);

        // update residue
        rarrUpdate(N, prev_realRes, (-1) * a, lfr, realRes);
        rarrUpdate(N, prev_imagRes, (-1) * a, lfi, imagRes);

        // record solution at this iteration to safely update
        rarrCopy(N, fr, prev_fr);
        rarrCopy(N, fi, prev_fi);

        // update solution
        rarrUpdate(N, prev_fr, a, realDir, fr);
        rarrUpdate(N, prev_fi, a, imagDir, fi);

        // EXTRA STEP FROM PRECONDITIONING acting on new residue
        compressComplex(nx,ny,realRes,imagRes,fmat);
        applyCond(nx,ny,Lx,Ly,mu,EQ->b,fmat);
        extractComplex(nx,ny,realCond,imagCond,fmat);

        // scalar to update direction
        beta = (rarrDot(N,realRes,realCond) + rarrDot(N,imagRes,imagCond)) \
               / condScalar;

        rarrUpdate(N,realCond,beta,realDir,realDir);
        rarrUpdate(N,imagCond,beta,imagDir,imagDir);

        condScalar = beta * condScalar;

        l = l + 1; // Update iteration counter

        if (l == maxiter)
        {
            printf("\n\nWARNING : exit before achieve desired residual ");
            printf("value in Conjugate Gradient method due to max number ");
            printf("of iterations given =  %d\n\n", maxiter);
            break;
        }
    }

    // Free function local memory
    free(lfr);
    free(lfi);
    free(realRes);
    free(imagRes);
    free(realDir);
    free(imagDir);
    free(prev_realRes);
    free(prev_imagRes);
    free(prev_fr);
    free(prev_fi);
    free(realCond);
    free(imagCond);

    free(fmat);

    return l;
}





void stationaryNewton(EqDataPkg EQ, Carray f, double err_tol, int iter_tol)
{

/** MAIN ROUTINE - NEWTON METHOD FOR NONLINEAR OPERATOR
    ***************************************************
    The Newton method is applied in the equation splited in real and imag
    parts, where the resulting linear operator is discretized  what yield
    a linear system which is solved by the iterative method  of conjugate
    gradients.

    Input parameter : 'tol' the maximum allowed error
    Input/Output parameter : function 'f'
    ********************************************************************/

    int
        i,
        nx,
        ny,
        Niter,
        CGiter;

    double
        E,
        b,
        g,
        hx,
        hy,
        mu,
        Ome,
        dev,
        norm,
        error,
        CGtol,
        error_check;

    Rarray
        x,
        y,
        V,
        cgr,
        cgi,
        fr,
        fi,
        abs2,
        Fcg_real,
        Fcg_imag;

    nx = EQ->nx;
    ny = EQ->ny;
    x = EQ->x;
    y = EQ->y;
    V = EQ->V;
    hx = EQ->hx;
    hy = EQ->hy;
    Ome = EQ->Ome;
    g = EQ->g;
    b = EQ->b;

    cgr = rarrDef(nx*ny);
    cgi = rarrDef(nx*ny);
    fr = rarrDef(nx*ny);
    fi = rarrDef(nx*ny);
    abs2 = rarrDef(nx*ny);
    Fcg_real = rarrDef(nx*ny);
    Fcg_imag = rarrDef(nx*ny);

    carrRealPart(nx*ny,f,fr);
    carrImagPart(nx*ny,f,fi);

    mu = creal(Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
    E = creal(Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
    carrAbs2(nx*ny,f,abs2);
    norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));

    stationaryOp(EQ,mu,fr,fi,Fcg_real,Fcg_imag);
    error = maxNorm(nx*ny,Fcg_real,Fcg_imag);
    error_check = error;

    // right hand side of linear operator passed to iteratice CG method
    for (i = 0; i < nx*ny; i++)
    {
        Fcg_real[i] = - Fcg_real[i];
        Fcg_imag[i] = - Fcg_imag[i];
    }

    printf("\nNewton It.      Error    CG It.     ");
    printf("Energy     mu        Norm");
    sepline();



    // NEWTON LOOP
    Niter = 0;
    CGiter = 0;
    while(1)
    {

        // In the Newton loop the norm can be lost.  If the divergence
        // from 1 is larger than 1E-5 correct the norm and start again
        while (error > err_tol)
        {
            printf("%5d       %10.5lf   %6d   ",Niter,error,CGiter);
            printf("%9.5lf  %9.5lf  %9.6lf\n",E,mu,norm);

            // Initial guess for conjugate-gradient method
            rarrFill(nx*ny,0.0,cgr);
            rarrFill(nx*ny,0.0,cgi);

            // Solve with conj. grad. according to current newton error
            if (1E-4 * error > 1E-2) { CGtol = 1E-2; }
            else                     { CGtol = 1E-4 * error; }
            CGiter = conjgrad(EQ,mu,CGtol,Fcg_real,Fcg_imag,fr,fi,cgr,cgi);

            // update solution from newton iteration
            rarrAdd(nx*ny,fr,cgr,fr);
            rarrAdd(nx*ny,fi,cgi,fi);

            for (i = 0; i < nx*ny; i++) f[i] = fr[i] + I * fi[i];
            carrAbs2(nx*ny,f,abs2);
            norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));
            // renormalize(nx,ny,f,hx,hy,1.0);
            // mu = creal(Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
            E = creal(Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
            // carrRealPart(nx*ny,f,fr);
            // carrImagPart(nx*ny,f,fi);

            stationaryOp(EQ,mu,fr,fi,Fcg_real,Fcg_imag);

            error = maxNorm(nx*ny,Fcg_real,Fcg_imag);

            for (i = 0; i < nx*ny; i++)
            {
                Fcg_real[i] = - Fcg_real[i];
                Fcg_imag[i] = - Fcg_imag[i];
            }

            Niter = Niter + 1;

            // check if some progress is being done
            if (Niter % 15 == 0)
            {
                dev = fabs(error - error_check);
                if (dev < err_tol)
                {
                    printf("\nWARNING : Progress Stopped\n");
                    break;
                }
                else error_check = error;
            }

            if (Niter > iter_tol) break;
        }

        // carrAbs2(nx*ny,f,abs2);
        // norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));
        // if (fabs(norm - 1.0) < 1E-4) break;

        printf("%5d       %10.5lf   %6d   ",Niter,error,CGiter);
        printf("%9.5lf  %9.5lf  %9.6lf\n",E,mu,norm);

        renormalize(nx,ny,f,hx,hy,1.0);
        mu = creal(Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
        E = creal(Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
        carrAbs2(nx*ny,f,abs2);
        norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));
        carrRealPart(nx*ny,f,fr);
        carrImagPart(nx*ny,f,fi);

        stationaryOp(EQ,mu,fr,fi,Fcg_real,Fcg_imag);
        error = maxNorm(nx*ny,Fcg_real,Fcg_imag);
        error_check = error;

        // right hand side of linear operator passed to iteratice CG method
        for (i = 0; i < nx*ny; i++)
        {
            Fcg_real[i] = - Fcg_real[i];
            Fcg_imag[i] = - Fcg_imag[i];
        }

        printf("Renormalizing");
        printf("  -------------------------");
        printf("-------------------------\n");
        //printf("\nNewton It.      Error      CG It.       ");
        //printf("Energy       mu        Norm");
        //sepline();

        if (error < err_tol || Niter > iter_tol) break;
    }

    printf("%5d       %10.5lf   %6d   ",Niter,error,CGiter);
    printf("%9.5lf  %9.5lf  %9.6lf",E,mu,norm);
    sepline();

    if (Niter > iter_tol)
    {
        printf("\nWARNING : Achieved maximum iterations allowed");
        printf(" by the user in job-ncg.conf file. Exiting\n");
    }

    free(fr);
    free(fi);
    free(Fcg_real);
    free(Fcg_imag);
    free(cgr);
    free(cgi);
    free(abs2);
}
