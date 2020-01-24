#include "newtoncg.h"



double maxNorm(int N, Rarray fr, Rarray fi)
{
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



void stationaryOp(EqDataPkg EQ, double mu, Rarray fr, Rarray fi, Rarray Fr,
                  Rarray Fi)
{

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

    DfDx_real(nx,ny,fr_dx,hx,fr_dxdx);
    DfDx_real(nx,ny,fi_dx,hx,fi_dxdx);

    DfDy_real(nx,ny,fr_dy,hy,fr_dydy);
    DfDy_real(nx,ny,fi_dy,hy,fi_dydy);

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



    INPUT AND OUTPUT PARAMETERS
    ***************************

    inf(r/i) are the input real and imag part of the function
    outf(r/i) are the result of the action of the operator in inf(r/i) **/


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

    DfDx_real(nx,ny,fr_dx,hx,fr_dxdx);
    DfDx_real(nx,ny,fi_dx,hx,fi_dxdx);

    DfDy_real(nx,ny,fr_dy,hy,fr_dydy);
    DfDy_real(nx,ny,fi_dy,hy,fi_dydy);

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
  * tion in the discretized domain. Although here we work separetely on the
  * real and imaginary parts to track the derivation steps.
  *
  * Therefore any time the vector is update, we call twice the underlying
  * routine for the first nx*ny elements corresponding to real part and
  * and the later nx*ny related to imaginary part
  *
  * INPUT/OUTPUT PARAMETERS
  * ***********************
  *
  * Fr and Fi are the real and imaginary part of the righ-hand-side of the
  * linear system to be solved
  *
  * fr and fi are the initial guess (good choice is zero), but finish with
  * the corresponding solution if the method success **/

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



void stationaryNewton(EqDataPkg EQ, Carray f, double tol)
{

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
        error,
        CGtol;

    Rarray
        x,
        y,
        V,
        cgr,
        cgi,
        fr,
        fi,
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
    Fcg_real = rarrDef(nx*ny);
    Fcg_imag = rarrDef(nx*ny);

    carrRealPart(nx*ny,f,fr);
    carrImagPart(nx*ny,f,fi);

    mu = creal(Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
    E = creal(Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,f));

    stationaryOp(EQ,mu,fr,fi,Fcg_real,Fcg_imag);

    error = maxNorm(nx*ny,Fcg_real,Fcg_imag);

    // right hand side of linear operator passed to iteratice CG method
    for (i = 0; i < nx*ny; i++)
    {
        Fcg_real[i] = - Fcg_real[i];
        Fcg_imag[i] = - Fcg_imag[i];
    }

    printf("\nNewton It.      Error      CG It.       ");
    printf("Energy       mu");



    // NEWTON LOOP
    Niter = 0;
    CGiter = 0;
    sepline();
    while (error > tol)
    {
        printf("%5d       %10.5lf     %6d     ",Niter,error,CGiter);
        printf("%9.5lf    %9.5lf\n",E,mu);

        // Initial guess for conjugate-gradient method
        rarrFill(nx*ny,0.0,cgr);
        rarrFill(nx*ny,0.0,cgi);

        // Solve with conj. grad. according to current newton error
        if (1E-5 * error > 1E-2) { CGtol = 1E-2; }
        else                     { CGtol = 1E-5 * error; }
        CGiter = conjgrad(EQ,mu,CGtol,Fcg_real,Fcg_imag,fr,fi,cgr,cgi);

        // update solution from newton iteration
        rarrAdd(nx*ny,fr,cgr,fr);
        rarrAdd(nx*ny,fi,cgi,fi);

        for (i = 0; i < nx*ny; i++) f[i] = fr[i] + I * fi[i];
        renormalize(nx,ny,f,hx,hy,1.0);
        mu = creal(Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
        E = creal(Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,f));
        carrRealPart(nx*ny,f,fr);
        carrImagPart(nx*ny,f,fi);

        stationaryOp(EQ,mu,fr,fi,Fcg_real,Fcg_imag);

        error = maxNorm(nx*ny,Fcg_real,Fcg_imag);

        for (i = 0; i < nx*ny; i++)
        {
            Fcg_real[i] = - Fcg_real[i];
            Fcg_imag[i] = - Fcg_imag[i];
        }

        Niter = Niter + 1;
    }

    printf("%5d       %10.5lf     %6d     ",Niter,error,CGiter);
    printf("%9.5lf    %9.5lf",E,mu);
    sepline();

    free(fr);
    free(fi);
    free(Fcg_real);
    free(Fcg_imag);
    free(cgr);
    free(cgi);
}
