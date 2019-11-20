#include "imagtimeIntegrator.h"



void ExplicitY(int nx, int ny, int j, doublec dt, double hy, double b,
     double Ome, Rarray x, Carray in, Carray out)
{
    int
        i;

    double complex
        mu,
        sig;

    mu = b * dt / hy / hy;
    sig = 0.5 * Ome * dt / hy;

    if (j == 0)
    {
        for (i = 0; i < nx; i++)
        {
            out[i] = (I - mu) * in[i + j*nx] + \
                     (0.5*I*x[i]*sig + 0.5*mu) * in[i + (j+1)*nx];
        }

        return;
    }

    if (j == ny - 1)
    {
        for (i = 0; i < nx; i++)
        {
            out[i] = (I - mu) * in[i + j*nx] + \
                     (0.5*mu - 0.5*I*x[i]*sig) * in[i + (j-1)*nx];
        }

        return;
    }

    for (i = 0; i < nx; i++)
    {
        out[i] = (I - mu) * in[i + j*nx] + \
                 (0.5*I*x[i]*sig + 0.5*mu) * in[i + (j+1)*nx] + \
                 (0.5*mu - 0.5*I*x[i]*sig) * in[i + (j-1)*nx];
    }
}





void ExplicitX(int nx, int ny, int i, doublec dt, double hx, double b,
     double Ome, Rarray y, Carray in, Carray out)
{
    int
        j;

    double complex
        mu,
        sig;

    mu = b * dt / hx / hx;
    sig = 0.5 * Ome * dt / hx;

    if (i == 0)
    {
        for (j = 0; j < ny; j++)
        {
            out[j] = (I - mu) * in[i + j*nx] + \
                     0.5 * (-I*y[j]*sig + mu) * in[i+1 + j*nx];
        }

        return;
    }

    if (i == nx - 1)
    {
        for (j = 0; j < ny; j++)
        {
            out[j] = (I - mu) * in[i + j*nx] + \
                     0.5 * (mu + I*y[j]*sig) * in[i-1 + j*nx];
        }

        return;
    }

    for (j = 0; j < ny; j++)
    {
        out[j] = (I - mu) * in[i + j*nx] + \
                 0.5 * (-I*y[j]*sig + mu) * in[i+1 + j*nx] + \
                 0.5 * (mu + I*y[j]*sig) * in[i-1 + j*nx];
    }
}





int SplitStepPR(EqDataPkg EQ, int N, double realDT, Carray S)
{

    unsigned int
        i,
        j,
        k,
        nx,
        ny;

    double
        b,
        g,
        Ome,
        hx,
        hy,
        Idt,
        norm,
        meanr,
        maxres,
        avgres;

    double complex
        E,
        mu,
        dt,
        midx,
        midy;

    Carray
        stepexp,
        linpart,
        upperx,
        uppery,
        lowerx,
        lowery,
        rhsx,
        rhsy,
        aux;

    Rarray
        x,
        y,
        V,
        pot,
        abs2;



    dt  = - I * realDT; // Wick rotation - attenuate e^(- i E T)
    Idt = - realDT;     // - I * dt : multiply operators after split-step

    nx = EQ->nx;    // discretize points in x-direction
    ny = EQ->ny;    // discretize points in y-direction
    hx = EQ->hx;    // grid spacing in x-direction
    hy = EQ->hy;    // grid spacing in y-direction
    b = EQ->b;      // second order derivative coefficient
    Ome = EQ->Ome;  // Rotation frequency
    g = EQ->g;      // contact interaction strength
    V = EQ->V;      // Potential in grid points
    x = EQ->x;      // x-direction grid points
    y = EQ->y;      // y-direction grid points

    pot = rarrDef(nx*ny); // potential - trap and interaction combined
    stepexp = carrDef(nx*ny); // Result after applied potential part
    linpart = carrDef(nx*ny); // result after solved differential part

    // diagonals of tridiagonal systems. When implicit in x-direction
    // there are ny tridiagonal systems to be solved, each one with
    // constant entries, but not equal from one system to another,
    // thus there are ny constant values to setup. Analogously
    // when y-direction is implicit
    upperx = carrDef(ny);
    uppery = carrDef(nx);
    lowerx = carrDef(ny);
    lowery = carrDef(nx);

    // Right hand size for each implicit direction
    rhsx = carrDef(nx);
    rhsy = carrDef(ny);
    aux = carrDef(ny);

    abs2 = rarrDef(nx*ny);

    // compute norm of initial guess to maintain in the propagation
    carrAbs2(nx*ny,S,abs2);
    // norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));
    norm = 1.0;



    // Configure the implicit part of Peaceman-Rachford ADI method, for both
    // x- and y-direction. In the presence of rotation, the step implicit in
    // x-direction depends on the value for y and analogously when y is  the
    // implicit direction, it depends on x value. Nonetheless for all system
    // of equations the entries are constant for each fixed value of the
    // explicit direction



    // Implicit on x-direction
    midx = I + b*dt/hx/hx;
    for (j = 0; j < ny; j++)
    {
        upperx[j] = -0.5*b*dt/hx/hx + I * 0.25*y[j]*Ome*dt/hx;
        lowerx[j] = -0.5*b*dt/hx/hx - I * 0.25*y[j]*Ome*dt/hx;
    }

    // Implicit on y-direction
    midy = I + b*dt/hy/hy;
    for (i = 0; i < nx; i++)
    {
        uppery[i] = -0.5*b*dt/hy/hy - I * 0.25*x[i]*Ome*dt/hy;
        lowery[i] = -0.5*b*dt/hy/hy + I * 0.25*x[i]*Ome*dt/hy;
    }



    printf("\n\n    time     Energy           mu");
    printf("               <r>          Max. Res.");
    printf("    Avg. Res.");
    sepline();



    for (k = 0; k < N; k++)
    {

        carrAbs2(nx*ny,S,abs2);

        // Evolution of potential part in half step
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,S,linpart);



        // Evolution of derivative part entire step
        for (j = 0; j < ny; j++)
        {
            ExplicitY(nx,ny,j,dt,hy,b,Ome,x,linpart,rhsx);
            triDiag(nx,upperx[j],lowerx[j],midx,rhsx,&linpart[j*nx]);
        }

        for (i = 0; i < nx; i++)
        {
            ExplicitX(nx,ny,i,dt,hx,b,Ome,y,linpart,rhsy);
            triDiag(ny,uppery[i],lowery[i],midy,rhsy,aux);
            for (j = 0; j < ny; j++)
            {
                linpart[i + j*nx] = aux[j];
            }
        }



        carrAbs2(nx*ny,linpart,abs2);

        // Evolution of potential part another half step
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,linpart,S);



        // Renormalize
        renormalize(nx,ny,S,hx,hy,norm);


        if ( (k+1) % 100 == 0 )
        {
            E = Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            mu = Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            meanr = MeanR(nx,ny,S,hx,hy,x,y);
            maxres = MaxResidue(nx,ny,hx,hy,b,Ome,g,V,x,y,S,mu);
            avgres = AvgResidue(nx,ny,hx,hy,b,Ome,g,V,x,y,S,mu);

            printf("\n%8.3lf   %15.7E",(k + 1)*realDT,creal(E));
            printf("  %15.7E",creal(mu));
            printf("    %9.7lf",meanr);
            printf("    %9.7lf",maxres);
            printf("    %9.7lf",avgres);
        }

    }

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upperx);
    free(uppery);
    free(lowerx);
    free(lowery);
    free(rhsx);
    free(rhsy);
    free(pot);
    free(aux);

    return N + 1;
}
