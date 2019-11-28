#include "realtimeIntegrator.h"



int RealSplitStepPR(EqDataPkg EQ, int N, double dt, Carray S, int skipFrames,
                    char fname [])
{

    unsigned int
        i,
        j,
        k,
        nx,
        ny,
        countFrames;

    double
        b,
        g,
        Ome,
        hx,
        hy,
        norm;

    double complex
        E,
        Idt,
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

    FILE * f;



    f = fopen(fname, "w");

    if (f == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file '%s'", fname);
        printf(" to write solution in time steps\n\n");
        return -1; 
    }



    Idt = - I * dt; // - I * dt : multiply operators after split-step

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
    norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));



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



    printf("\n\nProgrs     Energy       Norm");
    sepline();



    // Record initial condition
    countFrames = 1;
    fprintf(f,"# Solution at equally spaced time steps\n");
    carr_inline(f,nx*ny,S);

    for (k = 0; k < N; k++)
    {

        // Evolution of potential part in half step
        carrAbs2(nx*ny,S,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,S,linpart);
        carrCopy(nx*ny,linpart,S);



        // Evolution of derivative part entire step
        for (j = 0; j < ny; j++)
        {
            ExplicitY(nx,ny,j,dt,hy,b,Ome,x,S,rhsx);
            triDiag(nx,upperx[j],lowerx[j],midx,rhsx,&linpart[j*nx]);
        }
        carrCopy(nx*ny,linpart,S);

        for (i = 0; i < nx; i++)
        {
            ExplicitX(nx,ny,i,dt,hx,b,Ome,y,S,rhsy);
            triDiag(ny,uppery[i],lowery[i],midy,rhsy,aux);
            for (j = 0; j < ny; j++)
            {
                linpart[i + j*nx] = aux[j];
            }
        }



        // Evolution of potential part another half step
        carrAbs2(nx*ny,linpart,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,linpart,S);



        if ( (k + 1) % 100 == 0 )
        {
            E = Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            carrAbs2(nx*ny,S,abs2);
            norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));

            printf("%5.1lf%%   %11.7lf",(100.0*k)/N,creal(E));
            printf("    %9.7lf",norm);
            printf("\n");
        }

        if ( countFrames == skipFrames + 1 )
        {
            carr_inline(f,nx*ny,S);
            countFrames = 1;
        }
        else
        {
            countFrames = countFrames + 1;
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

    fclose(f);

    return N + 1;
}





int RealSplitStepDYakonov(EqDataPkg EQ, int N, double dt, Carray S,
                          int skipFrames, char fname [])
{

    unsigned int
        i,
        j,
        k,
        nx,
        ny,
        countFrames;

    double
        b,
        g,
        Ome,
        hx,
        hy,
        norm;

    double complex
        E,
        Idt,
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
        auxy,
        auxx;

    Rarray
        x,
        y,
        V,
        pot,
        abs2;

    FILE * f;



    f = fopen(fname, "w");

    if (f == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file '%s'", fname);
        printf(" to write solution in time steps\n\n");
        return -1; 
    }

    Idt = - I * dt; // - I * dt : multiply operators after split-step

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
    auxy = carrDef(ny);
    auxx = carrDef(nx);

    abs2 = rarrDef(nx*ny);

    // compute norm of initial guess to maintain in the propagation
    carrAbs2(nx*ny,S,abs2);
    norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));



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



    printf("\n\nProgrs     Energy       norm");
    sepline();



    // Record initial condition
    countFrames = 1;
    fprintf(f,"# Solution at equally spaced time steps\n");
    carr_inline(f,nx*ny,S);

    for (k = 0; k < N; k++)
    {

        // Evolution of potential part in half step
        carrAbs2(nx*ny,S,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,S,linpart);
        carrCopy(nx*ny,linpart,S);



        // Evolution of derivative part entire step
        for (j = 0; j < ny; j++)
        {
            ExplicitY(nx,ny,j,dt,hy,b,Ome,x,S,auxx);
            ExplicitX_alongX(nx,ny,dt,hx,b,Ome,y[j],auxx,rhsx);
            triDiag(nx,upperx[j],lowerx[j],midx,rhsx,&linpart[j*nx]);
        }
        carrCopy(nx*ny,linpart,S);

        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                rhsy[j] = S[i + j*nx];
            }
            triDiag(ny,uppery[i],lowery[i],midy,rhsy,auxy);
            for (j = 0; j < ny; j++)
            {
                linpart[i + j*nx] = auxy[j];
            }
        }



        // Evolution of potential part another half step
        carrAbs2(nx*ny,linpart,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,linpart,S);



        if ( (k + 1) % 100 == 0 )
        {
            E = Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            carrAbs2(nx*ny,S,abs2);
            norm = sqrt(Rsimps2D(nx,ny,abs2,hx,hy));

            printf("%5.1lf%%   %11.7lf",(100.0*k)/N,creal(E));
            printf("    %9.7lf",norm);
            printf("\n");
        }

        if ( countFrames == skipFrames + 1 )
        {
            carr_inline(f,nx*ny,S);
            countFrames = 1;
        }
        else
        {
            countFrames = countFrames + 1;
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
    free(auxx);
    free(auxy);

    fclose(f);

    return N + 1;
}
