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
                     0.5 * (I*x[i]*sig + mu) * in[i + (j+1)*nx];
        }

        return;
    }

    if (j == ny - 1)
    {
        for (i = 0; i < nx; i++)
        {
            out[i] = (I - mu) * in[i + j*nx] + \
                     0.5 * (mu - I*x[i]*sig) * in[i + (j-1)*nx];
        }

        return;
    }

    for (i = 0; i < nx; i++)
    {
        out[i] = (I - mu) * in[i + j*nx] + \
                 0.5 * (I*x[i]*sig + mu) * in[i + (j+1)*nx] + \
                 0.5 * (mu - I*x[i]*sig) * in[i + (j-1)*nx];
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



void ExplicitX_alongX(int nx, int ny, doublec dt, double hx, double b,
     double Ome, double yj, Carray in, Carray out)
{
    int
        i;

    double complex
        mu,
        sig;

    mu = b * dt / hx / hx;
    sig = 0.5 * Ome * dt / hx;

    out[0] = (I - mu) * in[0] +  0.5 * (-I*yj*sig + mu) * in[1];

    out[nx-1] = (I - mu) * in[nx - 1] + 0.5 * (mu + I*yj*sig) * in[nx - 2];

    for (i = 1; i < nx - 1; i++)
    {
        out[i] = (I - mu) * in[i] + \
                 0.5 * (-I*yj*sig + mu) * in[i+1] + \
                 0.5 * (mu + I*yj*sig) * in[i-1];
    }
}





int SplitStepPR(EqDataPkg EQ, int N, double realDT, Carray S)
{

    unsigned int
        i,
        j,
        k,
        nx,
        ny,
	tid,
	nthreads;

    double
        b,
        g,
        Ome,
        hx,
        hy,
        Idt,
        norm,
        vir,
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
        auxy,
	Ux,
	Lx,
	Uy,
	Ly;


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
    Ux = carrDef(nx*ny);
    Lx = carrDef(nx*ny);
    Uy = carrDef(nx*ny);
    Ly = carrDef(nx*ny);

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
	LU(nx,upperx[j],lowerx[j],midx,&Lx[j*nx],&Ux[j*nx]);
    }

    // Implicit on y-direction
    midy = I + b*dt/hy/hy;
    for (i = 0; i < nx; i++)
    {
        uppery[i] = -0.5*b*dt/hy/hy - I * 0.25*x[i]*Ome*dt/hy;
        lowery[i] = -0.5*b*dt/hy/hy + I * 0.25*x[i]*Ome*dt/hy;
	LU(ny,uppery[i],lowery[i],midy,&Ly[i*ny],&Uy[i*ny]);
    }



    printf("\n\nProgrs     Energy       mu");
    printf("           <r>          Max. Res.");
    printf("     Virial");
    sepline();



    for (k = 0; k < N; k++)
    {

        // Evolution of potential part in half step
        carrAbs2(nx*ny,S,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,S,linpart);
        carrCopy(nx*ny,linpart,S);


#pragma omp parallel private(j,i,rhsx,rhsy,auxy,tid,nthreads)
	{

	rhsx = carrDef(nx);
	rhsy = carrDef(ny);
	auxy = carrDef(ny);

        tid = omp_get_thread_num();        // thread Id
        nthreads = omp_get_num_threads();  // number of threads being used

        // Evolution of derivative part entire step
        for (j = tid; j < ny; j += nthreads)
        {
            ExplicitY(nx,ny,j,dt,hy,b,Ome,x,S,rhsx);
            triDiagLU(nx,&Lx[j*nx],&Ux[j*nx],upperx[j],rhsx,&linpart[j*nx]);
        }

	// wait for all threads to solve for y-direction implicitly

#pragma omp barrier

        for (i = tid; i < nx; i += nthreads)
        {
            ExplicitX(nx,ny,i,dt,hx,b,Ome,y,linpart,rhsy);
            triDiagLU(ny,&Ly[i*ny],&Uy[i*ny],uppery[i],rhsy,auxy);
            for (j = 0; j < ny; j++) S[i + j*nx] = auxy[j];
        }

	free(rhsy);
	free(rhsx);
	free(auxy);

	}

        carrCopy(nx*ny,S,linpart);



        // Evolution of potential part another half step
        carrAbs2(nx*ny,linpart,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,linpart,S);



        // Renormalize
        renormalize(nx,ny,S,hx,hy,norm);


        if ( (k + 1) % 100 == 0 )
        {
            E = Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            mu = Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            meanr = MeanR(nx,ny,S,hx,hy,x,y);
            maxres = MaxResidue(nx,ny,hx,hy,b,Ome,g,V,x,y,S,mu);
            avgres = AvgResidue(nx,ny,hx,hy,b,Ome,g,V,x,y,S,mu);
            vir = Virial(nx,ny,hx,hy,b,g,V,S);

            printf("%5.1lf%%   %11.7lf",(100.0*k)/N,creal(E));
            printf("  %11.7lf",creal(mu));
            printf("    %9.7lf",meanr);
            printf("    %8.5lf",maxres);
            printf("    %9.5lf",vir);
            printf("\n");
        }

    }

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upperx);
    free(uppery);
    free(lowerx);
    free(lowery);
    free(pot);
    free(Ux);
    free(Lx);
    free(Uy);
    free(Ly);

    return N + 1;
}





int SplitStepDYakonov(EqDataPkg EQ, int N, double realDT, Carray S)
{

    unsigned int
        i,
        j,
        k,
        nx,
        ny,
	tid,
	nthreads;

    double
        b,
        g,
        Ome,
        hx,
        hy,
        Idt,
        norm,
        vir,
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
        auxy,
        auxx,
	Ux,
	Lx,
	Uy,
	Ly;

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
    Ux = carrDef(nx*ny);
    Lx = carrDef(nx*ny);
    Uy = carrDef(nx*ny);
    Ly = carrDef(nx*ny);

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
	LU(nx,upperx[j],lowerx[j],midx,&Lx[j*nx],&Ux[j*nx]);
    }

    // Implicit on y-direction
    midy = I + b*dt/hy/hy;
    for (i = 0; i < nx; i++)
    {
        uppery[i] = -0.5*b*dt/hy/hy - I * 0.25*x[i]*Ome*dt/hy;
        lowery[i] = -0.5*b*dt/hy/hy + I * 0.25*x[i]*Ome*dt/hy;
	LU(ny,uppery[i],lowery[i],midy,&Ly[i*ny],&Uy[i*ny]);
    }



    printf("\n\nProgrs     Energy       mu");
    printf("           <r>          Max. Res.");
    printf("     Virial");
    sepline();



    for (k = 0; k < N; k++)
    {

        // Evolution of potential part in half step
        carrAbs2(nx*ny,S,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,S,linpart);
        carrCopy(nx*ny,linpart,S);


	// EVOLUTION OF LINEAR PART USING ADI-D'YAKONOV METHOD

#pragma omp parallel private(j,i,rhsx,rhsy,auxy,auxx,tid,nthreads)
	{

	rhsx = carrDef(nx);
	rhsy = carrDef(ny);
	auxy = carrDef(ny);
	auxx = carrDef(nx);

        tid = omp_get_thread_num();        // thread Id
        nthreads = omp_get_num_threads();  // number of threads being used

        for (j = tid; j < ny; j += nthreads)
        {
            ExplicitY(nx,ny,j,dt,hy,b,Ome,x,S,auxx);
            ExplicitX_alongX(nx,ny,dt,hx,b,Ome,y[j],auxx,rhsx);
            triDiagLU(nx,&Lx[j*nx],&Ux[j*nx],upperx[j],rhsx,&linpart[j*nx]);
        }

	// wait for all threads to solve for y-direction implicitly

#pragma omp barrier

        for (i = tid; i < nx; i += nthreads)
        {
            for (j = 0; j < ny; j++) rhsy[j] = linpart[i + j*nx];
            triDiagLU(ny,&Ly[i*ny],&Uy[i*ny],uppery[i],rhsy,auxy);
            for (j = 0; j < ny; j++) S[i + j*nx] = auxy[j];
        }

	free(rhsy);
	free(rhsx);
	free(auxy);
	free(auxx);

	}

        carrCopy(nx*ny,S,linpart);



        // Evolution of potential part another half step
        carrAbs2(nx*ny,linpart,abs2);
        rarrUpdate(nx*ny,V,g,abs2,pot);
        rcarrExp(nx*ny,0.5*Idt,pot,stepexp);
        carrMultiply(nx*ny,stepexp,linpart,S);



        // Renormalize
        renormalize(nx,ny,S,hx,hy,norm);


        if ( (k + 1) % 100 == 0 )
        {
            E = Energy(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            mu = Chem(nx,ny,hx,hy,b,Ome,g,V,x,y,S);
            meanr = MeanR(nx,ny,S,hx,hy,x,y);
            maxres = MaxResidue(nx,ny,hx,hy,b,Ome,g,V,x,y,S,mu);
            avgres = AvgResidue(nx,ny,hx,hy,b,Ome,g,V,x,y,S,mu);
            vir = Virial(nx,ny,hx,hy,b,g,V,S);

            printf("%5.1lf%%   %11.7lf",(100.0*k)/N,creal(E));
            printf("  %11.7lf",creal(mu));
            printf("    %9.7lf",meanr);
            printf("    %8.5lf",maxres);
            printf("    %9.5lf",vir);
            printf("\n");
        }

    }

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upperx);
    free(uppery);
    free(lowerx);
    free(lowery);
    free(pot);
    free(Ux);
    free(Lx);
    free(Uy);
    free(Ly);

    return N + 1;
}
