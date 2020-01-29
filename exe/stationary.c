#include <string.h>
#include <stdio.h>
#include <math.h>
#include "inout.h"
#include "newtoncg.h"
#include "linearPotential.h"



void TimePrint(double t)
{
    
    // format and print time in days / hours / minutes

    int
        tt = (int) t,
        hours = 0,
        mins  = 0;

    if ( tt / 3600  > 0 )
    { hours = tt / 3600;  tt = tt % 3600;  }

    if ( tt / 60    > 0 )
    { mins  = tt / 60;    tt = tt % 60;    }

    printf(" %d hour(s) %d minute(s)",hours,mins);
}



void ReachNewLine(FILE * f)
{

/** Read until get new line or 'end of file' in a opened file.
    Function to skip comment lines in job configuration file **/

    char
        sentinel;

    while (1)
    {
        fscanf(f, "%c", &sentinel);
        if (sentinel == '\n' || sentinel == EOF) return;
    }
}





void SaveConf(FILE * f, EqDataPkg EQ, double dt, int N)
{

    double
        xi,
        xf,
        yi,
        yf;

    xi = EQ->x[0];
    yi = EQ->y[0];
    xf = EQ->x[EQ->nx-1];
    yf = EQ->y[EQ->ny-1];

    // Grid domain
    fprintf(f, "%.10lf %.10lf %d ", xi, xf, EQ->nx);
    fprintf(f, "%.10lf %.10lf %d ", yi, yf, EQ->ny);

    // Equation parameters
    fprintf(f, "%.15lf %.15lf %.15lf ", EQ->b, EQ->Ome, EQ->g);
    fprintf(f, "%.15lf %.15lf %.15lf %.15lf\n",
            EQ->p[0], EQ->p[1], EQ->p[2], EQ->p[3]);

}





EqDataPkg SetupParams(FILE * paramFile, FILE * confFile,
          char Vname[], double * dt, int * N)
{

/** Read a line of _domain file and _eq to setup initial function
    in the grid, the trap potential and the equation parameters **/

    int
        k,
        nx,
        ny;

    double
        b,
        g,
        Ome,
        xi,
        xf,
        yi,
        yf,
        p[4];



    // Domain grid information
    k = fscanf(confFile, "%lf %lf %d ", &xi, &xf, &nx);
    k = fscanf(confFile, "%lf %lf %d ", &yi, &yf, &ny);
    // ISSUE - time step and Number of steps are useless  for
    // newton method but they are kept as notation convention
    k = fscanf(confFile, "%lf %d", dt, N);

    // Read a line of numbers corresponding to equation parameters
    k = fscanf(paramFile,"%lf %lf %lf %lf %lf %lf %lf", &b, &Ome, &g,
        &p[0], &p[1], &p[2], &p[3]);

    return PackEqData(nx,ny,xi,xf,yi,yf,b,g,Ome,Vname,p);

}





int main(int argc, char * argv[])
{

    // DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER ARCHITECTURE
    omp_set_num_threads(omp_get_max_threads() / 2);



    int
        N,
        i,
        j,
        k,
        nx, // number of grid points along x
        ny, // number of grid points along y
        iter_tol; // Maximum allowed iterations in newton method



    double
        err_tol,    // Maximum error allowed to stop newton iterations
        start,      // start trigger to measure time
        time_used,  // Time used in calling evolution routine
        dt,
        real,   // real part of initial guess read from file
        imag;   // imag part of initial guess read from file



    char
        c,
        fname[100],     // name of files to open
        potname[50],    // name of trap function in linearPotential.c
        infname[100],   // input file name prefix
        outfname[100];  // output file name prefix



    FILE
        * domain_file,
        * job_file,
        * f0_file,
        * eq_file;



    Rarray
        abs2;



    Carray
        S;



    EqDataPkg
        EQ;





    job_file = fopen("job-ncg.conf","r");
    if (job_file == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n","job-ncg.conf");
        exit(EXIT_FAILURE);
    }



    // Read fields of configuration file to setup the job to do
    i = 1;
    while ( (c  = getc(job_file)) != EOF)
    {

        // jump comment line initiated with #
        if (c == '#') { ReachNewLine(job_file); continue; }
        else          { fseek(job_file, -1, SEEK_CUR);    }

        switch (i)
        {
            case 1:
                fscanf(job_file, "%s", potname);
                i = i + 1;
                break;
            case 2:
                fscanf(job_file, "%s", infname);
                i = i + 1;
                break;
            case 3:
                fscanf(job_file, "%s", outfname);
                i = i + 1;
                break;
            case 4:
                fscanf(job_file, "%lf", &err_tol);
                i = i + 1;
                break;
            case 5:
                fscanf(job_file, "%d", &iter_tol);
                i = i + 1;
                break;
        }

        ReachNewLine(job_file);

    }

    fclose(job_file);










    /*********************************************************************
     *********************************************************************
     ***************     OPEN TO EXECUTE LIST OF JOBS      ***************
     *********************************************************************
     *********************************************************************/



    // open file with values of equation's parameters
    printf("\n\n");
    printf("\t\t*********************************************\n");
    printf("\t\t*                                           *\n");
    printf("\t\t*           SEARCHING SETUP FILES           *\n");
    printf("\t\t*                                           *\n");
    printf("\t\t*********************************************\n");

    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_eq.dat");
    printf("\nLooking for %s", fname);
    eq_file = fopen(fname, "r");
    if (eq_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ....... Found !");
    }



    // open file to configure grid domain
    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_domain.dat");
    printf("\nLooking for %s", fname);
    domain_file = fopen(fname, "r");
    if (domain_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ... Found !");
    }



    // open file with values of initial function in grid points
    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_init.dat");
    printf("\nLooking for %s", fname);
    f0_file = fopen(fname, "r");
    if (f0_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ..... Found !");
    }



    // open file to write parameters of domain and equation
    strcpy(fname,"output/");
    strcat(fname, outfname);
    strcat(fname,"_conf_ncg.dat");
    job_file = fopen(fname, "w");
    if (job_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    fprintf(job_file, "# Trap name : %s\n", potname);



    // Setup structure to hold all necessary input data that
    // define the equation to be solved
    EQ = SetupParams(eq_file, domain_file, potname, &dt, &N);
    nx = EQ->nx;
    ny = EQ->ny;



    // SETUP INITIAL FUNCTION IN DISCRETIZED POINTS
    S = carrDef(nx*ny);
    abs2 = rarrDef(nx*ny);
    for (i = 0; i < nx*ny; i++)
    {
        k = fscanf(f0_file, " (%lf%lfj)", &real, &imag);
        S[i] = real + I * imag;
        abs2[i] = real*real + imag*imag;
    }
    fclose(f0_file);

    printf("\nGot Initial condition with || . || =");
    printf(" %.6lf\n", sqrt(Rsimps2D(nx,ny,abs2,EQ->hx,EQ->hy)));










    printf("\n\n\n\n");
    printf("\t\t*********************************************\n");
    printf("\t\t*                                           *\n");
    printf("\t\t*            GRID SPECIFICATIONS            *\n");
    printf("\t\t*                                           *\n");
    printf("\t\t*********************************************\n");
    printf("\n");

    printf("x = [ %.2lf , %.2lf , ... , %.2lf , %.2lf ]\n",
	        EQ->x[0],EQ->x[1],EQ->x[nx-2],EQ->x[nx-1]);
    printf("%d points for x-direction | x-grid-spacing = %.3lf\n",nx,EQ->hx);

    printf("y = [ %.2lf , %.2lf , ... , %.2lf , %.2lf ]\n",
	        EQ->y[0],EQ->y[1],EQ->y[ny-2],EQ->y[ny-1]);
    printf("%d points for y-direction | y-grid-spacing = %.3lf\n",ny,EQ->hy);










    /*********************************************************************
     *********************************************************************
     ********************    CALL INTEGRATION ROUTINE    *****************
     *********************************************************************
     *********************************************************************/

    printf("\n\nStarting the iterations");
    printf(" : error_tol = %.1E | iter_tol = %d\n\n",err_tol,iter_tol);

    start = omp_get_wtime();

    stationaryNewton(EQ,S,err_tol,iter_tol);

    time_used = (double) (omp_get_wtime() - start);

    printf("\n\nTime taken for Newton-CG : %.0lf sec = ",time_used);
    TimePrint(time_used);

    // Record data
    strcpy(fname, "output/");
    strcat(fname, outfname);
    strcat(fname, "_ncg.dat");

    carr_txt(fname,nx*ny,S);
    SaveConf(job_file,EQ,dt,N);





    fclose(job_file);
    fclose(eq_file);
    fclose(domain_file);
    free(S);
    free(abs2);
    ReleaseEqDataPkg(EQ);





    printf("\n\nDone\n\n");
    return 0;
}
