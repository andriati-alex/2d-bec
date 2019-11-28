#include <string.h>
#include <stdio.h>
#include <math.h>
#include "imagtimeIntegrator.h"
#include "linearPotential.h"





void ReachNewLine(FILE * f)
{

    // Read until get new line or 'end of file' in a opened file.

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

    // Time domain
    fprintf(f, "%.15lf %d ", dt, N);

    // Equation parameters
    fprintf(f, "%.15lf %.15lf %.15lf ", EQ->b, EQ->Ome, EQ->g);
    fprintf(f, "%.15lf %.15lf %.15lf\n", EQ->p[0], EQ->p[1], EQ->p[2]);

}





EqDataPkg SetupParams(FILE * paramFile, FILE * confFile,
          char Vname[], double * dt, int * N)
{

/** Read line by line of _domain file and _eq to setup a new integration **/

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
        p[3];



    // Read spatial and time domain settings
    // -------------------------------------

    k = fscanf(confFile, "%lf %lf %d ", &xi, &xf, &nx);
    k = fscanf(confFile, "%lf %lf %d ", &yi, &yf, &ny);
    k = fscanf(confFile, "%lf %d", dt, N);

    // Read a line of numbers corresponding to equation parameters
    // -----------------------------------------------------------

    k = fscanf(paramFile,"%lf %lf %lf %lf %lf %lf", &b, &Ome, &g,
        &p[0], &p[1], &p[2]);

    return PackEqData(nx,ny,xi,xf,yi,yf,b,g,Ome,Vname,p);

}





int main(int argc, char * argv[])
{

    /*  DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER ARCHITECTURE
     *  --------------------------------------------------------------- */

    omp_set_num_threads(omp_get_max_threads() / 2);



    /* ==================================================================== *
     *
     *                         VARIABLES DEFINITION
     *
     * ==================================================================== */



    int
        N,
        i,
        j,
        k,
        nx,
        ny,
        method,
        Nlines,
        resetinit;



    double
        start,      // start trigger to measure time
        time_used,  // Time used in calling evolution routine
        dt,         // time step
        real,   // real part of initial guess read from file
        imag;   // imag part of initial guess read from file



    char
        c,
        timeinfo,       // 'r' for real and 'i' for imaginary
        fname[100],     // name of files to open
        potname[50],    // name of trap function in linearPotential.c
        infname[100],   // input file name prefix
        outfname[100];  // output file name prefix



    FILE
        * domain_file,
        * job_file,
        * orb_file,
        * eq_file;



    Rarray
        abs2;



    Carray
        S;



    EqDataPkg
        EQ;





    job_file = fopen("job.conf", "r");

    if (job_file == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", "job.conf");
        exit(EXIT_FAILURE);
    }

    i = 1;

    while ( (c  = getc(job_file)) != EOF)
    {

        // jump comment line marked as #
        if (c == '#') { ReachNewLine(job_file); continue; }
        else          { fseek(job_file, -1, SEEK_CUR);    }

        switch (i)
        {
            case 1:
                fscanf(job_file, "%s", fname);
                timeinfo = fname[0];
                i = i + 1;
                break;
            case 2:
                fscanf(job_file, "%s", potname);
                i = i + 1;
                break;
            case 3:
                fscanf(job_file, "%s", infname);
                i = i + 1;
                break;
            case 4:
                fscanf(job_file, "%s", outfname);
                i = i + 1;
                break;
            case 5:
                fscanf(job_file, "%d", &method);
                i = i + 1;
                break;
            case 6:
                fscanf(job_file, "%d", &Nlines);
                i = i + 1;
                break;
            case 7:
                fscanf(job_file, "%d", &resetinit);
                i = i + 1;
                break;
        }

        ReachNewLine(job_file);

    }

    fclose(job_file);





    if (timeinfo != 'r' && timeinfo != 'R')
    {
        if (timeinfo != 'i' && timeinfo != 'I')
        {
            printf("\n\nInvalid keyword for propagation type.\n");
            printf("Valid ones are either 'imag' or 'real'.\n\n");
            exit(EXIT_FAILURE);
        }
    }










    /*  ===============================================================
     
                   LET FILES OPENNED TO EXECUTE LIST OF JOBS
     
        ===============================================================  */



    // open file with values of equation's parameters

    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_eq.dat");

    printf("\nLooking for %s", fname);

    eq_file = fopen(fname, "r");

    if (eq_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        return -1;
    }
    else
    {
        printf(" ....... Found !\n");
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
        return -1;
    }
    else
    {
        printf(" ... Found !\n");
    }



    // open file with values of initial condition

    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_init.dat");

    printf("\nLooking for %s", fname);

    orb_file = fopen(fname, "r");

    if (orb_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        return -1;
    }
    else
    {
        printf(" ..... Found !\n");
    }



    // open file to write parameters of domain and equation

    strcpy(fname, "output/");
    strcat(fname, outfname);

    if (timeinfo == 'i' || timeinfo == 'I')
    {
        strcat(fname, "_conf_imagtime.dat");
    } else
    {
        strcat(fname, "_conf_realtime.dat");
    }

    job_file = fopen(fname, "w");

    if (job_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        return -1;
    }

    fprintf(job_file, "# Trap name : %s\n", potname);



    EQ = SetupParams(eq_file, domain_file, potname, &dt, &N);
    nx = EQ->nx;
    ny = EQ->ny;
    S = carrDef(nx*ny); // solution at discretized positions
    abs2 = rarrDef(nx*ny);










    /*  ===============================================================
     
                       SETUP INITIAL DATA TO PROPAGATE
     
        ===============================================================  */

    for (i = 0; i < nx*ny; i++)
    {
        k = fscanf(orb_file, " (%lf%lfj)", &real, &imag);
        S[i] = real + I * imag;
        abs2[i] = real*real + imag*imag;
    }

    fclose(orb_file);

    printf("\nGot Initial condition, with norm =");
    printf(" %.10lf", sqrt(Rsimps2D(nx,ny,abs2,EQ->hx,EQ->hy)));
    printf(". nx = %d | ny = %d | hx = %.10lf | hy = %.10lf\n",nx,ny,EQ->hx,EQ->hy);










    /*  ===============================================================
     
                             CALL INTEGRATION ROUTINE
     
        ===============================================================  */

    start = omp_get_wtime();

    if (timeinfo == 'i' || timeinfo == 'I')
    {
        sepline();
        printf("\nDoing imaginary time integration #%d\n\n", 1);
        switch (method)
        {
            case 1:
                N = SplitStepPR(EQ, N, dt, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\n\nTime taken to solve with Peaceman-Rachford");
                printf(" : %.0lf seconds\n", time_used);
                break;
            case 2:
                N = SplitStepDYakonov(EQ, N, dt, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\n\nTime taken to solve with D'Yakonov");
                printf(" : %.0lf seconds\n", time_used);
                break;
        }

        // Record data
        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_line-1_orb_imagtime.dat");

        carr_txt(fname, nx*ny, S);

    }

    SaveConf(job_file, EQ, dt, N);





    /* release memory
     * ------------------------------------------------------------------- */

    fclose(job_file);
    fclose(eq_file);
    fclose(domain_file);
    free(S);
    free(abs2);
    ReleaseEqDataPkg(EQ);
    /* ------------------------------------------------------------------- */



    /*   ==========================    END    ==========================   */

    printf("\n\n");
    return 0;
}
