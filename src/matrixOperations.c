#include "matrixOperations.h"



/*          ***********************************************

                        SETUP VALUES IN MATRIX

            ***********************************************          */



void cmatFill(int m, int n, double complex z, Cmatrix M)
{

/** Fill all matrix with a single constant value **/

    unsigned int
        i,
        j;

    for (i = 0; i < m; i++) { for (j = 0; j < n; j++) M[i][j] = z; }
}



void cmatFillDK(int n, int k, Carray z, Cmatrix M)
{

/** Fill a diagonal starting from k-column if k >= 0
  * Fill a diagonal starting from |k|-row  if k <  0
  *
  * of a square matrix of dimension 'n' with an array
  * that must have size of (n - |k|) at least.    **/

    unsigned int i;

    if (k >= 0) { for (i = 0; i < n - k; i++) M[i][i+k] = z[i]; }
    else        { for (i = 0; i < n + k; i++) M[i-k][i] = z[i]; }
}



void cmatFillTri(int n, Carray upper, Carray mid, Carray lower, Cmatrix M)
{

/** Fill a tridiagonal matrix **/

    cmatFill(n, n, 0, M);
    cmatFillDK(n, -1, lower, M);
    cmatFillDK(n,  1, upper, M);
    cmatFillDK(n,  0, mid, M);
}



void RowMajor(int m, int n, Cmatrix M, Carray v)
{

/** Copy data from matrix to vector using row major layout
  * The size of v is required to be at least m*n       **/

    int i,
        j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++) v[i * n + j] = M[i][j];
    }
}



void setValueCCScmat(int n, int i, int j, int col, doublec z, CCScmat M)
{

/** Set the element of row 'i' and column 'col' in a CCS matrix of  complex
  * numbers. It is required the number 'j' of how many nonzero elements are
  * behind it in the same row. 'j' cannot exceed the max number of nonzeros
  * pre-defined in the structure.                                       **/

    if (j >= M->m)
    {
        printf("\n\nSETUP ERROR : Trying to set element in column that");
        printf(" exceed pre-defined maximum number of nonzero elements");
        printf(" in a same row.\n\n");
        exit(EXIT_FAILURE);
    }

    M->vec[i + n * j] = z;
    M->col[i + n * j] = col;
}



void setValueCCSrmat(int n, int i, int j, int col, double x, CCSrmat M)
{

/** Set the element of row 'i' and column 'col' in a CCS matrix of  complex
  * numbers. It is required the number 'j' of how many nonzero elements are
  * behind it in the same row. 'j' cannot exceed the max number of nonzeros
  * pre-defined in the structure.                                       **/

    if (j >= M->m)
    {
        printf("\n\nSETUP ERROR : Trying to set element in column that");
        printf(" exceed pre-defined maximum number of nonzero elements");
        printf(" in a same row.\n\n");
        exit(EXIT_FAILURE);
    }

    M->vec[i + n * j] = x;
    M->col[i + n * j] = col;
}





/*          ***********************************************

             SPECIAL ROUTINES FOR FINITE DIFFERENCE SCHEME

            ***********************************************          */





CCScmat tri2CCS(int n, Carray upper, Carray lower, Carray mid)
{

/** Configure a CCS matrix from a tridigonal one, being passed through
  * 3  vectors  correponding to diagonals.  Return the address  of the
  * structure allocated **/

    unsigned int j;

    CCScmat M;

    M = ccscmatDef(n,3);

    // first and last rows have just 2 nonzero elements
    // then must be configured separately
    M->vec[0]         = mid[0];
    M->vec[n]         = upper[0];
    M->vec[2 * n]     = 0;
    M->vec[3 * n - 1] = 0;

    for (j = 1; j < n; j++)                 { M->vec[j] = lower[j - 1];     }
    for (j = n + 1; j < 2 * n; j++)         { M->vec[j] = mid[j - n];       }
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->vec[j] = upper[j - 2 * n]; }

    // first and last rows have just 2 nonzero elements
    // then must be configured separately
    M->col[0]         = 0;
    M->col[n]         = 1;
    M->col[2 * n]     = 0;
    M->col[3 * n - 1] = 0;

    for (j = 1; j < n; j++)                 { M->col[j] = j - 1; }
    for (j = n + 1; j < 2 * n; j++)         { M->col[j] = j - n; }
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->col[j] = j + 1 - 2 * n; }

    return M;
}





/*          **********************************************

            MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION

            **********************************************          */





void cmatvec(int m, int n, Cmatrix M, Carray v, Carray ans)
{

    unsigned int
        i,
        j;

    double complex
        z;

    for (i = 0; i < m; i++)
    {
        z = 0 + 0 * I;
        for (j = 0; j < n; j++) z = z + M[i][j] * v[j];
        ans[i] = z;
    }
}



void rmatvec(int m, int n, Rmatrix M, Rarray v, Rarray ans)
{

    unsigned int
        i,
        j;

    double
        x;

    for (i = 0; i < m; i++)
    {
        x = 0;
        for (j = 0; j < n; j++) x = x + M[i][j] * v[j];
        ans[i] = x;
    }
}



void cmatmat(int m, int n, int l, Cmatrix M, Cmatrix A, Cmatrix ans)
{

    unsigned int
        i,
        j,
        k;

    double complex
        z;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < l; j++)
        {
            z = 0 + 0 * I;
            for (k = 0; k < n; k++) z = z + M[i][k] * A[k][j];
            ans[i][j] = z;
        }
    }
}



void rmatmat(int m, int n, int l, Rmatrix M, Rmatrix A, Rmatrix ans)
{

    unsigned int
        i,
        j,
        k;

    double
        x;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < l; j++)
        {
            x = 0;
            for (k = 0; k < n; k++) x = x + M[i][k] * A[k][j];
            ans[i][j] = x;
        }
    }
}



void CCScmatvec(int n, Carray vals, int * cols, int m, Carray vec, Carray ans)
{

    unsigned int
        i,
        l;

    double complex
        z;

    #pragma omp parallel for private(l, i, z)
    for (i = 0; i < n; i++)
    {
        z = vals[i] * vec[cols[i]];
        for (l = 1; l < m; l++) z = z + vals[i + l*n] * vec[cols[i + l*n]];
        ans[i] = z;
    }
}



void CCSrmatvec(int n, Rarray vals, int * cols, int m, Rarray vec, Rarray ans)
{
    unsigned int
        i,
        l;

    double
        x;

    #pragma omp parallel for private(l, i, x)
    for (i = 0; i < n; i++)
    {
        x = vals[i] * vec[cols[i]];
        for (l = 1; l < m; l++) x = x + vals[i+l*n] * vec[cols[i+l*n]];
        ans[i] = x;
    }
}



/*          ***********************************************

                          Inversion of matrices

            ***********************************************          */



int HermitianInv(int M, Cmatrix A, Cmatrix A_inv)
{

/** Use Lapack routine to solve systems of equations with the
  * right-hand-side being identity matrix  to get the inverse **/

    int i, // counter
        j, // counter
        l; // lapack success parameter

    int
        * ipiv;

    CMKLarray
        ArrayForm, // To call zhesv routine use row major layout of Matrix
        Id;        // Identity matrix in row major layout



    ipiv = (int *) malloc(M * sizeof(int));

    ArrayForm = cmklDef(M * M);

    Id = cmklDef(M * M);



    for (i = 0; i < M; i++)
    {
        // Setup (L)ower triangular part as a Row-Major-Array to use lapack
        ArrayForm[i * M + i].real = creal(A[i][i]);
        ArrayForm[i * M + i].imag = 0;
        Id[i * M + i].real = 1;
        Id[i * M + i].imag = 0;

        for (j = 0; j < i; j++)
        {
            ArrayForm[i * M + j].real = creal(A[i][j]);
            ArrayForm[i * M + j].imag = cimag(A[i][j]);

            ArrayForm[j * M + i].real = 0; // symbolic values
            ArrayForm[j * M + i].imag = 0; // for upper triangular part

            Id[i * M + j].real = 0;
            Id[i * M + j].imag = 0;

            Id[j * M + i].real = 0;
            Id[j * M + i].imag = 0;
        }
    }

    l = LAPACKE_zhesv(LAPACK_ROW_MAJOR, 'L', M, M, ArrayForm, M, ipiv, Id, M);

    for (i = 0; i < M; i++)
    {
        // Transcript the result back to matrix form
        for (j = 0; j < M; j++)
        {
            A_inv[i][j] = Id[i * M + j].real + I * Id[i * M + j].imag;
        }
    }

    free(ipiv);
    free(Id);
    free(ArrayForm);

    return l;
}
