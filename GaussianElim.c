#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#pragma GCC optimize("O2")

#define OUTPUT_FORMAT  "%-8LG"
#define INPUT_FORMAT    "%Lf"

void printMatrix(long double *M, unsigned n, unsigned m);
void printMatrixAb(long double *, long double* , unsigned );


/* 
    expect: A have n^2 elements and b have n elements.
    make the simple Gaussian Elimination.
*/
void GaussianElim(long double *A, long double *b, unsigned n);

/* 
    expect: A have n^2 elements and b have n elements.
    make the Gaussian Elimination using partial pivoting.
*/
void GaussianElimParcialPivot(long double *A, long double* b, unsigned n);

/* 
    expect: A have n^2 elements and A is superior triangular matrix; b have n elements and X have n elements.
    compute the solution of using Backward Substitution method.
*/
void backwardSubstitution(long double *A, long double *b, long double *X, unsigned n);

unsigned main(void)
{
    register unsigned i, j;
    clock_t start, end;
    long double cpu_time_used = 0.0;
    unsigned n;

    scanf("%d", &n);

    long double *A = (long double*)malloc(sizeof(long double)*n*n);
    long double *b = (long double*)malloc(sizeof(long double)*n);
    long double *X = (long double*)malloc(sizeof(long double)*n);

    
    for(i = 0; i<n; i++)
        for(unsigned j = 0; j<n; j++){
            scanf(INPUT_FORMAT , &A[i * n + j]);
        }
    
    for(i = 0; i<n; i++)
        scanf(INPUT_FORMAT, &b[i]);

    

    start = clock();
    GaussianElimParcialPivot(A, b, n);
    //GaussianElim(A, b, n);
    end = clock();

    cpu_time_used = ((long double) (end - start)) / CLOCKS_PER_SEC;
    printf("***** %LG  seconds *****\n", cpu_time_used);

    
    printMatrixAb(A, b, n);


    start = clock();
    backwardSubstitution(A, b, X, n);
    end = clock();
    cpu_time_used = ((long double) (end - start)) / CLOCKS_PER_SEC;
    printf("***** %LG Seconds *****\n", cpu_time_used);

    printf("solution:\n");
    printMatrix(X, n, 1);


    free(A);
    free(b);
    free(X);

    return 0;
}

void printMatrix(long double *M, unsigned n, unsigned m)
{
    register int i;
    for(i = 0; i<n; i++)
    {
        for( unsigned j = 0; j<m; j++)
        {
            m==1?printf(OUTPUT_FORMAT, M[i] ) : printf(OUTPUT_FORMAT, M[n * i + j] );
        }
        printf("\n");
    }
    printf("\n");
}

void printMatrixAb(long double *A, long double *b, unsigned n)
{

    register int i, j, lin;

    for(i = 0; i<n; i++)
    {
        for(j = 0; j<n; j++)
        {
          printf(OUTPUT_FORMAT, A[n * i + j] );
        }
        printf("\t|\t");
        printf(OUTPUT_FORMAT, b[i] );
        printf("\n");
    }
    printf("\n");
}

void GaussianElim(long double *A, long double* b, unsigned n)
{
    register int i, j, k, lin;

    for(i = 0; i<n-1; i++)
    { 
        for(j = i+1; j<n; j++)
        {
            long double aux = A[ n * j + i ] / A[n * i + i];
            for(k = i+1; k<n; k++)
            {
                A[n * j + k] = A[n * j + k] - aux * A[n * i + k];
            }
            b[j] = b[j] - aux * b[i];
        }
    }

     for(i = 0; i<n; i++)
        for(lin = i+1; lin<n; lin++)
            A[ (n*lin) + i ]  = 0;
}

void GaussianElimParcialPivot(long double *A, long double *b, unsigned n)
{
    register int i, j, k, l, lin;

    for(i = 0; i<n-1; i++)
    {
        for(j = i+1; j<n; j++)
        {
            unsigned index = i + 1;
            long double pivot =  A[n * i + i];

            for( l = i+1; l<n; l++)
            {
                if( fabs(A[n * l + i]) > fabs(pivot) )
                {
                    pivot =  A[n * l + i]; 
                    index = l;
                }
            }
            if( pivot !=  A[n * i + i])
            {
                //Exchange of "lines" in "matrix" A
                long double *linhaAux = (long double*)malloc(sizeof(long double) * n );
                memcpy(linhaAux , A + n * i, sizeof(long double) * n );
                memcpy(A + n * i, A + n * index, sizeof(long double) * n );
                memcpy(A + n * index, linhaAux, sizeof(long double) * n );
                free(linhaAux);

                //Exchange of "lines" in vector b
                long double x = b[i];
                b[i] = b[index];
                b[index] = x;
            }

            long double aux = A[ n * j + i ] / A[n * i + i];
            for(k = i+1; k<n; k++)
            {
                A[n * j + k] = A[n * j + k] - aux * A[n * i + k];
            }
            b[j] = b[j] - aux * b[i];
        }
    }

     for(i = 0; i<n; i++)
        for(lin = i+1; lin<n; lin++)
            A[ (n*lin) + i ]  = 0;
}

void backwardSubstitution(long double* A, long double *b,  long double *X, unsigned n)
{
    register int lin, col;
    long double soma = 0.0;
    
    for( lin = n-1; lin >=0; lin--)
    {
        soma = 0.0;
        for(col = n-1; col>lin; col--)
        {
            soma = soma + A[n * lin + col] * X[col];
        }
        X[lin] = (b[lin] - soma) / A[lin * n + lin];
    }
}