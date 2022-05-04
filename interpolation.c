#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define OUTPUT_FORMAT  "%-8LG"
#define INPUT_FORMAT    "%Lf"

//Polynomial representation: 2x^3 + 8x + 5 <==> P[0] = 5, P[1] = 8, P[2] = 0, P[3] = 2;

void printMatrix(long double *M, unsigned n, unsigned m);
void printMatrixAb(long double *, long double* , unsigned );

/* 
    expect: A have n^2 elements and b have n elements.
    make the simple Gaussian Elimination.
*/
void GaussianElim(long double *A, long double *b, unsigned n);

/* 
    expect: A have n^2 elements and A is superior triangular matrix; b have n elements and X have n elements.
    compute the solution of using Backward Substitution method.
*/
void backwardSubstitution(long double *A, long double *b, long double *X, unsigned n);

void printpolynomial(long double *p, int n);
long double Valueofpolynomial(long double *p, long double x, int n);
void buildMatrix(long double *x, long double *y, 
                 long double *M, long double *b, int n);
int getDataFromfile(char *filename, long double **x, long double **y);
void printTable(long double *x, long double *y, int n);

int main(void)
{

   char str[256];
   long double *x = NULL;
   long double *y = NULL;
   long double k = 0.0;

   printf("inputfile: ");
   scanf("%s", str);
   int n = getDataFromfile( str, &x, &y);
   int degree = n-1;

   long double *P = (long double*)malloc(sizeof(long double)*n);
   long double *M = (long double*)malloc(sizeof(long double)*n*n);
   long double *b = (long double*)malloc(sizeof(long double)*n);
   
   printTable(x, y, n);

   buildMatrix(x, y, M, b, n);
   printMatrixAb( M, b , (unsigned)n );

   GaussianElim(M, b, n);
   printMatrixAb( M, b , (unsigned)n );

   backwardSubstitution(M, b, P, (unsigned)n);

   printpolynomial(P, (unsigned)n);

   printf(">");
   scanf("%Lf", &k);
   printf("P(%LG) = %LG\n", k, Valueofpolynomial(P, k, n) ); 


   free(x);
   free(y);
   free(P);
   free(M);
   free(b);

   return 0;
}

void buildMatrix(long double *x, long double *y, 
                 long double *M, long double *b, int n)
{

   for( int i = 0; i<n; i++)
   {
      b[i] = y[i];
      for( int j = 0; j<n; j++)
         M[n * i + j] = pow(x[i], j);
   }
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

long double Valueofpolynomial(long double *p, long double x, int n)
{
   long double y = 0.0;

   for( int i = 0; i<n; i++)
      y = y + p[i] * pow(x,i);

   return y;
}

void printpolynomial(long double *p, int n)
{
   printf("\n");

   printf("P(x) = ");
   printf("%+LG", p[0]);
   printf("%+LG*x", p[1]);
   for( int i = 2; i<n; i++)
      printf("%+LG*x^%d", p[i], i);

   printf("\n\n");

}

int getDataFromfile(char *filename, long double **x, long double **y)
{
   int n = 0;
   FILE* arq = fopen(filename, "r");

   if( arq == NULL)
      exit(1);

   fscanf(arq, "%d", &n);

   (*x) = (long double*)malloc(sizeof(long double)*n);
   (*y) = (long double*)malloc(sizeof(long double)*n);
  

   for(int i = 0; i<n; i++)
      fscanf(arq, "%Lf", &((*x)[i]));

   for(int i = 0; i<n; i++)
      fscanf(arq, "%Lf", &((*y)[i]));

   fclose(arq);

   return n;
   
}

void printTable(long double *x, long double *y, int n)
{
   printf("\n");
   printf("%5c|", 'x');
   for(int i = 0; i<n; i++)
      printf("%5LG|", x[i]);

   printf("\n");
   printf("%5s|", "f(x)");
   for( int i = 0; i<n; i++)
      printf("%5LG|", y[i]);

   printf("\n\n");
  
}
