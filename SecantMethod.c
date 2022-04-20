#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#define OUTPUT_FORMAT  "%-8LF"
#define INPUT_FORMAT    "%Lf"

#define TOL     1e-6L

/*
    Compute the "root" of f using the Secant method and return the "root"
*/
long double SecantMethod( long double (*f)(int n, ...), long double x0, long double x1, long double tol, unsigned int *niter);

/* f(x) =  x^5+x^4-3.3; root: 1.117329744559439*/
long double f(int n,...);


int main(void)
{
    long double x0, x1, x;
    unsigned int niter = 0;

    x0 = 2;
    x1 = 3; 

    
    x = SecantMethod(f, x0, x1, TOL, &niter);
    printf("%u iterations\n", niter);
    printf("x = ");
    printf(OUTPUT_FORMAT, x ); printf("\n");

    return 0;
}


long double SecantMethod( long double (*f)(int n, ...), long double x0, long double x1, long double tol, unsigned int *niter)
{
    long double k = 0.0;
    unsigned int cont = 0;

    while( fabs( f(1, x1) ) > tol )
    {
       
        k = x1 - f(1, x1) * (x0 - x1) / ( f(1, x0) - f(1, x1) );
        x0 = x1;
        x1 = k;
        cont++;
    }

    *niter = cont;
    return x1;
}

/* f(x) =  x^5 + x^4 - 3.3 */
long double f(int n,...)
{
    va_list arg_list;
    long double value = 0.0;
    long double x;
    
    va_start(arg_list, n);

    x = va_arg(arg_list, long double);

    value = powl(x, 5.0) + powl(x, 4.0) - 3.3;

    return value;
}
