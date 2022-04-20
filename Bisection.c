#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#define TOL     1e-6L

#define OUTPUT_FORMAT  "%-8LG"
#define INPUT_FORMAT    "%Lf"


/* 
    expect: f as function with variable number of arguments.
    compute the simple Bisection method.
*/
long double bissec( long double(f)(int, ...), long double a, long double b, long double tol, unsigned int *niter);

long double myfunction1(int, ...);
long double myfunction2(int, ...);

int main(void)
{
    long double a, b, x;
    unsigned int niter = 0;

    a = -4;
    b = 4;

    x = bissec(myfunction2, a, b, TOL, &niter);

    printf("%u iterations\n", niter);
    printf("x = %LG\n", x);

    return 0;
}

/* f(x) = x^5 + x^4 - 3.3*/
long double myfunction1(int n, ...)
{
    long double value = 0.0;
    va_list arg_list;

    va_start(arg_list, n);

    long double x = va_arg(arg_list, long double);

    value = powl(x,5.0) + powl(x, 4.0) - 3.3;

    return value;
}

/* f(x) = (x−2)^2 * (x+1) * (x+5)*/
long double myfunction2(int n, ...)
{
    long double value = 0.0;
    va_list arg_list;

    va_start(arg_list, n);

    long double x = va_arg(arg_list, long double);

    value = powl(x-2, 2.0) * (x + 1) * ( x + 5);

    va_end(arg_list);

    return value;
}

long double bissec( long double(f)(int, ...), long double a, long double b, long double tol, unsigned int *niter)
{
    long double error = LONG_MAX;
    long double c = 0.0;
    unsigned int cont = 0;

    while( error > tol)
    {
        c = (a+b)/2;
        if( f(1, a)*f(1,c) < 0 )
            b = c;
        else if( f(1, b)*f(1, c) < 0)
            a = c;
        error = fabs(f(1, c));
        
        cont++;
    }
    
    *niter = cont;

    return c;
}
