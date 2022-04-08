#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#define OUTPUT_FORMAT  "%-8LG"
#define INPUT_FORMAT    "%Lf"

#define TOL     1e-14L

/*
    expect: f and fl as function with variable number of arguments. furthermore fl is mathematically the derivative of f
    Compute the "root" of f using the Newton's method and return the "root"
*/
long double NewtonMethod( long double (*f)(int n, ...), long double (*df)(int n, ...), long double x0, long double tol);

/* f(x) =  x^5+x^4-3.3; root: 1.117329744559439*/
long double f(int n,...);

/* df(x) =  3*x^2 - 5 */
long double df(int n, ...);

int main(void)
{
    long double x0;

    //x0 = 1.1173; 
    x0 = 0.5;

    printf("root: ");
    printf(OUTPUT_FORMAT, NewtonMethod(f, df, x0, TOL) ); printf("\n");

    return 0;
}


long double NewtonMethod( long double (*f)(int n, ...), long double (*df)(int n, ...), long double x0, long double tol)
{
    long double k = 0.0;

    while ( fabs( f(1, x0) ) > tol )
    {
        k = x0 - f(1, x0) / df(1, x0);
        x0 = k;
    }

    return x0;
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

/* df(x) =  5*x^4 + 4*x^3 */
long double df(int n, ...)
{
    va_list arg_list;
    long double value = 0.0;
    long double x;

    va_start(arg_list, n);

    x = va_arg(arg_list, long double);
    value = 5 * powl(x, 4.0) + 4 * powl(x, 3.0);

    return value;
}