#include <assert.h>
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
    compute the simple Regula Falsi method.
*/
long double Regulafalsi( long double(f)(int, ...), long double a, long double b, long double tol);

long double myfunction1(int, ...);
long double myfunction2(int, ...);

int main(void)
{
    long double a, b, x;

    a = -4;
    b = 4;

    x = Regulafalsi(myfunction2, a, b, TOL);

    printf("%LG\n", x);

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

long double Regulafalsi( long double(f)(int,...), long double a, long double b, long double tol)
{
    long double error = LONG_MAX;
    long double c = 0.0;

    while( error > tol)
    {
        c = ( a * f(1, b) - b * f(1,a) ) / ( f(1, b) - f(1, a) );
        if( f(1, a)*f(1,c) < 0 )
            b = c;
        else if( f(1, b)*f(1, c) < 0)
            a = c;
        error = fabs(f(1, c));
    }

    return c;
}