#include <math.h>
#include "rnl.h"
double factor( int n )
{
    if ( n == 0 )
        return 1.0;
    else{
        double y;
        int i;
        y = 1.0;
        for( i = 1; i < n + 1; i ++ )
            y = y * (double)(i);
        return y;
    }
}

double ReNormFactor( int Z, int n, int l )
{
    double A, B, C;
    A = 2.0 / ( (double)(n*n) * factor(2*l+1));
    B = sqrt( factor(n+l) / factor(n-l-1) );
    C = sqrt( (double)(Z*Z*Z) );
    return A * B * C;
}

double KummerFunc( int alpha, int gamma, double ita )
{
    double sum, ck;
    int k;
    sum = 1.0;
    ck = 1.0;
    for( k = 1; k < -alpha+1; k ++ )
    {
        ck = ck * (double)(alpha+k-1) / (double)(gamma+k-1) * (1.0/(double)(k));
        sum = sum+ ck * pow( ita, (double)(k) );
    }
    return sum;
}

double RadialFunc( int Z, int n, int l, double r )
{
    double ita, NormFactor, A, B, C;
    ita = 2.0 * (double)(Z) / (double)(n) * r;
    NormFactor = ReNormFactor( Z, n, l );
    A = exp( -0.5 * ita  );
    B = pow( ita, (double)(l) );
    C = KummerFunc( -n+l+1, 2*l+2, ita );
    return NormFactor * A * B * C;
}
