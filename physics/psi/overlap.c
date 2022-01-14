/*
 * This program calculates the overlap of two atomic wave functions. The overlap
 * is defined by Sij = Integrate[ psi_i(r) * psi_j(r), dr ]. The atomic wave
 * functions are defined by the radial function and CombSphereHarm defined in
 * radial.c and ylm.c, so they are always real.
 *
 * This program requires an input file named overlap.inp which must be in the 
 * following form.
 * -----------------------------------------------------------------------------
 * 1 1 0 0     // Z1, n1, l1, m1
 * 1 1 0 0     // Z2. n2, l2, m2
 * 0. 0. 0.    // x1, y1, z1
 * 0. 0. 0.    // x2, y2, z2
 * 20. 0.05    // range of integration and step
 * 1.0E-9      // thresthold to avoid r = 0
 * -----------------------------------------------------------------------------
 *
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include "rnl.h"
#include "ylm.h"


int main( void )
{
    // flags in overlap.inp
    FILE *inp;
    int Zv1, n1, l1, m1;
    int Zv2, n2, l2, m2;
    double x1, y1, z1;
    double x2, y2, z2;
    double l, h;
    double eps;
    
    // variables realted to the integral
    int NG;
    double dv, integral;
    double xi, yi, *zi;
    double d1x, d1y, *d1z, d2x, d2y, *d2z, *r1, *r2;
    double *inti;

    // loop conters
    int i, j, k;


    // read overlap.inp
    inp = fopen( "overlap.inp", "r" );
    fscanf( inp, "%d%d%d%d\n",  &Zv1, &n1, &l1, &m1 );
    fscanf( inp, "%d%d%d%d\n",  &Zv2, &n2, &l2, &m2 );
    fscanf( inp, "%lf%lf%lf\n", &x1,  &y1, &z1 );
    fscanf( inp, "%lf%lf%lf\n", &x2,  &y2, &z2 );
    fscanf( inp, "%lf%lf\n",    &l,   &h );
    fscanf( inp, "%lf\n",       &eps );
    fclose( inp );
    
    // shift the coordinate system such that atom1 at (0,0,0)
    x2 = x2 - x1;
    y2 = y2 - y1;
    z2 = z2 - z1;
    x1 = x1 - x1;
    y1 = y1 - y1;
    z1 = z1 - z1;

    // determine meshgrid information
    NG = (int)(l/h);
    dv = pow( h, 3.0 );
    
    // the main loop
    integral = 0.0;
     zi  = (double *)malloc( sizeof(double) * (2*NG+1) );
    d1z  = (double *)malloc( sizeof(double) * (2*NG+1) );
    d2z  = (double *)malloc( sizeof(double) * (2*NG+1) );
     r1  = (double *)malloc( sizeof(double) * (2*NG+1) );
     r2  = (double *)malloc( sizeof(double) * (2*NG+1) );
    inti = (double *)malloc( sizeof(double) * (2*NG+1) );
    for( i = -NG; i <= NG; i ++ )
    {
        xi = h * i;
        d1x = xi - x1;
        d2x = xi - x2;
        for( j = -NG; j <= NG; j ++ )
        {
            yi = h * j;
            d1y = yi - y1;
            d2y = yi - y2;
#pragma omp parallel for
            for( k = 0; k <= 2*NG; k ++ )
            {
                zi[k]  = h * (k-NG);
                d1z[k] = zi[k] - z1;
                d2z[k] = zi[k] - z2;
                r1[k]  = sqrt( d1x*d1x + d1y*d1y + d1z[k]*d1z[k] );
                r2[k]  = sqrt( d2x*d2x + d2y*d2y + d2z[k]*d2z[k] );
                if( r1[k] <= eps || r2[k] <= eps )
                    inti[k] = 0.0;
                else
                    inti[k] = RadialFunc( Zv1, n1, l1, r1[k] ) * CombSphereHarm( l1, m1, d1x, d1y, d1z[k], r1[k] )
                            * RadialFunc( Zv2, n2, l2, r2[k] ) * CombSphereHarm( l2, m2, d2x, d2y, d2z[k], r2[k] );
            }
            for( k = 0; k <= 2*NG; k ++ )
                integral = integral + inti[k];
        }
    }
    free( zi   );
    free( d1z  );
    free( d2z  );
    free( r1   );
    free( r2   );
    free( inti );

    integral = integral * dv;
    printf( "%12s%24.16e\n", "overlap =", integral );
    return 0;
}
