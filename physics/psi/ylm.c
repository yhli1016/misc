#include <math.h>
#include "ylm.h"

double coeff( int a, int b, int c, int d )
{
    const double pi = 3.14159265358979323846;
    return (double)(a)/(double)(b) * sqrt(c/(d*pi));
}

double RealSphereHarm( int l, int m, double theta, double phi )
{
    double ct = cos( theta );
    double st = sin( theta );
    double cmp = cos( m * phi );
    
    if (l == 0)
        return coeff( 1, 2, 1, 1 );
    else if (l == 1)
        switch(m)
        {
            case -1: return coeff(  1, 2, 3, 2 ) * cmp * st; break;
            case  0: return coeff(  1, 2, 3, 1 ) * cmp * ct; break;
            case  1: return coeff( -1, 2, 3, 2 ) * cmp * st; break;
            default: return 0.0; break;
        }
    else if (l == 2)
        switch(m)
        {
            case -2: return coeff(  1, 4, 15, 2 ) * cmp * st * st; break;
            case -1: return coeff(  1, 2, 15, 2 ) * cmp * st * ct; break;
            case  0: return coeff(  1, 4,  5, 1 ) * cmp * ( 3.0*ct*ct - 1.0 ); break;
            case  1: return coeff( -1, 2, 15, 2 ) * cmp * st * ct; break;
            case  2: return coeff(  1, 4, 15, 2 ) * cmp * st * st; break;
            default: return 0.0; break;
        }
    else if (l == 3)
        switch(m)
        {
            case -3: return coeff(  1, 8,  35, 1 ) * cmp * st * st * st; break;
            case -2: return coeff(  1, 4, 105, 2 ) * cmp * st * st * ct; break;
            case -1: return coeff(  1, 8,  21, 1 ) * cmp * st * ( 5.0*ct*ct - 1.0 ); break;
            case  0: return coeff(  1, 4,   7, 1 ) * cmp * ct * ( 5.0*ct*ct - 3.0 ); break;
            case  1: return coeff( -1, 8,  21, 1 ) * cmp * st * ( 5.0*ct*ct - 1.0 ); break;
            case  2: return coeff(  1, 4, 105, 2 ) * cmp * st * st * ct; break;
            case  3: return coeff( -1, 8,  35, 1 ) * cmp * st * st * st; break;
            default: return 0.0; break;
        }
    else
        return 0.0;
}

double ImagSphereHarm( int l, int m, double theta, double phi )
{
    double ct = cos( theta );
    double st = sin( theta );
    double smp = sin( m * phi );
    
    if (l == 0)
        return coeff( 1, 2, 1, 1 );
    else if (l == 1)
        switch(m)
        {
            case -1: return coeff(  1, 2, 3, 2 ) * smp * st; break;
            case  0: return coeff(  1, 2, 3, 1 ) * smp * ct; break;
            case  1: return coeff( -1, 2, 3, 2 ) * smp * st; break;
            default: return 0.0; break;
        }
    else if (l == 2)
        switch(m)
        {
            case -2: return coeff(  1, 4, 15, 2 ) * smp * st * st; break;
            case -1: return coeff(  1, 2, 15, 2 ) * smp * st * ct; break;
            case  0: return coeff(  1, 4,  5, 1 ) * smp * ( 3.0*ct*ct - 1.0 ); break;
            case  1: return coeff( -1, 2, 15, 2 ) * smp * st * ct; break;
            case  2: return coeff(  1, 4, 15, 2 ) * smp * st * st; break;
            default: return 0.0; break;
        }
    else if (l == 3)
        switch(m)
        {
            case -3: return coeff(  1, 8,  35, 1 ) * smp * st * st * st; break;
            case -2: return coeff(  1, 4, 105, 2 ) * smp * st * st * ct; break;
            case -1: return coeff(  1, 8,  21, 1 ) * smp * st * ( 5.0*ct*ct - 1.0 ); break;
            case  0: return coeff(  1, 4,   7, 1 ) * smp * ct * ( 5.0*ct*ct - 3.0 ); break;
            case  1: return coeff( -1, 8,  21, 1 ) * smp * st * ( 5.0*ct*ct - 1.0 ); break;
            case  2: return coeff(  1, 4, 105, 2 ) * smp * st * st * ct; break;
            case  3: return coeff( -1, 8,  35, 1 ) * smp * st * st * st; break;
            default: return 0.0; break;
        }
    else
        return 0.0;
}

// Formulas obtained from the URL: 
// https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
double CombSphereHarm( int l, int m, double x, double y, double z, double r )
{
    if (l==0)
        return coeff( 1, 2, 1, 1 );
    else if (l==1)
    {
        double denorm = r;
        switch (m)
        {
            case -1: return coeff( 1, 1, 3, 4 ) / denorm * y; break; // py
            case  0: return coeff( 1, 1, 3, 4 ) / denorm * z; break; // pz
            case  1: return coeff( 1, 1, 3, 4 ) / denorm * x; break; // px
            default: return 0.0; break;
        }
    }
    else if (l==2)
    {
        double denorm = pow( r, 2.0 );
        switch (m)
        {
            case -2: return coeff( 1, 2, 15, 1 ) / denorm * ( x*y ); break; // dxy
            case -1: return coeff( 1, 2, 15, 1 ) / denorm * ( y*z ); break; // dyz
            case  0: return coeff( 1, 4,  5, 1 ) / denorm * ( 2*z*z - x*x - y*y ); break; // dz^2
            case  1: return coeff( 1, 2, 15, 1 ) / denorm * ( z*x ); break; // dxz
            case  2: return coeff( 1, 4, 15, 1 ) / denorm * ( x*x - y*y ); break; // d(x^2-y^2)
            default: return 0.0; break;
        }
    }
    else if (l==3)
    {
        double denorm = pow( r, 3.0 );
        switch (m)
        {
            case -3: return coeff( 1, 4,  35, 2 ) / denorm * y * ( 3*x*x -y*y ); break; // fy(3x^2-y^2)
            case -2: return coeff( 1, 2, 105, 1 ) / denorm * ( x*y*z ); break; // fxyz
            case -1: return coeff( 1, 4,  21, 2 ) / denorm * y * ( 4*z*z - x*x - y*y ); break; // fyz^2
            case  0: return coeff( 1, 4,   7, 1 ) / denorm * z * ( 2*z*z - 3*x*x - 3*y*y ); break; // fz^3
            case  1: return coeff( 1, 4,  21, 2 ) / denorm * x * ( 4*z*z - x*x - y*y ); break; // fxz^2
            case  2: return coeff( 1, 4, 105, 1 ) / denorm * z * ( x*x - y*y ); break; // fz(x^2-y^2)
            case  3: return coeff( 1, 4,  35, 2 ) / denorm * x * ( x*x - 3*y*y ); break; // fx(x^2-3y^2)
            default: return 0.0; break;
        }
    }
    else
        return 0.0;
}
