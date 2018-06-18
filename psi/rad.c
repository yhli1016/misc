/*
 * This program outputs the radial part of the wave function of hrdrogen-like
 * atom with specified Zval and quantum numbers n and l.
 *
 * This program requires an input file named rad.inp which must be in the 
 * following form.
 * -----------------------------------------------------------------------------
 * H21.dat    // output file name
 * 1          // Zval
 * 2 1        // n l
 * 50 0.1     // rmax and dr
 * -----------------------------------------------------------------------------
 *
 * The output has 3 columns: ri, Rnl, (ri*Rnl)^2
 *
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "rnl.h"


int main( void )
{
    // rad.inp
    char datname[32];
    int Zval;
    int n, l;
    double rmax, dr;

    // radial function
    double ri, Rnl;

    // files used through the program
    FILE *inp, *dat;
    
    
    // read rad.inp
    inp = fopen("rad.inp", "r");
    fscanf(inp, "%s\n", &datname);
    fscanf(inp, "%d\n", &Zval);
    fscanf(inp, "%d%d\n", &n, &l);
    fscanf(inp, "%lf%lf\n", &rmax, &dr);
    fclose(inp);
    
    // open output file
    dat = fopen(datname, "w");
    fprintf(dat, "# Z = %4d, n = %4d, l = %4d\n", Zval, n, l);
    
    // write radial distribution to output
    for (ri = 0.0; ri <= rmax; ri = ri + dr)
    {
        Rnl = RadialFunc(Zval, n, l, ri);
        fprintf(dat, "%14.4e%14.4e%14.4e\n", ri, Rnl, pow(ri * Rnl, 2.0));
    }
    
    // close output file
    fclose(dat);

    // exit
    return 0;
}
