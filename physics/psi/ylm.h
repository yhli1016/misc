#ifndef _YLM_H
#define _YLM_H
double RealSphereHarm( int l, int m, double theta, double phi );
double ImagSphereHarm( int l, int m, double theta, double phi );
double CombSphereHarm( int l, int m, double x, double y, double z, double r );
#endif
