#include "checkinp.h"

int checkinp( int l, int m, int flavor )
{
    if ( l < 0 )
        return -1;
    else if ( l > 3 )
        return -2;
    else if ( m < 0 && -m > l )
        return -3;
    else if ( m >= 0 && m > l )
        return -3;
    else if ( flavor != -1 && flavor != 1 && flavor != 0 )
        return -4;
    else
        return 0;
}
