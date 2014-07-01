/* This function is called by the Ada program ts_cpu2norm_dd. */

#include <stdio.h>
#include <math.h>
#include "double_double.h"

int cpu2norm_dd_in_c ( int n, double *x, double *nrm )
{
    int i;
    double dd_nrm[2], re_acc[2], im_acc[2];

    printf("The C function received as dimension %d\n", n);
    printf("The components of the vector :\n");
    
    dd_nrm[0] = 0.0;
    dd_nrm[1] = 0.0;
    for(i=0; i<n; i++)
    {
        dd_write(&x[4*i], 32);    printf("  "); 
        dd_write(&x[4*i+2], 32);  printf("\n");
        dd_copy(&x[4*i], re_acc);   dd_sqr(re_acc, re_acc);
        dd_copy(&x[4*i+2], im_acc); dd_sqr(im_acc, im_acc);
        dd_inc(re_acc, im_acc);
        dd_inc(dd_nrm, re_acc);
    }
    printf("The square of the 2-norm is ");
    dd_write(dd_nrm, 32); printf("\n");

    nrm[0] = dd_nrm[0];
    nrm[1] = dd_nrm[1];

    return 0;
}
