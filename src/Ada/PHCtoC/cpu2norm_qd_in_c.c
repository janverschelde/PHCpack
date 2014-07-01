/* This function is called by the Ada program ts_cpu2norm_qd. */

#include <stdio.h>
#include <math.h>
#include "quad_double.h"

int cpu2norm_qd_in_c ( int n, double *x, double *nrm )
{
    int i;
    double qd_nrm[4], re_acc[4], im_acc[4];

    printf("The C function received as dimension %d\n", n);
    printf("The components of the vector :\n");
    
    qd_nrm[0] = 0.0;
    qd_nrm[1] = 0.0;
    qd_nrm[2] = 0.0;
    qd_nrm[3] = 0.0;
    for(i=0; i<n; i++)
    {
        qd_write(&x[8*i], 64);    printf("  "); 
        qd_write(&x[8*i+4], 64);  printf("\n");
        qd_copy(&x[8*i], re_acc);   qd_sqr(re_acc, re_acc);
        qd_copy(&x[8*i+4], im_acc); qd_sqr(im_acc, im_acc);
        qd_add(re_acc, im_acc, re_acc);
        qd_add(qd_nrm, re_acc, qd_nrm);
    }
    printf("The square of the 2-norm is ");
    qd_write(qd_nrm, 64); printf("\n");

    nrm[0] = qd_nrm[0];
    nrm[1] = qd_nrm[1];
    nrm[2] = qd_nrm[2];
    nrm[3] = qd_nrm[3];

    return 0;
}
