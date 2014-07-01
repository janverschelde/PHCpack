/* This function is called by the Ada program ts_cpu2norm_d. */

#include <stdio.h>
#include <math.h>

int cpu2norm_d_in_c ( int n, double *x, double *nrm )
{
    int i;
    *nrm = 0.0;

    printf("The C function received as dimension %d\n", n);
    printf("The components of the vector :\n");
    for(i=0; i<n; i++)
    {
        printf("%22.15le  %22.15le\n",x[2*i],x[2*i+1]);
        *nrm = *nrm + x[2*i]*x[2*i] + x[2*i+1]*x[2*i+1];
    }
    *nrm = sqrt(*nrm);
    printf("The 2-norms is %.15le\n", *nrm);

    return 0;
}
