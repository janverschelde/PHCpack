#include <stdio.h>
#include <cstdlib>
#include <gqd_type.h>
#include "double_double.h"
#include "DefineType.h"
#include "complexH.h"
#include "norm_host.h"
#include "norm_kernels.h"

extern "C" int gpu2norm_dd ( int dim, double *x, double *nrm )
/*
 * DESCRIPTION :
 *   A C++ function to compute the 2-norm of a vector,
 *   encapsulated as a C function for double doubles.
 *
 * ON ENTRY :
 *   dim      dimension of the complex vector in x;
 *   x        contains 4*dim doubles with the real and imaginary parts
 *            of the complex numbers of the vector.
 *
 * ON RETURN :
 *   x        x[0] and x[1] store 2-norm of the vector x;
 *   nrm      2-norm of the vector defined by x. */
{
    printf("Passing vector of dimension %d to the GPU...\n", dim);

    complexH<T1> v[dim];
    complex<T> v_h[dim];
    T1 twonorm;
    T twonorm_h;

    for(int i=0; i<dim; i++)
    {
       v[i].init(x[4*i],x[4*i+1],x[4*i+2],x[4*i+3]);
       v_h[i].initH(x[4*i],x[4*i+1],x[4*i+2],x[4*i+3]);
    }

    {
       double* nrm_on_cpu;
       nrm_on_cpu = (double*)calloc(2,sizeof(double));
       CPU_norm(v,dim,&twonorm);
       nrm_on_cpu = (double*) &twonorm;
       printf("Norm computed on host   : ");
       printf("%.15e\n", nrm_on_cpu[0]);
       // dd_write(nrm_on_cpu,32); printf("\n");
    }

    {
       double* nrm_on_gpu;
       nrm_on_gpu = (double*)calloc(2,sizeof(double));
       GPU_norm(v_h,dim,1,32,&twonorm_h);
       nrm_on_gpu = (double*) &twonorm_h;
       printf("Norm computed on device : ");
       printf("%.15e\n", nrm_on_gpu[0]);
       // dd_write(nrm_on_gpu,32); printf("\n");
    }

    // printf("twonorm_h[0] = %.15e\n", twonorm_h.x);
    // printf("twonorm_h[1] = %.15e\n", twonorm_h.y);
    nrm = (double*) (&twonorm_h);
    // printf("nrm[0] = %.15e\n", nrm[0]);
    // printf("nrm[1] = %.15e\n", nrm[1]);
    nrm[0] = (double) twonorm_h.x;
    nrm[1] = (double) twonorm_h.y;
    x[0] = nrm[0];
    x[1] = nrm[1];

    return 0;
}
