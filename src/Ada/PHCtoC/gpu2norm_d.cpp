#include <stdio.h>
#include <gqd_type.h>
#include "DefineType.h"
#include "complexH.h"
#include "norm_host.h"
#include "norm_kernels.h"

extern "C" int gpu2norm_d ( int dim, double *x, double *nrm )
/*
 * DESCRIPTION :
 *   A C++ function to compute the 2-norm of a vector,
 *   encapsulated as a C function for ordinary doubles.
 *
 * ON ENTRY :
 *   dim      dimension of the complex vector in x;
 *   x        contains 2*dim doubles with the real and imaginary parts
 *            of the complex numbers of the vector.
 *
 * ON RETURN :
 *   nrm      2-norm of the vector defined by x. */
{
    printf("Passing vector of dimension %d to the GPU...\n", dim);

    complexH<T1> v[dim];
    complex<T> v_h[dim];
    T1 twonorm;
    T twonorm_h;

    for(int i=0; i<dim; i++)
    {
       v[i].init(x[2*i],x[2*i+1]);
       v_h[i].initH(x[2*i],x[2*i+1]);
    }

    {
       double nrm_on_cpu;
       CPU_norm(v,dim,&twonorm);
       nrm_on_cpu = (double) twonorm;
       printf("Norm computed on host : %.15le\n", nrm_on_cpu);
    }

    {
       double nrm_on_gpu;
       GPU_norm(v_h,dim,1,32,&twonorm_h);
       nrm_on_gpu = (double) twonorm_h;
       printf("Norm computed on device : %.15le\n", nrm_on_gpu);
    }

    *nrm = (double) twonorm_h;

    return 0;
}
