#include <stdio.h>
#include <cstdlib>
#include <gqd_type.h>
#include "DefineType.h"
#include "complexH.h"
#include "norm_host.h"
#include "norm_kernels.h"

extern "C" int gpu2norm_qd ( int dim, double *x, double *nrm )
/*
 * DESCRIPTION :
 *   A C++ function to compute the 2-norm of a vector,
 *   encapsulated as a C function for quad doubles.
 *
 * ON ENTRY :
 *   dim      dimension of the complex vector in x;
 *   x        contains 8*dim doubles with the real and imaginary parts
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

    // printf("x[0] = %.15e\n", x[0]);
    // printf("x[1] = %.15e\n", x[1]);
    // printf("x[2] = %.15e\n", x[2]);
    // printf("x[3] = %.15e\n", x[3]);

    for(int i=0; i<dim; i++)
    {
       v[i].init(x[8*i],x[8*i+1],x[8*i+2],x[8*i+3],
                 x[8*i+4],x[8*i+5],x[8*i+6],x[8*i+7]);
       v_h[i].initH(x[8*i],x[8*i+1],x[8*i+2],x[8*i+3],
                    x[8*i+4],x[8*i+5],x[8*i+6],x[8*i+7]);
       // v[i].init(x[8*i],x[8*i+4]);
       // v_h[i].initH(x[8*i],x[8*i+4]);
    }

    {
       double* nrm_on_cpu;
       nrm_on_cpu = (double*)calloc(4,sizeof(double));
       CPU_norm(v,dim,&twonorm);
       nrm_on_cpu = (double*) &twonorm;
       printf("Norm computed on host   : ");
       printf("cpu[0] = %.15e\n", nrm_on_cpu[0]);
       printf("cpu[1] = %.15e\n", nrm_on_cpu[1]);
       printf("cpu[2] = %.15e\n", nrm_on_cpu[2]);
       printf("cpu[3] = %.15e\n", nrm_on_cpu[3]);
       // x[0] = nrm_on_cpu[0];
       // x[1] = nrm_on_cpu[1];
       // x[2] = nrm_on_cpu[2];
       // x[3] = nrm_on_cpu[3];
       // dd_write(nrm_on_cpu,32); printf("\n");
    }

    {
       double* nrm_on_gpu;
       nrm_on_gpu = (double*)calloc(4,sizeof(double));
       GPU_norm(v_h,dim,1,32,&twonorm_h);
       nrm_on_gpu = (double*) &twonorm_h;
       printf("Norm computed on device : ");
       printf("gpu[0] = %.15e\n", nrm_on_gpu[0]);
       printf("gpu[1] = %.15e\n", nrm_on_gpu[1]);
       printf("gpu[2] = %.15e\n", nrm_on_gpu[2]);
       printf("gpu[3] = %.15e\n", nrm_on_gpu[3]);
       // dd_write(nrm_on_gpu,32); printf("\n");
    }

    // printf("twonorm_h[0] = %.15e\n", twonorm_h.x);
    // printf("twonorm_h[1] = %.15e\n", twonorm_h.y);
    // printf("twonorm_h[2] = %.15e\n", twonorm_h.z);
    // printf("twonorm_h[3] = %.15e\n", twonorm_h.w);
    nrm = (double*) (&twonorm_h);
    // printf("nrm[0] = %.15e\n", nrm[0]);
    // printf("nrm[1] = %.15e\n", nrm[1]);
    // printf("nrm[2] = %.15e\n", nrm[2]);
    // printf("nrm[3] = %.15e\n", nrm[3]);
    nrm[0] = (double) twonorm_h.x;
    nrm[1] = (double) twonorm_h.y;
    nrm[2] = (double) twonorm_h.z;
    nrm[3] = (double) twonorm_h.w;
    x[0] = nrm[0];
    x[1] = nrm[1];
    x[2] = nrm[2];
    x[3] = nrm[3];

    return 0;
}
