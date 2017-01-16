// This file contains all PED kernels for all three precisions.
// The corresponding prototypes are in "ped_kernels.h".
// In order to isolate the separate compilation by nvcc from the linking
// by gcc with the main test program, this file contains all three instances
// of the GPU_evaldiff function, for all three precisions.
// The definitions that are common in ped_kernels.cu and ped_kernels_qd.cu
// are placed at beginning of this file.
// For the double and double double precision, the code in ped_kernels.cu
// is copied twice, once with realD replaced by double, and once with realD
// replaced by gdd_real.  For quad double precision, realD is replaced by
// gqd_real in the code copied from ped_kernels_qd.cu.

#include <iostream>
#include <gqd_type.h>
#include "gqd.cu"
#include "complexD.h"

using namespace std;

__constant__ char positions[20000];
__constant__ char exponents[20000];

inline int number_of_blocks ( int dim, int BS )
/* Computes the number of blocks based on the dimension
   and the number of threads in a block. */
{
   int nblocks;

   if((dim*dim + dim) % BS == 0)
      nblocks = (dim*dim + dim)/BS;
   else
      nblocks = (dim*dim + dim)/BS + 1;

   return nblocks;
}

// for double precision

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<double> *xval, complexD<double> *factors );

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<double> *xval, complexD<double> *factors );

__global__ void sum_monoms
 ( int dim_s, complexD<double> *monvalues, complexD<double> *polvalues,
   int act_n_threads );

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<double> *xval, complexD<double> *roots, complexD<double> *coefs,
   complexD<double> *monvalues, complexD<double> *monderivatives );

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<double> *x_h, complexD<double> *c_h,
   complexD<double> *factors_h, complexD<double> *polvalues_h )
{
   int nblocks = number_of_blocks(dim,BS);
   complexD<double> *derivatives_d;
   complexD<double> *factors_d;
   complexD<double> *monvalues_d;
   complexD<double> *polvalues_d;
   complexD<double> *x_d;
   complexD<double> *c_d;
   // allocate space for output
   int ant = ((dim*dim+dim)/BS + 1)*BS;
   int aas = ant*m;
   //complexD<double> *derivatives_h = new complexD<double>[aas];
   complexD<double> derivatives_h[aas];
   // cudaMalloc((void**)&derivatives_h,aas);
   for(int i=0; i<aas; i++)
      derivatives_h[i].initH(0.0,0.0);
   size_t size_c = ncoefs*sizeof(complexD<double>);
   size_t size_d = aas*sizeof(complexD<double>);
   size_t size_pols = ant*sizeof(complexD<double>);
   // copy positions and exponents to constant memory
   cudaMemcpyToSymbol(positions,pos,NM*NV*sizeof(char));
   cudaMemcpyToSymbol(exponents,exp,NM*NV*sizeof(char));
   size_t size_m = NM*sizeof(complexD<double>);
   cudaMalloc((void**)&monvalues_d,size_m);

   size_t size = dim*sizeof(complexD<double>);
   cudaMalloc((void**)&x_d,size);
   cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
   size_t size_NM = NM*sizeof(complexD<double>);
   cudaMalloc((void**)&factors_d,size_NM);
   cudaMalloc((void**)&c_d,size_c);
   cudaMemcpy(c_d,c_h,size_c,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&derivatives_d,size_d);
   cudaMalloc((void**)&polvalues_d,size_pols);
   cudaMemcpy(derivatives_d,derivatives_h,size_d,cudaMemcpyHostToDevice);
   for(int j=0; j<r; j++)
   {
       mult1_sw_ind_for_shar<<<NM/BS,BS>>>(dim,deg,NV,x_d,factors_d);
       speeldif<<<NM/BS,BS>>>(dim,NM,ant,m,NV,x_d,factors_d,c_d,
                              monvalues_d,derivatives_d);
       sum_monoms<<<nblocks,BS>>>(m,derivatives_d,polvalues_d,ant);
   }
   cudaMemcpy(factors_h,factors_d,size_NM,cudaMemcpyDeviceToHost);
   cudaMemcpy(derivatives_h,derivatives_d,size_d,cudaMemcpyDeviceToHost);
   cudaMemcpy(polvalues_h,polvalues_d,size_pols,cudaMemcpyDeviceToHost);
}

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<double> *xval, complexD<double> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   // up to 32 variables, each in up to degree 16
   __shared__ complexD<double> degrees[14][32];
   // precomputing degrees
   // first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<double> a;
      a = xval[k];
      complexD<double> b(1.0,0.0);
      degrees[0][k] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[j][k] = b;
      }
   }
   __syncthreads(); 
   complexD<double> factor_i(1.0,0.0);
   int pos,exp;
   for(int ind=0; ind<NV; ind++)
   {
      pos = (int)positions[NV*i+ind];
      exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[exp][pos];
   }
   factors[i] = factor_i;
}

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<double> *xval, complexD<double> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   //up to 32 variables, each in up to degree 16
   __shared__ complexD<double> degrees[32][14];
   //precomputing degrees
   //first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<double> a;
      a = xval[k];
      complexD<double> b(1.0,0.0);
      degrees[k][0] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[k][j] = b;
      }
   }
   __syncthreads();
   complexD<double> factor_i(1.0,0.0);
   for(int ind=0; ind<NV; ind++)
   {
      int pos = (int)positions[NV*i+ind];
      int exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[pos][exp];
   }
   factors[i] = factor_i;
}

__global__ void sum_monoms
 ( int dim_s, complexD<double> *monvalues, complexD<double> *polvalues,
   int act_n_threads )
{
   int i=blockIdx.x*blockDim.x + threadIdx.x;
   complexD<double> polvalue_regs(0.0,0.0);
   for(int ind=0;ind<dim_s; ind++)
      polvalue_regs=polvalue_regs + monvalues[ind*act_n_threads + i];
   polvalues[i]=polvalue_regs;
}

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<double> *xval, complexD<double> *roots, complexD<double> *coefs,
   complexD<double> *monvalues, complexD<double> *monderivatives )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int j = threadIdx.x;
   int mn = i % m;
   int pn = i/m;
   int sp = ant*mn;
   __shared__ complexD<double> xv_sh[32];
   __shared__ complexD<double> derivatives[14][32];
   int n_download_rounds=dim/blockDim.x;
   for(int ind=0; ind < n_download_rounds; ind++)
      xv_sh[blockDim.x*ind + j]=xval[blockDim.x*ind +j];
   if(j < dim % blockDim.x )
      xv_sh[blockDim.x * n_download_rounds + j]
         = xval[blockDim.x * n_download_rounds + j];
   __syncthreads();
   // calculate forward products
   // put the first variable in the second cell of the array
   int pos;
   pos = (int)positions[nvarm*i];
   derivatives[1][j]=xv_sh[pos];
   // sequentially compute and store forward products 
   // in the last (# of vars - 2)  cells of the array
   // those are forward products for the last (# of vars) partial derivatives
   for(int ind=0;ind<(nvarm-2);ind++)
   {
      pos = (int)positions[nvarm*i+ind+1];
      derivatives[ind+2][j] = derivatives[ind+1][j]*xv_sh[pos];
   }
   // compute in registers the current backward product
   complexD<double> curr_back_prod;
   pos = (int) positions[nvarm*i + nvarm -1];
   curr_back_prod=xv_sh[pos];
   derivatives[nvarm-2][j]=derivatives[nvarm-2][j] * curr_back_prod;
   derivatives[0][j].init(1.0,0.0);
   for(int ind=0; ind < nvarm - 2; ind ++)
   {
      pos = (int) positions[nvarm*i +nvarm-2-ind];
      curr_back_prod=curr_back_prod * xv_sh[pos];
      derivatives[nvarm-3-ind][j]
         = derivatives[nvarm-3-ind][j] * curr_back_prod;
   }
   for(int ind=0; ind < nvarm; ind ++)
      derivatives[ind][j]= derivatives[ind][j] * roots[i];
   // Computing the monomial value using the second varibale
   // derivatives[nvarm][j]=derivatives[1][j]*xv_sh[pos];
   pos = (int) positions[nvarm*i];
   derivatives[nvarm][j]=derivatives[0][j]*xv_sh[pos];
   // multiply each monomial by its coefficient
   for(int ind=0 ; ind < nvarm+1; ind ++)
      derivatives[ind][j] = derivatives[ind][j]*coefs[tot_n_mons*ind + i];
   // writing to the global memory of derivatives values
   for(int ind=0; ind < nvarm; ind ++)
   {
      pos = (int)positions[nvarm*i+ind];
      monderivatives[sp+dim*(pos+1)+pn] = derivatives[ind][j];
   }
   monderivatives[sp+pn]= derivatives[nvarm][j];
}

// for double double precision 

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<gdd_real> *xval, complexD<gdd_real> *factors );

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<gdd_real> *xval, complexD<gdd_real> *factors );

__global__ void sum_monoms
 ( int dim_s, complexD<gdd_real> *monvalues, complexD<gdd_real> *polvalues,
   int act_n_threads );

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<gdd_real> *xval, complexD<gdd_real> *roots, complexD<gdd_real> *coefs,
   complexD<gdd_real> *monvalues, complexD<gdd_real> *monderivatives );

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<gdd_real> *x_h, complexD<gdd_real> *c_h,
   complexD<gdd_real> *factors_h, complexD<gdd_real> *polvalues_h )
{
   int nblocks = number_of_blocks(dim,BS);
   complexD<gdd_real> *derivatives_d;
   complexD<gdd_real> *factors_d;
   complexD<gdd_real> *monvalues_d;
   complexD<gdd_real> *polvalues_d;
   complexD<gdd_real> *x_d;
   complexD<gdd_real> *c_d;
   // allocate space for output
   int ant = ((dim*dim+dim)/BS + 1)*BS;
   int aas = ant*m;
   //complexD<gdd_real> *derivatives_h = new complexD<gdd_real>[aas];
   complexD<gdd_real> derivatives_h[aas];
   // cudaMalloc((void**)&derivatives_h,aas);
   for(int i=0; i<aas; i++)
      derivatives_h[i].initH(0.0,0.0);
   size_t size_c = ncoefs*sizeof(complexD<gdd_real>);
   size_t size_d = aas*sizeof(complexD<gdd_real>);
   size_t size_pols = ant*sizeof(complexD<gdd_real>);
   // copy positions and exponents to constant memory
   cudaMemcpyToSymbol(positions,pos,NM*NV*sizeof(char));
   cudaMemcpyToSymbol(exponents,exp,NM*NV*sizeof(char));
   size_t size_m = NM*sizeof(complexD<gdd_real>);
   cudaMalloc((void**)&monvalues_d,size_m);

   size_t size = dim*sizeof(complexD<gdd_real>);
   cudaMalloc((void**)&x_d,size);
   cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
   size_t size_NM = NM*sizeof(complexD<gdd_real>);
   cudaMalloc((void**)&factors_d,size_NM);
   cudaMalloc((void**)&c_d,size_c);
   cudaMemcpy(c_d,c_h,size_c,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&derivatives_d,size_d);
   cudaMalloc((void**)&polvalues_d,size_pols);
   cudaMemcpy(derivatives_d,derivatives_h,size_d,cudaMemcpyHostToDevice);
   for(int j=0; j<r; j++)
   {
       mult1_sw_ind_for_shar<<<NM/BS,BS>>>(dim,deg,NV,x_d,factors_d);
       speeldif<<<NM/BS,BS>>>(dim,NM,ant,m,NV,x_d,factors_d,c_d,
                              monvalues_d,derivatives_d);
       sum_monoms<<<nblocks,BS>>>(m,derivatives_d,polvalues_d,ant);
   }
   cudaMemcpy(factors_h,factors_d,size_NM,cudaMemcpyDeviceToHost);
   cudaMemcpy(derivatives_h,derivatives_d,size_d,cudaMemcpyDeviceToHost);
   cudaMemcpy(polvalues_h,polvalues_d,size_pols,cudaMemcpyDeviceToHost);
}

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<gdd_real> *xval, complexD<gdd_real> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   // up to 32 variables, each in up to degree 16
   __shared__ complexD<gdd_real> degrees[14][32];
   // precomputing degrees
   // first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<gdd_real> a;
      a = xval[k];
      complexD<gdd_real> b(1.0,0.0);
      degrees[0][k] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[j][k] = b;
      }
   }
   __syncthreads(); 
   complexD<gdd_real> factor_i(1.0,0.0);
   int pos,exp;
   for(int ind=0; ind<NV; ind++)
   {
      pos = (int)positions[NV*i+ind];
      exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[exp][pos];
   }
   factors[i] = factor_i;
}

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<gdd_real> *xval, complexD<gdd_real> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   //up to 32 variables, each in up to degree 16
   __shared__ complexD<gdd_real> degrees[32][14];
   //precomputing degrees
   //first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<gdd_real> a;
      a = xval[k];
      complexD<gdd_real> b(1.0,0.0);
      degrees[k][0] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[k][j] = b;
      }
   }
   __syncthreads();
   complexD<gdd_real> factor_i(1.0,0.0);
   for(int ind=0; ind<NV; ind++)
   {
      int pos = (int)positions[NV*i+ind];
      int exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[pos][exp];
   }
   factors[i] = factor_i;
}

__global__ void sum_monoms
 ( int dim_s, complexD<gdd_real> *monvalues, complexD<gdd_real> *polvalues,
   int act_n_threads )
{
   int i=blockIdx.x*blockDim.x + threadIdx.x;
   complexD<gdd_real> polvalue_regs(0.0,0.0);
   for(int ind=0;ind<dim_s; ind++)
      polvalue_regs=polvalue_regs + monvalues[ind*act_n_threads + i];
   polvalues[i]=polvalue_regs;
}

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<gdd_real> *xval, complexD<gdd_real> *roots,
   complexD<gdd_real> *coefs,
   complexD<gdd_real> *monvalues, complexD<gdd_real> *monderivatives )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int j = threadIdx.x;
   int mn = i % m;
   int pn = i/m;
   int sp = ant*mn;
   __shared__ complexD<gdd_real> xv_sh[32];
   __shared__ complexD<gdd_real> derivatives[14][32];
   int n_download_rounds=dim/blockDim.x;
   for(int ind=0; ind < n_download_rounds; ind++)
      xv_sh[blockDim.x*ind + j]=xval[blockDim.x*ind +j];
   if(j < dim % blockDim.x )
      xv_sh[blockDim.x * n_download_rounds + j]
         = xval[blockDim.x * n_download_rounds + j];
   __syncthreads();
   // calculate forward products
   // put the first variable in the second cell of the array
   int pos;
   pos = (int)positions[nvarm*i];
   derivatives[1][j]=xv_sh[pos];
   // sequentially compute and store forward products 
   // in the last (# of vars - 2)  cells of the array
   // those are forward products for the last (# of vars) partial derivatives
   for(int ind=0;ind<(nvarm-2);ind++)
   {
      pos = (int)positions[nvarm*i+ind+1];
      derivatives[ind+2][j] = derivatives[ind+1][j]*xv_sh[pos];
   }
   // compute in registers the current backward product
   complexD<gdd_real> curr_back_prod;
   pos = (int) positions[nvarm*i + nvarm -1];
   curr_back_prod=xv_sh[pos];
   derivatives[nvarm-2][j]=derivatives[nvarm-2][j] * curr_back_prod;
   derivatives[0][j].init(1.0,0.0);
   for(int ind=0; ind < nvarm - 2; ind ++)
   {
      pos = (int) positions[nvarm*i +nvarm-2-ind];
      curr_back_prod=curr_back_prod * xv_sh[pos];
      derivatives[nvarm-3-ind][j]
         = derivatives[nvarm-3-ind][j] * curr_back_prod;
   }
   for(int ind=0; ind < nvarm; ind ++)
      derivatives[ind][j]= derivatives[ind][j] * roots[i];
   // Computing the monomial value using the second varibale
   // derivatives[nvarm][j]=derivatives[1][j]*xv_sh[pos];
   pos = (int) positions[nvarm*i];
   derivatives[nvarm][j]=derivatives[0][j]*xv_sh[pos];
   // multiply each monomial by its coefficient
   for(int ind=0 ; ind < nvarm+1; ind ++)
      derivatives[ind][j] = derivatives[ind][j]*coefs[tot_n_mons*ind + i];
   // writing to the global memory of derivatives values
   for(int ind=0; ind < nvarm; ind ++)
   {
      pos = (int)positions[nvarm*i+ind];
      monderivatives[sp+dim*(pos+1)+pn] = derivatives[ind][j];
   }
   monderivatives[sp+pn]= derivatives[nvarm][j];
}

// for quad double precision

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<gqd_real> *xval, complexD<gqd_real> *factors );

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<gqd_real> *xval, complexD<gqd_real> *factors );

__global__ void sum_monoms
 ( int dim_s, complexD<gqd_real> *monvalues, complexD<gqd_real> *polvalues,
   int act_n_threads );

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<gqd_real> *xval, complexD<gqd_real> *roots,
   complexD<gqd_real> *coefs,
   complexD<gqd_real> *monvalues, complexD<gqd_real> *monderivatives );

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<gqd_real> *x_h, complexD<gqd_real> *c_h,
   complexD<gqd_real> *factors_h, complexD<gqd_real> *polvalues_h )
{
   int nblocks = number_of_blocks(dim,BS);
   complexD<gqd_real> *derivatives_d;
   complexD<gqd_real> *factors_d;
   complexD<gqd_real> *monvalues_d;
   complexD<gqd_real> *polvalues_d;
   complexD<gqd_real> *x_d;
   complexD<gqd_real> *c_d;
   // allocate space for output
   int ant = ((dim*dim+dim)/BS + 1)*BS;
   int aas = ant*m;
   // complexD<gqd_real> *derivatives_h = new complexD<gqd_real>[aas]; // illegal
   complexD<gqd_real> derivatives_h[aas]; // replaces the above allocation
   for(int i=0; i<aas; i++)
      derivatives_h[i].initH(0.0,0.0);
   size_t size_c = ncoefs*sizeof(complexD<gqd_real>);
   size_t size_d = aas*sizeof(complexD<gqd_real>);
   size_t size_pols = ant*sizeof(complexD<gqd_real>);
   // copy positions and exponents to constant memory
   cudaMemcpyToSymbol(positions,pos,NM*NV*sizeof(char));
   cudaMemcpyToSymbol(exponents,exp,NM*NV*sizeof(char));
   size_t size_m = NM*sizeof(complexD<gqd_real>);
   cudaMalloc((void**)&monvalues_d,size_m);

   size_t size = dim*sizeof(complexD<gqd_real>);
   cudaMalloc((void**)&x_d,size);
   cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
   size_t size_NM = NM*sizeof(complexD<gqd_real>);
   cudaMalloc((void**)&factors_d,size_NM);
   cudaMalloc((void**)&c_d,size_c);
   cudaMemcpy(c_d,c_h,size_c,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&derivatives_d,size_d);
   cudaMalloc((void**)&polvalues_d,size_pols);
   cudaMemcpy(derivatives_d,derivatives_h,size_d,cudaMemcpyHostToDevice);
   for(int j=0; j<r; j++)
   {
      mult1_sw_ind_for_shar<<<NM/BS,BS>>>(dim,deg,NV,x_d,factors_d);
      speeldif<<<NM/BS,BS>>>(dim,NM,ant,m,NV,x_d,factors_d,c_d,
                             monvalues_d,derivatives_d);
      sum_monoms<<<nblocks,BS>>>(m,derivatives_d,polvalues_d,ant);
   }
   cudaMemcpy(factors_h,factors_d,size_NM,cudaMemcpyDeviceToHost);
   cudaMemcpy(derivatives_h,derivatives_d,size_d,cudaMemcpyDeviceToHost);
   cudaMemcpy(polvalues_h,polvalues_d,size_pols,cudaMemcpyDeviceToHost);
}

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<gqd_real> *xval, complexD<gqd_real> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   // up to 32 variables, each in up to degree 4 ?
   __shared__ complexD<gqd_real> degrees[4][32];
   // precomputing degrees
   // first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<gqd_real> a;
      a = xval[k];
      complexD<gqd_real> b(1.0,0.0);
      degrees[0][k] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[j][k] = b;
      }
   }
   __syncthreads(); 
   complexD<gqd_real> factor_i(1.0,0.0);
   int pos,exp;
   for(int ind=0; ind<NV; ind++)
   {
      pos = (int)positions[NV*i+ind];
      exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[exp][pos];
   }
   factors[i] = factor_i;
}

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<gqd_real> *xval, complexD<gqd_real> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   //up to 32 variables, each in up to degree 16
   __shared__ complexD<gqd_real> degrees[32][4];
   //precomputing degrees
   //first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<gqd_real> a;
      a = xval[k];
      complexD<gqd_real> b(1.0,0.0);
      degrees[k][0] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[k][j] = b;
      }
   }
   __syncthreads();
   complexD<gqd_real> factor_i(1.0,0.0);
   for(int ind=0; ind<NV; ind++)
   {
      int pos = (int)positions[NV*i+ind];
      int exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[pos][exp];
   }
   factors[i] = factor_i;
}

__global__ void sum_monoms
 ( int dim_s, complexD<gqd_real> *monvalues, complexD<gqd_real> *polvalues,
   int act_n_threads)
{
   int i=blockIdx.x*blockDim.x + threadIdx.x;
   complexD<gqd_real> polvalue_regs(0.0,0.0);
   for(int ind=0;ind<dim_s; ind++)
      polvalue_regs=polvalue_regs + monvalues[ind*act_n_threads + i];
   polvalues[i]=polvalue_regs;
}

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<gqd_real> *xval, complexD<gqd_real> *roots, complexD<gqd_real> *coefs,
   complexD<gqd_real> *monvalues, complexD<gqd_real> *monderivatives )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int j = threadIdx.x;
   int mn = i % m;
   int pn = i/m;
   int sp = ant*mn;
   __shared__ complexD<gqd_real>  xv_sh[32];
   __shared__ complexD<gqd_real>  derivatives[4][32];
   int n_download_rounds=dim/blockDim.x;
   for(int ind=0; ind < n_download_rounds; ind++)
      xv_sh[blockDim.x*ind + j]=xval[blockDim.x*ind +j];
   if(j < dim % blockDim.x )
      xv_sh[blockDim.x * n_download_rounds + j]
         = xval[blockDim.x * n_download_rounds + j];
   __syncthreads();
   // calculate forward products
   // put the first variable in the second cell of the array
   int pos;
   pos = (int)positions[nvarm*i];
   derivatives[1][j]=xv_sh[pos];
   // sequentially compute and store forward products 
   // in the last (# of vars - 2)  cells of the array
   // those are forward products for the last (# of vars) partial derivatives
   for(int ind=0;ind<(nvarm-2);ind++)
   {
      pos = (int)positions[nvarm*i+ind+1];
      derivatives[ind+2][j] = derivatives[ind+1][j]*xv_sh[pos];
   }
   // compute in registers the current backward product
   complexD<gqd_real> curr_back_prod;
   pos = (int) positions[nvarm*i + nvarm -1];
   curr_back_prod=xv_sh[pos];
   derivatives[nvarm-2][j]=derivatives[nvarm-2][j] * curr_back_prod;
   derivatives[0][j].init(1.0,0.0);
   for(int ind=0; ind < nvarm - 2; ind ++)
   {
      pos = (int) positions[nvarm*i +nvarm-2-ind];
      curr_back_prod=curr_back_prod * xv_sh[pos];
      derivatives[nvarm-3-ind][j]
         = derivatives[nvarm-3-ind][j] * curr_back_prod;
   }
   for(int ind=0; ind < nvarm; ind ++)
      derivatives[ind][j]= derivatives[ind][j] * roots[i];
   // Computing the monomial value using the second varibale
   // derivatives[nvarm][j]=derivatives[1][j]*xv_sh[pos];
   pos = (int) positions[nvarm*i];
   derivatives[nvarm][j]=derivatives[0][j]*xv_sh[pos];
   // multiply each monomial by its coefficient
   for(int ind=0 ; ind < nvarm+1; ind ++)
      derivatives[ind][j] = derivatives[ind][j]*coefs[tot_n_mons*ind + i];
   // writing to the global memory of derivatives values
   for(int ind=0; ind < nvarm; ind ++)
   {
      pos = (int)positions[nvarm*i+ind];
      monderivatives[sp+dim*(pos+1)+pn] = derivatives[ind][j];
   }
   monderivatives[sp+pn]= derivatives[nvarm][j];
}
