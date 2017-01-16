// Kernels to evaluate and differentiate polynomials in several variables,
// for quad double precision.

#include<iostream>
#include"gqd.cu"
#include"complexD.h"

#define d  0 
#define dd 1
#define qd 2

#ifdef precision
#define p precision
#else
#define p 0
#endif

#if(p == 0)
typedef double realD;
#elif(p == 1)
typedef gdd_real realD;
#else
typedef gqd_real realD;
#endif

using namespace std;

__constant__ char positions[20000];
__constant__ char exponents[20000];

__global__ void mult1_sw_ind_for_shar
 ( int dim, int Mdegr, int NV,
   complexD<realD> *xval, complexD<realD> *factors );

__global__ void mult1
 ( int dim, int Mdegr, int NV,
   complexD<realD> *xval, complexD<realD> *factors );

__global__ void sum_monoms
 ( int dim_s, complexD<realD> *monvalues, complexD<realD> *polvalues,
   int act_n_threads );

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<realD> *xval, complexD<realD> *roots, complexD<realD> *coefs,
   complexD<realD> *monvalues, complexD<realD> *monderivatives );

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

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<realD> *x_h, complexD<realD> *c_h,
   complexD<realD> *factors_h, complexD<realD> *polvalues_h )
{
   int nblocks = number_of_blocks(dim,BS);
   complexD<realD> *derivatives_d;
   complexD<realD> *factors_d;
   complexD<realD> *monvalues_d;
   complexD<realD> *polvalues_d;
   complexD<realD> *x_d;
   complexD<realD> *c_d;
   // allocate space for output
   int ant = ((dim*dim+dim)/BS + 1)*BS;
   int aas = ant*m;
   // complexD<realD> *derivatives_h = new complexD<realD>[aas]; // illegal
   complexD<realD> derivatives_h[aas]; // replaces the above allocation
   for(int i=0; i<aas; i++)
      derivatives_h[i].initH(0.0,0.0);
   size_t size_c = ncoefs*sizeof(complexD<realD>);
   size_t size_d = aas*sizeof(complexD<realD>);
   size_t size_pols = ant*sizeof(complexD<realD>);
   // copy positions and exponents to constant memory
   cudaMemcpyToSymbol(positions,pos,NM*NV*sizeof(char));
   cudaMemcpyToSymbol(exponents,exp,NM*NV*sizeof(char));
   size_t size_m = NM*sizeof(complexD<realD>);
   cudaMalloc((void**)&monvalues_d,size_m);

   size_t size = dim*sizeof(complexD<realD>);
   cudaMalloc((void**)&x_d,size);
   cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
   size_t size_NM = NM*sizeof(complexD<realD>);
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
   complexD<realD> *xval, complexD<realD> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   // up to 32 variables, each in up to degree 4 ?
   __shared__ complexD<realD> degrees[4][32];
   // precomputing degrees
   // first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<realD> a;
      a = xval[k];
      complexD<realD> b(1.0,0.0);
      degrees[0][k] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[j][k] = b;
      }
   }
   __syncthreads(); 
   complexD<realD> factor_i(1.0,0.0);
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
   complexD<realD> *xval, complexD<realD> *factors )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int k = threadIdx.x;
   //up to 32 variables, each in up to degree 16
   __shared__ complexD<realD> degrees[32][4];
   //precomputing degrees
   //first N threads in each block compute degrees of variables
   if(k<dim)
   {
      complexD<realD> a;
      a = xval[k];
      complexD<realD> b(1.0,0.0);
      degrees[k][0] = b;
      for(int j=1; j<Mdegr; j++)
      {
         b = b*a;
         degrees[k][j] = b;
      }
   }
   __syncthreads();
   complexD<realD> factor_i(1.0,0.0);
   for(int ind=0; ind<NV; ind++)
   {
      int pos = (int)positions[NV*i+ind];
      int exp = (int)exponents[NV*i+ind];
      factor_i = factor_i * degrees[pos][exp];
   }
   factors[i] = factor_i;
}

__global__ void sum_monoms
 ( int dim_s, complexD<realD> *monvalues, complexD<realD> *polvalues,
   int act_n_threads)
{
   int i=blockIdx.x*blockDim.x + threadIdx.x;
   complexD<realD> polvalue_regs(0.0,0.0);
   for(int ind=0;ind<dim_s; ind++)
      polvalue_regs=polvalue_regs + monvalues[ind*act_n_threads + i];
   polvalues[i]=polvalue_regs;
}

__global__ void speeldif
 ( int dim, int tot_n_mons, int ant, int m, int nvarm, 
   complexD<realD> *xval, complexD<realD> *roots, complexD<realD> *coefs,
   complexD<realD> *monvalues, complexD<realD> *monderivatives )
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int j = threadIdx.x;
   int mn = i % m;
   int pn = i/m;
   int sp = ant*mn;
   __shared__ complexD<realD>  xv_sh[32];
   __shared__ complexD<realD>  derivatives[4][32];
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
   complexD<realD> curr_back_prod;
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
