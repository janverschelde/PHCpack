/* The file dbl2_tabs_kernels.cu defines the functions specified in
 * the file dbl2_tabs_kernels.h. */

#include <iostream>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_tabs_kernels.h"
#include "dbl_tabs_flopcounts.h"

using namespace std;

__global__ void dbl2_small_invert_upper 
( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolhi[tabsdd_shmemsize];
   __shared__ double Ucollo[tabsdd_shmemsize];
   __shared__ double invUrowshi[tabsdd_shmemsize];
   __shared__ double invUrowslo[tabsdd_shmemsize];

   double rhshi,rhslo,xvalhi,xvallo,acchi,acclo;

   int colidx = dim*(dim-1);          // start with the last column

   Ucolhi[k] = Uhi[colidx+k];         // load the last column
   Ucollo[k] = Ulo[colidx+k];
   rhshi = ((double) int(k == dim-1));  // right hand side for each thread
   rhslo = 0.0;
   int rowidx = (dim - 1)*dim + k;      // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k]; // last row of the inverse
   ddg_div(rhshi,rhslo,Ucolhi[k],Ucollo[k],
           &invUrowshi[rowidx],&invUrowslo[rowidx]);

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolhi[k] = Uhi[colidx+k];
         Ucollo[k] = Ulo[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalhi = invUrowshi[rowidx];
         xvallo = invUrowslo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         ddg_mul(Ucolhi[i],Ucollo[i],xvalhi,xvallo,&acchi,&acclo);
         ddg_dec(&rhshi,&rhslo,acchi,acclo);
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolhi[k] = Uhi[colidx+k];
      Ucollo[k] = Ulo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      ddg_div(rhshi,rhslo,Ucolhi[i],Ucollo[i],
              &invUrowshi[rowidx],&invUrowslo[rowidx]);
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUhi[rowidx+k] = invUrowshi[rowidx+k];
      invUlo[rowidx+k] = invUrowslo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void cmplx2_small_invert_upper
 ( int dim, double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolrehi[tabsdd_shmemsize];
   __shared__ double Ucolrelo[tabsdd_shmemsize];
   __shared__ double Ucolimhi[tabsdd_shmemsize];
   __shared__ double Ucolimlo[tabsdd_shmemsize];
   __shared__ double invUrowsrehi[tabsdd_shmemsize];
   __shared__ double invUrowsrelo[tabsdd_shmemsize];
   __shared__ double invUrowsimhi[tabsdd_shmemsize];
   __shared__ double invUrowsimlo[tabsdd_shmemsize];

   double rhsrehi,rhsrelo,rhsimhi,rhsimlo;
   double xvalrehi,xvalrelo,xvalimhi,xvalimlo;
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;
   double invrehi,invrelo,invimhi,invimlo,denhi,denlo;

   int colidx = dim*(dim-1);             // start with the last column

   Ucolrehi[k] = Urehi[colidx+k];        // load the last column
   Ucolrelo[k] = Urelo[colidx+k];
   Ucolimhi[k] = Uimhi[colidx+k];
   Ucolimlo[k] = Uimlo[colidx+k];
   rhsrehi = ((double) int(k == dim-1)); // right hand side for each thread
   rhsrelo = 0.0;
   rhsimhi = 0.0;
   rhsimlo = 0.0;
   int rowidx = (dim - 1)*dim + k;       // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse
   ddg_mul(Ucolrehi[k],Ucolrelo[k],Ucolrehi[k],Ucolrelo[k],&denhi,&denlo);
   ddg_mul(Ucolimhi[k],Ucolimlo[k],Ucolimhi[k],Ucolimlo[k],&acc1hi,&acc1lo);
   ddg_inc(&denhi,&denlo,acc1hi,acc1lo);
   ddg_div(Ucolrehi[k],Ucolrelo[k],denhi,denlo,&invrehi,&invrelo);
   ddg_div(Ucolimhi[k],Ucolimlo[k],denhi,denlo,&invimhi,&invimlo);
   ddg_minus(&invimhi,&invimlo);
   ddg_mul(rhsrehi,rhsrelo,invrehi,invrelo,&acc1hi,&acc1lo);
   ddg_mul(rhsimhi,rhsimlo,invimhi,invimlo,&acc2hi,&acc2lo);
   ddg_mul(rhsimhi,rhsimlo,invrehi,invrelo,&acc3hi,&acc3lo);
   ddg_mul(rhsrehi,rhsrelo,invimhi,invimlo,&acc4hi,&acc4lo);
   ddg_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
   invUrowsrehi[rowidx] = acc1hi;
   invUrowsrelo[rowidx] = acc1lo;
   ddg_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
   invUrowsimhi[rowidx] = acc3hi;
   invUrowsimlo[rowidx] = acc3lo;

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhsrelo = 0.0;
      rhsimhi = 0.0;
      rhsimlo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolrehi[k] = Urehi[colidx+k];
         Ucolrelo[k] = Urelo[colidx+k];
         Ucolimhi[k] = Uimhi[colidx+k];
         Ucolimlo[k] = Uimlo[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalrehi = invUrowsrehi[rowidx];
         xvalrelo = invUrowsrelo[rowidx];
         xvalimhi = invUrowsimhi[rowidx];
         xvalimlo = invUrowsimlo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         ddg_mul(Ucolrehi[i],Ucolrelo[i],xvalrehi,xvalrelo,&acc1hi,&acc1lo);
         ddg_mul(Ucolimhi[i],Ucolimlo[i],xvalimhi,xvalimlo,&acc2hi,&acc2lo);
         ddg_mul(Ucolimhi[i],Ucolimlo[i],xvalrehi,xvalrelo,&acc3hi,&acc3lo);
         ddg_mul(Ucolrehi[i],Ucolrelo[i],xvalimhi,xvalimlo,&acc4hi,&acc4lo);
         ddg_dec(&rhsrehi,&rhsrelo,acc1hi,acc1lo);
         ddg_inc(&rhsrehi,&rhsrelo,acc2hi,acc2lo);
         ddg_dec(&rhsimhi,&rhsimlo,acc3hi,acc3lo);
         ddg_dec(&rhsimhi,&rhsimlo,acc4hi,acc4lo);
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolrehi[k] = Urehi[colidx+k];
      Ucolrelo[k] = Urelo[colidx+k];
      Ucolimhi[k] = Uimhi[colidx+k];
      Ucolimlo[k] = Uimlo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      ddg_mul(Ucolrehi[i],Ucolrelo[i],Ucolrehi[i],Ucolrelo[i],&denhi,&denlo);
      ddg_mul(Ucolimhi[i],Ucolimlo[i],Ucolimhi[i],Ucolimlo[i],&acc1hi,&acc1lo);
      __syncthreads();
      ddg_inc(&denhi,&denlo,acc1hi,acc1lo);
      ddg_div(Ucolrehi[i],Ucolrelo[i],denhi,denlo,&invrehi,&invrelo);
      ddg_div(Ucolimhi[i],Ucolimlo[i],denhi,denlo,&invimhi,&invimlo);
      ddg_minus(&invimhi,&invimlo);
      ddg_mul(rhsrehi,rhsrelo,invrehi,invrelo,&acc1hi,&acc1lo);
      ddg_mul(rhsimhi,rhsimlo,invimhi,invimlo,&acc2hi,&acc2lo);
      ddg_mul(rhsimhi,rhsimlo,invrehi,invrelo,&acc3hi,&acc3lo);
      ddg_mul(rhsrehi,rhsrelo,invimhi,invimlo,&acc4hi,&acc4lo);
      __syncthreads();
      ddg_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
      invUrowsrehi[rowidx] = acc1hi;
      invUrowsrelo[rowidx] = acc1lo;
      ddg_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
      invUrowsimhi[rowidx] = acc3hi;
      invUrowsimlo[rowidx] = acc3lo;
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUrehi[rowidx+k] = invUrowsrehi[rowidx+k];
      invUrelo[rowidx+k] = invUrowsrelo[rowidx+k];
      invUimhi[rowidx+k] = invUrowsimhi[rowidx+k];
      invUimlo[rowidx+k] = invUrowsimlo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void dbl2_medium_invert_upper
 ( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo)
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolhi[tabsdd_shmemsize];      // one column of U
   __shared__ double Ucollo[tabsdd_shmemsize];      // one column of U
   __shared__ double invUrowhi[tabsdd_shmemsize];   // one row of invU
   __shared__ double invUrowlo[tabsdd_shmemsize];   // one row of invU

   double rhshi,rhslo,xvalhi,xvallo,acchi,acclo;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolhi[k] = Uhi[colidx+k];          // load the last column
   Ucollo[k] = Ulo[colidx+k];
   rhshi = ((double) int(k == dim-1)); // right hand side for each thread
   rhslo = 0.0;
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   ddg_div(rhshi,rhslo,Ucolhi[k],Ucollo[k],&invUrowhi[k],&invUrowlo[k]);
   invUhi[rowidx] = invUrowhi[k];     // store the last row into invU
   invUlo[rowidx] = invUrowlo[k]; 

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucolhi[k] = Uhi[colidx+k];
         Ucollo[k] = Ulo[colidx+k];

         rowidx = j*dim + k;            // need solution value
         invUrowhi[k] = invUhi[rowidx]; // load invU row into invUrow
         invUrowlo[k] = invUlo[rowidx];
         xvalhi = invUrowhi[k];
         xvallo = invUrowlo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         ddg_mul(Ucolhi[i],Ucollo[i],xvalhi,xvallo,&acchi,&acclo);
         ddg_dec(&rhshi,&rhslo,acchi,acclo);
      }
      colidx = dim*i;                 // need column i of U
      Ucolhi[k] = Uhi[colidx+k];
      Ucollo[k] = Ulo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      ddg_div(rhshi,rhslo,Ucolhi[i],Ucollo[i],&invUrowhi[k],&invUrowlo[k]);
      invUhi[rowidx] = invUrowhi[k];
      invUlo[rowidx] = invUrowlo[k];
   }
}

__global__ void cmplx2_medium_invert_upper
 ( int dim, double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolrehi[tabsdd_shmemsize];    // one column of U
   __shared__ double Ucolrelo[tabsdd_shmemsize]; 
   __shared__ double Ucolimhi[tabsdd_shmemsize];
   __shared__ double Ucolimlo[tabsdd_shmemsize]; 
   __shared__ double invUrowrehi[tabsdd_shmemsize]; // one row of invU
   __shared__ double invUrowrelo[tabsdd_shmemsize]; 
   __shared__ double invUrowimhi[tabsdd_shmemsize]; 
   __shared__ double invUrowimlo[tabsdd_shmemsize]; 

   double rhsrehi,rhsrelo,rhsimhi,rhsimlo;
   double xvalrehi,xvalrelo,xvalimhi,xvalimlo;
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;
   double invrehi,invrelo,invimhi,invimlo,denhi,denlo;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolrehi[k] = Urehi[colidx+k];      // load the last column
   Ucolrelo[k] = Urelo[colidx+k];
   Ucolimhi[k] = Uimhi[colidx+k];
   Ucolimlo[k] = Uimlo[colidx+k];
   rhsrehi = ((double) int(k == dim-1)); // right hand side for each thread
   rhsrelo = 0.0;
   rhsimhi = 0.0;
   rhsimlo = 0.0;
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   __syncthreads();
   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   ddg_mul(Ucolrehi[k],Ucolrelo[k],Ucolrehi[k],Ucolrelo[k],&denhi,&denlo);
   ddg_mul(Ucolimhi[k],Ucolimlo[k],Ucolimhi[k],Ucolimlo[k],&acc1hi,&acc1lo);
   ddg_inc(&denhi,&denlo,acc1hi,acc1lo);
   ddg_div(Ucolrehi[k],Ucolrelo[k],denhi,denlo,&invrehi,&invrelo);
   ddg_div(Ucolimhi[k],Ucolimlo[k],denhi,denlo,&invimhi,&invimlo);
   ddg_minus(&invimhi,&invimlo);
   ddg_mul(rhsrehi,rhsrelo,invrehi,invrelo,&acc1hi,&acc1lo);
   ddg_mul(rhsimhi,rhsimlo,invimhi,invimlo,&acc2hi,&acc2lo);
   ddg_mul(rhsimhi,rhsimlo,invrehi,invrelo,&acc3hi,&acc3lo);
   ddg_mul(rhsrehi,rhsrelo,invimhi,invimlo,&acc4hi,&acc4lo);
   ddg_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
   invUrowrehi[k] = acc1hi;
   invUrowrelo[k] = acc1lo;
   ddg_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
   invUrowimhi[k] = acc3hi;
   invUrowimlo[k] = acc3lo;
   invUrehi[rowidx] = invUrowrehi[k];     // store the last row into invU
   invUrelo[rowidx] = invUrowrelo[k]; 
   invUimhi[rowidx] = invUrowimhi[k];
   invUimlo[rowidx] = invUrowimlo[k]; 

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhsrelo = 0.0;
      rhsimhi = 0.0;
      rhsimlo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucolrehi[k] = Urehi[colidx+k];
         Ucolrelo[k] = Urelo[colidx+k];
         Ucolimhi[k] = Uimhi[colidx+k];
         Ucolimlo[k] = Uimlo[colidx+k];

         rowidx = j*dim + k;            // need solution value
         invUrowrehi[k] = invUrehi[rowidx]; // load invU row into invUrow
         invUrowrelo[k] = invUrelo[rowidx];
         invUrowimhi[k] = invUimhi[rowidx];
         invUrowimlo[k] = invUimlo[rowidx];
         xvalrehi = invUrowrehi[k];
         xvalrelo = invUrowrelo[k];
         xvalimhi = invUrowimhi[k];
         xvalimlo = invUrowimlo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         ddg_mul(Ucolrehi[i],Ucolrelo[i],xvalrehi,xvalrelo,&acc1hi,&acc1lo);
         ddg_mul(Ucolimhi[i],Ucolimlo[i],xvalimhi,xvalimlo,&acc2hi,&acc2lo);
         ddg_mul(Ucolimhi[i],Ucolimlo[i],xvalrehi,xvalrelo,&acc3hi,&acc3lo);
         ddg_mul(Ucolrehi[i],Ucolrelo[i],xvalimhi,xvalimlo,&acc4hi,&acc4lo);
         ddg_dec(&rhsrehi,&rhsrelo,acc1hi,acc1lo);
         ddg_inc(&rhsrehi,&rhsrelo,acc2hi,acc2lo);
         ddg_dec(&rhsimhi,&rhsimlo,acc3hi,acc3lo);
         ddg_dec(&rhsimhi,&rhsimlo,acc4hi,acc4lo);
      }
      colidx = dim*i;                 // need column i of U
      Ucolrehi[k] = Urehi[colidx+k];
      Ucolrelo[k] = Urelo[colidx+k];
      Ucolimhi[k] = Uimhi[colidx+k];
      Ucolimlo[k] = Uimlo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      ddg_mul(Ucolrehi[i],Ucolrelo[i],Ucolrehi[i],Ucolrelo[i],&denhi,&denlo);
      ddg_mul(Ucolimhi[i],Ucolimlo[i],Ucolimhi[i],Ucolimlo[i],&acc1hi,&acc1lo);
      ddg_inc(&denhi,&denlo,acc1hi,acc1lo);
      ddg_div(Ucolrehi[i],Ucolrelo[i],denhi,denlo,&invrehi,&invrelo);
      ddg_div(Ucolimhi[i],Ucolimlo[i],denhi,denlo,&invimhi,&invimlo);
      ddg_minus(&invimhi,&invimlo);
      ddg_mul(rhsrehi,rhsrelo,invrehi,invrelo,&acc1hi,&acc1lo);
      ddg_mul(rhsimhi,rhsimlo,invimhi,invimlo,&acc2hi,&acc2lo);
      ddg_mul(rhsimhi,rhsimlo,invrehi,invrelo,&acc3hi,&acc3lo);
      ddg_mul(rhsrehi,rhsrelo,invimhi,invimlo,&acc4hi,&acc4lo);
      ddg_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
      invUrowrehi[k] = acc1hi;
      invUrowrelo[k] = acc1lo;
      ddg_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
      invUrowimhi[k] = acc3hi;
      invUrowimlo[k] = acc3lo;
      invUrehi[rowidx] = invUrowrehi[k];
      invUrelo[rowidx] = invUrowrelo[k];
      invUimhi[rowidx] = invUrowimhi[k];
      invUimlo[rowidx] = invUrowimlo[k];
   }
}

__global__ void  dbl2_invert_tiles
 ( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolhi[tabsdd_shmemsize];      // one column of U
   __shared__ double Ucollo[tabsdd_shmemsize];
   __shared__ double invUrowhi[tabsdd_shmemsize];   // one row of invU
   __shared__ double invUrowlo[tabsdd_shmemsize]; 

   double rhshi,rhslo,xvalhi,xvallo,acchi,acclo;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolhi[k] = Uhi[colidx+k];         // load the last column
   Ucollo[k] = Ulo[colidx+k];
   rhshi = ((double) int(k == dim-1));  // right hand side for each thread
   rhslo = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
   invUhi[rowidx] = 0.0;      // initialize in case of zero divisor
   invUlo[rowidx] = 0.0;
   if(1.0 + Ucolhi[k] != 1.0)
   {
      ddg_div(rhshi,rhslo,Ucolhi[k],Ucollo[k],&invUrowhi[k],&invUrowlo[k]);
      invUhi[rowidx] = invUrowhi[k];     // store the last row into invU
      invUlo[rowidx] = invUrowlo[k];
   }
   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshi = ((double) int(k == i));   // set rhs for i-th unit vector
      rhslo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;     // need column j of U
         Ucolhi[k] = Uhi[colidx+k];
         Ucollo[k] = Ulo[colidx+k];

         rowidx = offset + j*dim + k; // need solution value
         invUrowhi[k] = invUhi[rowidx]; // load invU row into invUrow
         invUrowlo[k] = invUlo[rowidx]; // load invU row into invUrow
         xvalhi = invUrowhi[k];
         xvallo = invUrowlo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         ddg_mul(Ucolhi[i],Ucollo[i],xvalhi,xvallo,&acchi,&acclo);
         ddg_dec(&rhshi,&rhslo,acchi,acclo);
      }
      colidx = offset + dim*i;        // need column i of U
      Ucolhi[k] = Uhi[colidx+k];
      Ucollo[k] = Ulo[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      invUhi[rowidx] = 0.0;
      invUlo[rowidx] = 0.0;
      if(1.0 + Ucolhi[i] != 1.0)
      {
         ddg_div(rhshi,rhslo,Ucolhi[i],Ucollo[i],&invUrowhi[k],&invUrowlo[k]);
         invUhi[rowidx] = invUrowhi[k];
         invUlo[rowidx] = invUrowlo[k];
      }
   }
}

__global__ void  cmplx2_invert_tiles
 ( int dim, double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolrehi[tabsdd_shmemsize];    // one column of U
   __shared__ double Ucolrelo[tabsdd_shmemsize];
   __shared__ double Ucolimhi[tabsdd_shmemsize]; 
   __shared__ double Ucolimlo[tabsdd_shmemsize];
   __shared__ double invUrowrehi[tabsdd_shmemsize];   // one row of invU
   __shared__ double invUrowrelo[tabsdd_shmemsize]; 
   __shared__ double invUrowimhi[tabsdd_shmemsize];
   __shared__ double invUrowimlo[tabsdd_shmemsize]; 

   double rhsrehi,rhsrelo,rhsimhi,rhsimlo;
   double xvalrehi,xvalrelo,xvalimhi,xvalimlo;
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;
   double invrehi,invrelo,invimhi,invimlo,denhi,denlo;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolrehi[k] = Urehi[colidx+k];       // load the last column
   Ucolrelo[k] = Urelo[colidx+k];
   Ucolimhi[k] = Uimhi[colidx+k];
   Ucolimlo[k] = Uimlo[colidx+k];
   rhsrehi = ((double) int(k == dim-1));  // right hand side for each thread
   rhsrelo = 0.0;
   rhsimhi = 0.0;
   rhsimlo = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
   ddg_mul(Ucolrehi[k],Ucolrelo[k],Ucolrehi[k],Ucolrelo[k],&denhi,&denlo);
   ddg_mul(Ucolimhi[k],Ucolimlo[k],Ucolimhi[k],Ucolimlo[k],&acc1hi,&acc1lo);
   ddg_inc(&denhi,&denlo,acc1hi,acc1lo);

   invUrehi[rowidx] = 0.0;  // initialize in case of zero denominator
   invUrelo[rowidx] = 0.0;
   invUimhi[rowidx] = 0.0;
   invUimlo[rowidx] = 0.0;
   
   if(1.0 + denhi != 1.0)
   {
      ddg_div(Ucolrehi[k],Ucolrelo[k],denhi,denlo,&invrehi,&invrelo);
      ddg_div(Ucolimhi[k],Ucolimlo[k],denhi,denlo,&invimhi,&invimlo);

      ddg_minus(&invimhi,&invimlo);
      ddg_mul(rhsrehi,rhsrelo,invrehi,invrelo,&acc1hi,&acc1lo);
      ddg_mul(rhsimhi,rhsimlo,invimhi,invimlo,&acc2hi,&acc2lo);
      ddg_mul(rhsimhi,rhsimlo,invrehi,invrelo,&acc3hi,&acc3lo);
      ddg_mul(rhsrehi,rhsrelo,invimhi,invimlo,&acc4hi,&acc4lo);
      ddg_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
      invUrowrehi[k] = acc1hi;
      invUrowrelo[k] = acc1lo;
      ddg_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
      invUrowimhi[k] = acc3hi;
      invUrowimlo[k] = acc3lo;
      invUrehi[rowidx] = invUrowrehi[k];     // store the last row into invU
      invUrelo[rowidx] = invUrowrelo[k];
      invUimhi[rowidx] = invUrowimhi[k];
      invUimlo[rowidx] = invUrowimlo[k];
   }
   __syncthreads();
   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehi = ((double) int(k == i));   // set rhs for i-th unit vector
      rhsrelo = 0.0;
      rhsimhi = 0.0;
      rhsimlo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;        // need column j of U
         Ucolrehi[k] = Urehi[colidx+k];
         Ucolrelo[k] = Urelo[colidx+k];
         Ucolimhi[k] = Uimhi[colidx+k];
         Ucolimlo[k] = Uimlo[colidx+k];

         rowidx = offset + j*dim + k;       // need solution value
         invUrowrehi[k] = invUrehi[rowidx]; // load invU row into invUrow
         invUrowrelo[k] = invUrelo[rowidx];
         invUrowimhi[k] = invUimhi[rowidx];
         invUrowimlo[k] = invUimlo[rowidx];
         xvalrehi = invUrowrehi[k];
         xvalrelo = invUrowrelo[k];
         xvalimhi = invUrowimhi[k];
         xvalimlo = invUrowimlo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         ddg_mul(Ucolrehi[i],Ucolrelo[i],xvalrehi,xvalrelo,&acc1hi,&acc1lo);
         ddg_mul(Ucolimhi[i],Ucolimlo[i],xvalimhi,xvalimlo,&acc2hi,&acc2lo);
         ddg_mul(Ucolimhi[i],Ucolimlo[i],xvalrehi,xvalrelo,&acc3hi,&acc3lo);
         ddg_mul(Ucolrehi[i],Ucolrelo[i],xvalimhi,xvalimlo,&acc4hi,&acc4lo);
         ddg_dec(&rhsrehi,&rhsrelo,acc1hi,acc1lo);
         ddg_inc(&rhsrehi,&rhsrelo,acc2hi,acc2lo);
         ddg_dec(&rhsimhi,&rhsimlo,acc3hi,acc3lo);
         ddg_dec(&rhsimhi,&rhsimlo,acc4hi,acc4lo);
      }
      colidx = offset + dim*i;        // need column i of U
      Ucolrehi[k] = Urehi[colidx+k];
      Ucolrelo[k] = Urelo[colidx+k];
      Ucolimhi[k] = Uimhi[colidx+k];
      Ucolimlo[k] = Uimlo[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      ddg_mul(Ucolrehi[i],Ucolrelo[i],Ucolrehi[i],Ucolrelo[i],&denhi,&denlo);
      ddg_mul(Ucolimhi[i],Ucolimlo[i],Ucolimhi[i],Ucolimlo[i],&acc1hi,&acc1lo);
      ddg_inc(&denhi,&denlo,acc1hi,acc1lo);

      invUrehi[rowidx] = 0.0; // initialize in case of zero denominator
      invUrelo[rowidx] = 0.0;
      invUimhi[rowidx] = 0.0;
      invUimlo[rowidx] = 0.0;

      if(1.0 + denhi != 1.0)
      {
         ddg_div(Ucolrehi[i],Ucolrelo[i],denhi,denlo,&invrehi,&invrelo);
         ddg_div(Ucolimhi[i],Ucolimlo[i],denhi,denlo,&invimhi,&invimlo);
         ddg_minus(&invimhi,&invimlo);
         ddg_mul(rhsrehi,rhsrelo,invrehi,invrelo,&acc1hi,&acc1lo);
         ddg_mul(rhsimhi,rhsimlo,invimhi,invimlo,&acc2hi,&acc2lo);
         ddg_mul(rhsimhi,rhsimlo,invrehi,invrelo,&acc3hi,&acc3lo);
         ddg_mul(rhsrehi,rhsrelo,invimhi,invimlo,&acc4hi,&acc4lo);
         ddg_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
         invUrowrehi[k] = acc1hi;
         invUrowrelo[k] = acc1lo;
         ddg_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
         invUrowimhi[k] = acc3hi;
         invUrowimlo[k] = acc3lo;
         invUrehi[rowidx] = invUrowrehi[k];
         invUrelo[rowidx] = invUrowrelo[k];
         invUimhi[rowidx] = invUrowimhi[k];
         invUimlo[rowidx] = invUrowimlo[k];
      }
   }
}

__global__ void dbl2_multiply_inverse
 ( int dim, int idx, double *invUhi, double *invUlo,
   double *whi, double *wlo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workhi[tabsdd_shmemsize];      // copy of w
   __shared__ double worklo[tabsdd_shmemsize];      // copy of w

   workhi[k] = whi[rhsoff+k];
   worklo[k] = wlo[rhsoff+k];

   double resulthi = 0.0; // each thread stores its product in result
   double resultlo = 0.0;
   double coeffhi,coefflo,acchi,acclo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhi = invUhi[offset+k*dim+j]; // thread k does row k
      coefflo = invUlo[offset+k*dim+j];
      // result = result + coeff*work[j];
      ddg_mul(coeffhi,coefflo,workhi[j],worklo[j],&acchi,&acclo);
      ddg_inc(&resulthi,&resultlo,acchi,acclo);
   }
   whi[rhsoff+k] = resulthi;
   wlo[rhsoff+k] = resultlo;
}

__global__ void cmplx2_multiply_inverse
 ( int dim, int idx,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo,
   double *wrehi, double *wrelo, double *wimhi, double *wimlo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workrehi[tabsdd_shmemsize];      // copy of w
   __shared__ double workrelo[tabsdd_shmemsize]; 
   __shared__ double workimhi[tabsdd_shmemsize];
   __shared__ double workimlo[tabsdd_shmemsize];

   workrehi[k] = wrehi[rhsoff+k];
   workrelo[k] = wrelo[rhsoff+k];
   workimhi[k] = wimhi[rhsoff+k];
   workimlo[k] = wimlo[rhsoff+k];

   double resultrehi = 0.0; // each thread stores its product in result
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;
   double coeffrehi,coeffrelo,coeffimhi,coeffimlo;
   double acc1hi,acc1lo,acc2hi,acc2lo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffrehi = invUrehi[offset+k*dim+j]; // thread k does row k
      coeffrelo = invUrelo[offset+k*dim+j];
      coeffimhi = invUimhi[offset+k*dim+j];
      coeffimlo = invUimlo[offset+k*dim+j];
      // result = result + coeff*work[j];
      ddg_mul(coeffrehi,coeffrelo,workrehi[j],workrelo[j],&acc1hi,&acc1lo);
      ddg_mul(coeffimhi,coeffimlo,workimhi[j],workimlo[j],&acc2hi,&acc2lo);
      ddg_inc(&resultrehi,&resultrelo,acc1hi,acc1lo);
      ddg_dec(&resultrehi,&resultrelo,acc2hi,acc2lo);
      ddg_mul(coeffimhi,coeffimlo,workrehi[j],workrelo[j],&acc1hi,&acc1lo);
      ddg_mul(coeffrehi,coeffrelo,workimhi[j],workimlo[j],&acc2hi,&acc2lo);
      ddg_inc(&resultimhi,&resultimlo,acc1hi,acc1lo);
      ddg_inc(&resultimhi,&resultimlo,acc2hi,acc2lo);
   }
   wrehi[rhsoff+k] = resultrehi; wrelo[rhsoff+k] = resultrelo;
   wimhi[rhsoff+k] = resultimhi; wimlo[rhsoff+k] = resultimlo;
}

__global__ void dbl2_back_substitute
 ( int dim, int idx, double *Uhi, double *Ulo, double *whi, double *wlo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkhi[tabsdd_shmemsize];   // copy of w
   __shared__ double wrklo[tabsdd_shmemsize]; 
   __shared__ double solhi[tabsdd_shmemsize];    // solution to update with
   __shared__ double sollo[tabsdd_shmemsize];

   wrkhi[k] = whi[B*dim+k];    // block B updates B-th slice of w
   wrklo[k] = wlo[B*dim+k];
   solhi[k] = whi[idx*dim+k];  // solution that is back substituted
   sollo[k] = wlo[idx*dim+k];

   double resulthi = 0.0; // each thread stores its product in result
   double resultlo = 0.0;
   double coeffhi,coefflo,acchi,acclo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhi = Uhi[offset+k*dim+j];
      coefflo = Ulo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      ddg_mul(coeffhi,coefflo,solhi[j],sollo[j],&acchi,&acclo);
      ddg_inc(&resulthi,&resultlo,acchi,acclo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   ddg_dec(&wrkhi[k],&wrklo[k],resulthi,resultlo);
   whi[B*dim+k] = wrkhi[k];
   wlo[B*dim+k] = wrklo[k];
}

__global__ void cmplx2_back_substitute
 ( int dim, int idx,
   double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *wrehi, double *wrelo, double *wimhi, double *wimlo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkrehi[tabsdd_shmemsize];   // copy of w
   __shared__ double wrkrelo[tabsdd_shmemsize]; 
   __shared__ double wrkimhi[tabsdd_shmemsize];
   __shared__ double wrkimlo[tabsdd_shmemsize]; 
   __shared__ double solrehi[tabsdd_shmemsize];    // solution to update with
   __shared__ double solrelo[tabsdd_shmemsize];
   __shared__ double solimhi[tabsdd_shmemsize];
   __shared__ double solimlo[tabsdd_shmemsize];

   wrkrehi[k] = wrehi[B*dim+k];    // block B updates B-th slice of w
   wrkrelo[k] = wrelo[B*dim+k];
   wrkimhi[k] = wimhi[B*dim+k];
   wrkimlo[k] = wimlo[B*dim+k];
   solrehi[k] = wrehi[idx*dim+k];  // solution that is back substituted
   solrelo[k] = wrelo[idx*dim+k];
   solimhi[k] = wimhi[idx*dim+k];
   solimlo[k] = wimlo[idx*dim+k];

   double resultrehi = 0.0; // each thread stores its product in result
   double resultrelo = 0.0;
   double resultimhi = 0.0;
   double resultimlo = 0.0;
   double coeffrehi,coeffrelo,coeffimhi,coeffimlo;
   double acc1hi,acc1lo,acc2hi,acc2lo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffrehi = Urehi[offset+k*dim+j];
      coeffrelo = Urelo[offset+k*dim+j];
      coeffimhi = Uimhi[offset+k*dim+j];
      coeffimlo = Uimlo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      ddg_mul(coeffrehi,coeffrelo,solrehi[j],solrelo[j],&acc1hi,&acc1lo);
      ddg_mul(coeffimhi,coeffimlo,solimhi[j],solimlo[j],&acc2hi,&acc2lo);
      ddg_inc(&resultrehi,&resultrelo,acc1hi,acc1lo);
      ddg_dec(&resultrehi,&resultrelo,acc2hi,acc2lo);
      ddg_mul(coeffimhi,coeffimlo,solrehi[j],solrelo[j],&acc1hi,&acc1lo);
      ddg_mul(coeffrehi,coeffrelo,solimhi[j],solimlo[j],&acc2hi,&acc2lo);
      ddg_inc(&resultimhi,&resultimlo,acc1hi,acc1lo);
      ddg_inc(&resultimhi,&resultimlo,acc2hi,acc2lo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   ddg_dec(&wrkrehi[k],&wrkrelo[k],resultrehi,resultrelo);
   ddg_dec(&wrkimhi[k],&wrkimlo[k],resultimhi,resultimlo);
   wrehi[B*dim+k] = wrkrehi[k];
   wrelo[B*dim+k] = wrkrelo[k];
   wimhi[B*dim+k] = wrkimhi[k];
   wimlo[B*dim+k] = wrkimlo[k];
}

void GPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo,
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Uhi_h = new double[szU];     // Uhi_h stores the columns of Uhi
   double *Ulo_h = new double[szU];     // Ulo_h stores the columns of Ulo 
   double *Uhi_d;                       // Uhi_d is Uhi_h on the device
   double *Ulo_d;                       // Ulo_d is Ulo_h on the device
   double *invUhi_h = new double[szU];  // high doubles of the inverse
   double *invUlo_h = new double[szU];  // low doubles of the inverse
   double *invUhi_d;                    // invUhi_d is invUhi_h on the device
   double *invUlo_d;                    // invUlo_d is invUlo_h on the device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Uhi_h[ix]   = Uhi[i][j];
         Ulo_h[ix++] = Ulo[i][j];
      }

   // only for debugging
   // test_dbl2_small_invert_upper(dim,Uhi_h,Ulo_h,invUhi,invUlo_h);

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Uhi_d,szmat);
   cudaMalloc((void**)&Ulo_d,szmat);
   cudaMalloc((void**)&invUhi_d,szmat);
   cudaMalloc((void**)&invUlo_d,szmat);
   cudaMemcpy(Uhi_d,Uhi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulo_d,Ulo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      dbl2_small_invert_upper<<<1,dim>>>(dim,Uhi_d,Ulo_d,invUhi_d,invUlo_d);
   else
      dbl2_medium_invert_upper<<<1,dim>>>(dim,Uhi_d,Ulo_d,invUhi_d,invUlo_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUhi_h,invUhi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlo_h,invUlo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUhi[i][j] = invUhi_h[ix];
         invUlo[i][j] = invUlo_h[ix++];
      }

   free(Uhi_h); free(invUhi_h);
   free(Ulo_h); free(invUlo_h);
}

void GPU_cmplx2_upper_inverse
 ( int dim, double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double **invUrehi, double **invUrelo, double **invUimhi, double **invUimlo,
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Urehi_h = new double[szU];    // Urehi_h has high real parts
   double *Urelo_h = new double[szU];    // Urelo_h has low real parts
   double *Uimhi_h = new double[szU];    // Uimhi_h has high imag parts
   double *Uimlo_h = new double[szU];    // Uimlo_h has low imag parts
   double *Urehi_d;                      // Urehi_d is Urehi_h on the device
   double *Urelo_d;                      // Urelo_d is Urelo_h on the device
   double *Uimhi_d;                      // Uimhi_d is Uimhi_h on the device
   double *Uimlo_d;                      // Uimlo_d is Uimlo_h on the device
   double *invUrehi_h = new double[szU]; // high real parts of the inverse
   double *invUrelo_h = new double[szU]; // low real parts of the inverse
   double *invUimhi_h = new double[szU]; // high imag parts of the inverse
   double *invUimlo_h = new double[szU]; // low imag parts of the inverse
   double *invUrehi_d;                   // invUrehi_d ~ invUrehi_h on device
   double *invUrelo_d;                   // invUrelo_d ~ invUrelo_h on device
   double *invUimhi_d;                   // invUimhi_d ~ invUimhi_h on device
   double *invUimlo_d;                   // invUimlo_d ~ invUimlo_h on device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Urehi_h[ix] = Urehi[i][j]; Urelo_h[ix]   = Urelo[i][j];
         Uimhi_h[ix] = Uimhi[i][j]; Uimlo_h[ix++] = Uimlo[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Urehi_d,szmat);
   cudaMalloc((void**)&Urelo_d,szmat);
   cudaMalloc((void**)&Uimhi_d,szmat);
   cudaMalloc((void**)&Uimlo_d,szmat);
   cudaMalloc((void**)&invUrehi_d,szmat);
   cudaMalloc((void**)&invUrelo_d,szmat);
   cudaMalloc((void**)&invUimhi_d,szmat);
   cudaMalloc((void**)&invUimlo_d,szmat);
   cudaMemcpy(Urehi_d,Urehi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelo_d,Urelo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhi_d,Uimhi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlo_d,Uimlo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      cmplx2_small_invert_upper<<<1,dim>>>
         (dim,   Urehi_d,   Urelo_d,   Uimhi_d,   Uimlo_d,
              invUrehi_d,invUrelo_d,invUimhi_d,invUimlo_d);
   else
      cmplx2_medium_invert_upper<<<1,dim>>>
         (dim,   Urehi_d,   Urelo_d,   Uimhi_d,   Uimlo_d,
              invUrehi_d,invUrelo_d,invUimhi_d,invUimlo_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUrehi_h,invUrehi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelo_h,invUrelo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhi_h,invUimhi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlo_h,invUimlo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUrehi[i][j] = invUrehi_h[ix];
         invUrelo[i][j] = invUrelo_h[ix];
         invUimhi[i][j] = invUimhi_h[ix];
         invUimlo[i][j] = invUimlo_h[ix++];
      }

   free(Urehi_h); free(Urelo_h); free(invUrehi_h); free(invUrelo_h);
   free(Uimhi_h); free(Uimlo_h); free(invUimhi_h); free(invUimlo_h);
}

void GPU_dbl2_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Uhi, double **Ulo,
   double *bhi, double *blo, double *xhi, double *xlo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *Dhi_h = new double[nbr];    // the diagonal tiles on the host
   double *Dlo_h = new double[nbr];    // low doubles of diagonal tiles
   double *Dhi_d;                      // diagonal tiles on the device
   double *Dlo_d;                      // low doubles of diagonal tiles
   double *invDhi_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invDlo_h = new double[nbr]; // low doubles of inverse tiles
   double *invDhi_d;                   // invDhi_d is invDhi_h on device
   double *invDlo_d;                   // invDlo_d is invDlo_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Dhi_h[ix]   = Uhi[offset+i][offset+j];
            Dlo_h[ix++] = Ulo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Dhi_d,sznum);
   cudaMalloc((void**)&Dlo_d,sznum);
   cudaMalloc((void**)&invDhi_d,sznum);
   cudaMalloc((void**)&invDlo_d,sznum);
   cudaMemcpy(Dhi_d,Dhi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlo_d,Dlo_h,sznum,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *invlapms = 0.0;
   *mullapms = 0.0;
   *sublapms = 0.0;
   *totlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   dbl2_invert_tiles<<<nbt,szt>>>(szt,Dhi_d,Dlo_d,invDhi_d,invDlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhshi_d;                    // right hand side on device
   double *rhslo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhshi_d,szrhs);
   cudaMalloc((void**)&rhslo_d,szrhs);
   cudaMemcpy(rhshi_d,bhi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslo_d,blo,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl2_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDhi_d,invDlo_d,rhshi_d,rhslo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;           // #doubles in column of U
   double *Ucolhi_h = new double[nbrUcol];  // column of U on host
   double *Ucollo_h = new double[nbrUcol];  // column of U on host
   double *Ucolhi_d;
   double *Ucollo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolhi_d,szUcol);
   cudaMalloc((void**)&Ucollo_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++)
            {
               Ucolhi_h[ix]   = Uhi[rowoff+i][coloff+j];
               Ucollo_h[ix++] = Ulo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolhi_d,Ucolhi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollo_d,Ucollo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      dbl2_back_substitute<<<k,szt>>>
         (szt,k,Ucolhi_d,Ucollo_d,rhshi_d,rhslo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_dbl_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      dbl2_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDhi_d,invDlo_d,rhshi_d,rhslo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *mullapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);

      nbrUcol = nbrUcol - szt*szt; // one tile less used in update
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(xhi,rhshi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlo,rhslo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDhi_h,invDhi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlo_h,invDlo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhi[offset+i][offset+j] = invDhi_h[ix];
            Ulo[offset+i][offset+j] = invDlo_h[ix++];
         }
   }
   free(Dhi_h); free(invDhi_h); free(Ucolhi_h);
   free(Dlo_h); free(invDlo_h); free(Ucollo_h);
}

void GPU_cmplx2_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;       // number of doubles on diagonal tiles
   double *Drehi_h = new double[nbr];    // the diagonal tiles on the host
   double *Drelo_h = new double[nbr];    // low doubles of real parts
   double *Dimhi_h = new double[nbr];    // high doubles of imaginary parts
   double *Dimlo_h = new double[nbr];    // low doubles of imaginary parts
   double *Drehi_d;                      // diagonal tiles on the device
   double *Drelo_d;                      // low doubles of real parts
   double *Dimhi_d;                      // high doubles of imaginary parts
   double *Dimlo_d;                      // low doubles of imaginary parts
   double *invDrehi_h = new double[nbr]; // inverse of tiles on host 
   double *invDrelo_h = new double[nbr]; // low doubles of inverse tiles
   double *invDimhi_h = new double[nbr]; // high doubles of imaginary parts
   double *invDimlo_h = new double[nbr]; // low doubles of imaginary parts
   double *invDrehi_d;                   // invDrehi_d ~ invDrehi_h on device
   double *invDrelo_d;                   // invDrelo_d ~ invDrelo_h on device
   double *invDimhi_d;                   // invDimhi_d ~ invDimhi_h on device
   double *invDimlo_d;                   // invDimlo_d ~ invDimlo_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Drehi_h[ix]   = Urehi[offset+i][offset+j];
            Drelo_h[ix]   = Urelo[offset+i][offset+j];
            Dimhi_h[ix]   = Uimhi[offset+i][offset+j];
            Dimlo_h[ix++] = Uimlo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Drehi_d,sznum);
   cudaMalloc((void**)&Drelo_d,sznum);
   cudaMalloc((void**)&Dimhi_d,sznum);
   cudaMalloc((void**)&Dimlo_d,sznum);
   cudaMalloc((void**)&invDrehi_d,sznum);
   cudaMalloc((void**)&invDrelo_d,sznum);
   cudaMalloc((void**)&invDimhi_d,sznum);
   cudaMalloc((void**)&invDimlo_d,sznum);
   cudaMemcpy(Drehi_d,Drehi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelo_d,Drelo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhi_d,Dimhi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlo_d,Dimlo_h,sznum,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *invlapms = 0.0;
   *mullapms = 0.0;
   *sublapms = 0.0;
   *totlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   cmplx2_invert_tiles<<<nbt,szt>>>
      (szt,   Drehi_d,   Drelo_d,   Dimhi_d,   Dimlo_d,
           invDrehi_d,invDrelo_d,invDimhi_d,invDimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhsrehi_d;                    // right hand side on device
   double *rhsrelo_d;
   double *rhsimhi_d;
   double *rhsimlo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhsrehi_d,szrhs);
   cudaMalloc((void**)&rhsrelo_d,szrhs);
   cudaMalloc((void**)&rhsimhi_d,szrhs);
   cudaMalloc((void**)&rhsimlo_d,szrhs);
   cudaMemcpy(rhsrehi_d,brehi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelo_d,brelo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhi_d,bimhi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlo_d,bimlo,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   cmplx2_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDrehi_d,invDrelo_d,invDimhi_d,invDimlo_d,
                  rhsrehi_d, rhsrelo_d, rhsimhi_d, rhsimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;             // #doubles in column of U
   double *Ucolrehi_h = new double[nbrUcol];  // column of U on host
   double *Ucolrelo_h = new double[nbrUcol];
   double *Ucolimhi_h = new double[nbrUcol];
   double *Ucolimlo_h = new double[nbrUcol];
   double *Ucolrehi_d;
   double *Ucolrelo_d;
   double *Ucolimhi_d;
   double *Ucolimlo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolrehi_d,szUcol);
   cudaMalloc((void**)&Ucolrelo_d,szUcol);
   cudaMalloc((void**)&Ucolimhi_d,szUcol);
   cudaMalloc((void**)&Ucolimlo_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++)
            {
               Ucolrehi_h[ix]   = Urehi[rowoff+i][coloff+j];
               Ucolrelo_h[ix]   = Urelo[rowoff+i][coloff+j];
               Ucolimhi_h[ix]   = Uimhi[rowoff+i][coloff+j];
               Ucolimlo_h[ix++] = Uimlo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolrehi_d,Ucolrehi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelo_d,Ucolrelo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhi_d,Ucolimhi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlo_d,Ucolimlo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      cmplx2_back_substitute<<<k,szt>>>
         (szt,k ,Ucolrehi_d,Ucolrelo_d,Ucolimhi_d,Ucolimlo_d,
                  rhsrehi_d, rhsrelo_d, rhsimhi_d, rhsimlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_cmplx_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      cmplx2_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDrehi_d,invDrelo_d,invDimhi_d,invDimlo_d,
                   rhsrehi_d, rhsrelo_d, rhsimhi_d, rhsimlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *mullapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);

      nbrUcol = nbrUcol - szt*szt; // one tile less used in update
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(xrehi,rhsrehi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelo,rhsrelo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhi,rhsimhi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlo,rhsimlo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDrehi_h,invDrehi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelo_h,invDrelo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhi_h,invDimhi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlo_h,invDimlo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Urehi[offset+i][offset+j] = invDrehi_h[ix];
            Urelo[offset+i][offset+j] = invDrelo_h[ix];
            Uimhi[offset+i][offset+j] = invDimhi_h[ix];
            Uimlo[offset+i][offset+j] = invDimlo_h[ix++];
         }
   }
   free(Drehi_h); free(invDrehi_h); free(Ucolrehi_h);
   free(Drelo_h); free(invDrelo_h); free(Ucolrelo_h);
   free(Dimhi_h); free(invDimhi_h); free(Ucolimhi_h);
   free(Dimlo_h); free(invDimlo_h); free(Ucolimlo_h);
}
