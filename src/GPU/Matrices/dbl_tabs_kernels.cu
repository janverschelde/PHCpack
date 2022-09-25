/* The file dbl_tabs_kernels.cu defines the functions with prototypes in
 * the file dbl_tabs_kernels.h. */

#include <iostream>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "dbl_tabs_kernels.h"
#include "dbl_tabs_flopcounts.h"

using namespace std;

__global__ void dbl_small_invert_upper ( int dim, double *U, double *invU )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucol[tabsd_shmemsize];
   __shared__ double invUrows[tabsd_shmemsize];

   double rhs,xval;

   int colidx = dim*(dim-1);          // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column
   rhs = ((double) int(k == dim-1));  // right hand side for each thread
   int rowidx = (dim - 1)*dim + k;    // the row index in the inverse

   __syncthreads();
   if(Ucol[k] != 0.0)
      invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse
   __syncthreads();
   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucol[k] = U[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xval = invUrows[rowidx];

         __syncthreads();
         rhs = rhs - Ucol[i]*xval;    // update right hand side
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucol[k] = U[colidx+k];

      __syncthreads();
      if(Ucol[i] != 0.0) invUrows[rowidx] = rhs/Ucol[i];
      __syncthreads();
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invU[rowidx+k] = invUrows[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void cmplx_small_invert_upper
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolre[tabsd_shmemsize];
   __shared__ double Ucolim[tabsd_shmemsize];
   __shared__ double invUrowsre[tabsd_shmemsize];
   __shared__ double invUrowsim[tabsd_shmemsize];

   double rhsre,rhsim,xvalre,xvalim,accre,accim,det;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolre[k] = Ure[colidx+k];          // load the last column
   Ucolim[k] = Uim[colidx+k];
   rhsre = ((double) int(k == dim-1)); // right hand side for each thread
   rhsim = 0.0;
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k];  // last row of the inverse
   det = Ucolre[k]*Ucolre[k] + Ucolim[k]*Ucolim[k];
   accre = Ucolre[k]/det;
   accim = -Ucolim[k]/det;
   invUrowsre[rowidx] = rhsre*accre - rhsim*accim;
   invUrowsim[rowidx] = rhsim*accre + rhsre*accim;

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsre = ((double) int(k == i)); // set rhs for i-th unit vector
      rhsim = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolre[k] = Ure[colidx+k];
         Ucolim[k] = Uim[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalre = invUrowsre[rowidx];
         xvalim = invUrowsim[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         accre = Ucolre[i]*xvalre - Ucolim[i]*xvalim;
         accim = Ucolim[i]*xvalre + Ucolre[i]*xvalim;
         rhsre = rhsre - accre;
         rhsim = rhsim - accim;
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolre[k] = Ure[colidx+k];
      Ucolim[k] = Uim[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      det = Ucolre[i]*Ucolre[i] + Ucolim[i]*Ucolim[i];
      accre = Ucolre[i]/det;
      accim = -Ucolim[i]/det;
      invUrowsre[rowidx] = rhsre*accre - rhsim*accim;
      invUrowsim[rowidx] = rhsim*accre + rhsre*accim;
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUre[rowidx+k] = invUrowsre[rowidx+k];
      invUim[rowidx+k] = invUrowsim[rowidx+k];
      rowidx = rowidx + dim;
   }
}

void test_dbl_small_invert_upper ( int dim, double *U, double *invU )
{
   double *Ucol = new double[tabsd_shmemsize];
   double *invUrows = new double[tabsd_shmemsize];

   for(int i=0; i<dim*dim; i++)
      cout << "U[" << i << "] : " << U[i] << endl;

   for(int k=0; k<dim; k++)  // compute k-th column of inverse
   {
      double rhs,xval;

      cout << "******** step " << k << " *********" << endl;

      int colidx = dim*(dim-1);          // start with the last column

      Ucol[k] = U[colidx+k];             // load the last column
      rhs = ((double) int(k == dim-1));  // right hand side for each thread
      int rowidx = (dim - 1)*dim + k;    // the row index in the inverse

      cout << "  k : " << k
           << "  i : " << dim-1
           << "  rhs : " << rhs << endl;
      cout << "  k : " << k
           << "  i : " << dim-1
           << "  assigning to rowidx : " << rowidx << endl;
      invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse
      cout << "  value : " << invUrows[rowidx] << endl;
    
      cout << "            entering loop" << endl;

      for(int i=dim-2; i>=0; i--)        // compute row with index i
      {
         rhs = ((double) int(k == i));   // set rhs for i-th unit vector
 
         cout << "  k : " << k
              << "  i : " << i
              << "  rhs : " << rhs << endl;

         for(int j=i+1; j<dim; j++)
         {
            colidx = dim*j;              // need column j of U
            cout << "i : " << i << "  j : " << j << "  k : " << k
                 << "  colidx : " << colidx << endl;
            // Ucol[k] = U[colidx+k];
            for(int L=0; L<dim; L++) Ucol[L] = U[colidx+L];
            rowidx = j*dim + k;          // need solution value
            xval = invUrows[rowidx];
            cout << "Ucol[" << i << "] : " << Ucol[i] << endl;
            cout << "xval : " << xval << endl;
            rhs = rhs - Ucol[i]*xval;    // update right hand side
         }
         rowidx = i*dim + k;             // save in i-th row of inverse

         colidx = dim*i;                 // need column i of U
         // Ucol[k] = U[colidx+k];
         for(int L=0; L<dim; L++) Ucol[L] = U[colidx+L];
         cout << "  k : " << k
              << "  i : " << i
              << "  assigning to rowidx : " << rowidx << endl;
         invUrows[rowidx] = rhs/Ucol[i];
         cout << "  value : " << invUrows[rowidx] << endl;
      }
      rowidx = 0;
      for(int i=0; i<dim; i++)
      {
         invU[rowidx+k] = invUrows[rowidx+k];
         rowidx = rowidx + dim;
      }
   }
   for(int i=0; i<dim*dim; i++)
      cout << "invUrows[" << i << "] : " << invUrows[i] << endl;
}

__global__ void dbl_medium_invert_upper ( int dim, double *U, double *invU )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucol[tabsd_shmemsize];      // one column of U
   __shared__ double invUrow[tabsd_shmemsize];   // one row of invU

   double rhs,xval;

   int colidx = dim*(dim-1);          // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column
   rhs = ((double) int(k == dim-1));  // right hand side for each thread
   int rowidx = (dim - 1)*dim + k;    // the row index in the inverse

   if(Ucol[k] != 0.0)
   {
      invUrow[k] = rhs/Ucol[k];       // last row of the inverse
      invU[rowidx] = invUrow[k];      // store the last row into invU
   }
   __syncthreads();

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucol[k] = U[colidx+k];

         rowidx = j*dim + k;          // need solution value
         invUrow[k] = invU[rowidx];   // load invU row into invUrow
         xval = invUrow[k];

         __syncthreads();
         rhs = rhs - Ucol[i]*xval;    // update right hand side
      }
      colidx = dim*i;                 // need column i of U
      Ucol[k] = U[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      if(Ucol[i] != 0.0)
      {
         invUrow[k] = rhs/Ucol[i];
         invU[rowidx] = invUrow[k];
      }
   }
}

__global__ void cmplx_medium_invert_upper
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolre[tabsd_shmemsize];    // one column of U
   __shared__ double Ucolim[tabsd_shmemsize];    // imaginary parts of U
   __shared__ double invUrowre[tabsd_shmemsize]; // one row of invU
   __shared__ double invUrowim[tabsd_shmemsize]; // imaginary parts of invU

   double rhsre,rhsim,xvalre,xvalim,det,accre,accim;

   int colidx = dim*(dim-1);            // start with the last column

   Ucolre[k] = Ure[colidx+k];           // load the last column
   Ucolim[k] = Uim[colidx+k]; 
   rhsre = ((double) int(k == dim-1));  // right hand side for each thread
   rhsim = 0.0;
   int rowidx = (dim - 1)*dim + k;      // the row index in the inverse

   // invUrow[k] = rhs/Ucol[k];         // last row of the inverse
   det = Ucolre[k]*Ucolre[k] + Ucolim[k]*Ucolim[k];
   accre = Ucolre[k]/det;
   accim = -Ucolim[k]/det;
   invUrowre[k] = rhsre*accre - rhsim*accim;
   invUrowim[k] = rhsim*accre + rhsre*accim;
   invUre[rowidx] = invUrowre[k];       // store the last row into invU
   invUim[rowidx] = invUrowim[k];

   for(int i=dim-2; i>=0; i--)          // compute row with index i
   {
      rhsre = ((double) int(k == i));   // set rhs for i-th unit vector
      rhsim = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;                // need column j of U
         Ucolre[k] = Ure[colidx+k];
         Ucolim[k] = Uim[colidx+k];

         rowidx = j*dim + k;            // need solution value
         invUrowre[k] = invUre[rowidx]; // load invU row into invUrow
         invUrowim[k] = invUim[rowidx];
         xvalre = invUrowre[k];
         xvalim = invUrowim[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;   // update right hand side
         accre = Ucolre[i]*xvalre - Ucolim[i]*xvalim;
         accim = Ucolim[i]*xvalre + Ucolre[i]*xvalim;
         rhsre = rhsre - accre;
         rhsim = rhsim - accim;
      }
      colidx = dim*i;                   // need column i of U
      Ucolre[k] = Ure[colidx+k];
      Ucolim[k] = Uim[colidx+k];
      rowidx = i*dim + k;               // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      det = Ucolre[i]*Ucolre[i] + Ucolim[i]*Ucolim[i];
      accre = Ucolre[i]/det;
      accim = -Ucolim[i]/det;
      invUrowre[k] = rhsre*accre - rhsim*accim;
      invUrowim[k] = rhsim*accre + rhsre*accim;
      invUre[rowidx] = invUrowre[k];
      invUim[rowidx] = invUrowim[k];
   }
}

__global__ void  dbl_invert_tiles ( int dim, double *U, double *invU )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucol[tabsd_shmemsize];      // one column of U
   __shared__ double invUrow[tabsd_shmemsize];   // one row of invU

   double rhs,xval;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column
   rhs = ((double) int(k == dim-1));  // right hand side for each thread
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   invU[rowidx] = 0.0;                // initialize in case Ucol[k] == 0.0
   if(Ucol[k] != 0.0)
   {
      invUrow[k] = rhs/Ucol[k];          // last row of the inverse
      invU[rowidx] = invUrow[k];         // store the last row into invU
   }
   __syncthreads();

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;     // need column j of U
         Ucol[k] = U[colidx+k];

         rowidx = offset + j*dim + k; // need solution value
         invUrow[k] = invU[rowidx];   // load invU row into invUrow
         xval = invUrow[k];

         __syncthreads();
         rhs = rhs - Ucol[i]*xval;    // update right hand side
      }
      colidx = offset + dim*i;        // need column i of U
      Ucol[k] = U[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      invU[rowidx] = 0.0;             // initialize in case Ucol[i] == 0.0
      __syncthreads();
      if(Ucol[i] != 0.0)
      {
         invUrow[k] = rhs/Ucol[i];
         invU[rowidx] = invUrow[k];
      }
   }
}

__global__ void  cmplx_invert_tiles
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolre[tabsd_shmemsize];    // one column of U
   __shared__ double Ucolim[tabsd_shmemsize];    // imaginary parts
   __shared__ double invUrowre[tabsd_shmemsize]; // one row of invU
   __shared__ double invUrowim[tabsd_shmemsize]; // imaginary parts

   double rhsre,rhsim,xvalre,xvalim,det,accre,accim;

   int colidx = offset + dim*(dim-1);   // start with the last column

   Ucolre[k] = Ure[colidx+k];           // load the last column
   Ucolim[k] = Uim[colidx+k];
   rhsre = ((double) int(k == dim-1));  // right hand side for each thread
   rhsim = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];         // last row of the inverse
   det = Ucolre[k]*Ucolre[k] + Ucolim[k]*Ucolim[k];
   invUre[rowidx] = 0.0;
   invUim[rowidx] = 0.0;             // initialize in case det == 0.0
   if(1.0 + det != 1.0)
   {
      accre = Ucolre[k]/det;
      accim = -Ucolim[k]/det;
      invUrowre[k] = rhsre*accre - rhsim*accim;
      invUrowim[k] = rhsim*accre + rhsre*accim;
      invUre[rowidx] = invUrowre[k];       // store the last row into invU
      invUim[rowidx] = invUrowim[k]; 
   }
   __syncthreads();
   for(int i=dim-2; i>=0; i--)          // compute row with index i
   {
      rhsre = ((double) int(k == i));   // set rhs for i-th unit vector
      rhsim = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;       // need column j of U
         Ucolre[k] = Ure[colidx+k];
         Ucolim[k] = Uim[colidx+k];

         rowidx = offset + j*dim + k;   // need solution value
         invUrowre[k] = invUre[rowidx]; // load invU row into invUrow
         invUrowim[k] = invUim[rowidx];
         xvalre = invUrowre[k];
         xvalim = invUrowim[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;   // update right hand side
         accre = Ucolre[i]*xvalre - Ucolim[i]*xvalim;
         accim = Ucolim[i]*xvalre + Ucolre[i]*xvalim;
         rhsre = rhsre - accre;
         rhsim = rhsim - accim;
      }
      colidx = offset + dim*i;        // need column i of U
      Ucolre[k] = Ure[colidx+k];
      Ucolim[k] = Uim[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      det = Ucolre[i]*Ucolre[i] + Ucolim[i]*Ucolim[i];
      invUre[rowidx] = 0.0;
      invUim[rowidx] = 0.0;          // initialize in case det == 0.0
      if(1.0 + det != 1.0)
      {
         accre = Ucolre[i]/det;
         accim = -Ucolim[i]/det;
         invUrowre[k] = rhsre*accre - rhsim*accim;
         invUrowim[k] = rhsim*accre + rhsre*accim;
         invUre[rowidx] = invUrowre[k];
         invUim[rowidx] = invUrowim[k];
      }
   }
}

__global__ void dbl_multiply_inverse
 ( int dim, int idx, double *invU, double *w )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double work[tabsd_shmemsize];      // copy of w

   work[k] = w[rhsoff+k];

   double result = 0.0; // each thread stores its product in result
   double coeff;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeff = invU[offset+k*dim+j]; // thread k does row k
      result = result + coeff*work[j];
   }
   w[rhsoff+k] = result;
}

__global__ void cmplx_multiply_inverse
 ( int dim, int idx, double *invUre, double *invUim,
   double *wre, double *wim )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workre[tabsd_shmemsize];      // copy of wre
   __shared__ double workim[tabsd_shmemsize];      // copy of wim

   workre[k] = wre[rhsoff+k];
   workim[k] = wim[rhsoff+k];

   double resultre = 0.0; // each thread stores its product in result
   double resultim = 0.0;
   double coeffre,coeffim,accre,accim;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffre = invUre[offset+k*dim+j]; // thread k does row k
      coeffim = invUim[offset+k*dim+j];
      // result = result + coeff*work[j];
      accre = coeffre*workre[j] - coeffim*workim[j];
      accim = coeffim*workre[j] + coeffre*workim[j];
      resultre = resultre + accre;
      resultim = resultim + accim;
   }
   wre[rhsoff+k] = resultre;
   wim[rhsoff+k] = resultim;
}

__global__ void dbl_back_substitute
 ( int dim, int idx, double *U, double *w )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrk[tabsd_shmemsize];    // copy of w
   __shared__ double sol[tabsd_shmemsize];    // solution to update with

   wrk[k] = w[B*dim+k];    // block B updates B-th slice of w
   sol[k] = w[idx*dim+k];  // solution that is back substituted

   double result = 0.0; // each thread stores its product in result
   double coeff;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeff = U[offset+k*dim+j];
      result = result + coeff*sol[j];
   }
   wrk[k] = wrk[k] - result; // subtract product
   w[B*dim+k] = wrk[k];
}

__global__ void cmplx_back_substitute
 ( int dim, int idx, double *Ure, double *Uim, double *wre, double *wim )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkre[tabsd_shmemsize];    // copy of wre
   __shared__ double wrkim[tabsd_shmemsize];    // copy of wim
   __shared__ double solre[tabsd_shmemsize];    // solution to update with
   __shared__ double solim[tabsd_shmemsize];    // imaginary parts

   wrkre[k] = wre[B*dim+k];    // block B updates B-th slice of w
   wrkim[k] = wim[B*dim+k];
   solre[k] = wre[idx*dim+k];  // solution that is back substituted
   solim[k] = wim[idx*dim+k];

   double resultre = 0.0; // each thread stores its product in result
   double resultim = 0.0;
   double coeffre,coeffim,accre,accim;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffre = Ure[offset+k*dim+j];
      coeffim = Uim[offset+k*dim+j];
      // result = result + coeff*sol[j];
      accre = coeffre*solre[j] - coeffim*solim[j];
      accim = coeffim*solre[j] + coeffre*solim[j];
      resultre = resultre + accre;
      resultim = resultim + accim;
   }
   wrkre[k] = wrkre[k] - resultre; // subtract product
   wrkim[k] = wrkim[k] - resultim;
   wre[B*dim+k] = wrkre[k];
   wim[B*dim+k] = wrkim[k];
}

void GPU_dbl_upper_inverse
 ( int dim, double **U, double **invU, double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *U_h = new double[szU];     // U_h stores the columns of U 
   double *U_d;                       // U_d is U_h on the device
   double *invU_h = new double[szU];  // the columns of the inverse
   double *invU_d;                    // invU_d is invU_h on the device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++) U_h[ix++] = U[i][j];

   // test_small_invert_upper(dim,U_h,invU_h); // only for debugging

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&U_d,szmat);
   cudaMalloc((void**)&invU_d,szmat);
   cudaMemcpy(U_d,U_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      dbl_small_invert_upper<<<1,dim>>>(dim,U_d,invU_d);
   else
      dbl_medium_invert_upper<<<1,dim>>>(dim,U_d,invU_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invU_h,invU_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) invU[i][j] = invU_h[ix++];

   free(U_h); free(invU_h);
}

void GPU_cmplx_upper_inverse
 ( int dim, double **Ure, double **Uim, double **invUre, double **invUim,
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Ure_h = new double[szU];    // Ure_h stores real U 
   double *Uim_h = new double[szU];    // Uim_h stores imaginary U 
   double *Ure_d;                      // Ure_d is Ure_h on the device
   double *Uim_d;                      // Uim_d is Uim_h on the device
   double *invUre_h = new double[szU]; // real parts of the inverse
   double *invUim_h = new double[szU]; // imaginary parts of the inverse
   double *invUre_d;                   // invUre_d is invUre_h on the device
   double *invUim_d;                   // invUim_d is invUim_h on the device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Ure_h[ix]   = Ure[i][j];
         Uim_h[ix++] = Uim[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Ure_d,szmat);
   cudaMalloc((void**)&Uim_d,szmat);
   cudaMalloc((void**)&invUre_d,szmat);
   cudaMalloc((void**)&invUim_d,szmat);
   cudaMemcpy(Ure_d,Ure_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uim_d,Uim_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      cmplx_small_invert_upper<<<1,dim>>>(dim,Ure_d,Uim_d,invUre_d,invUim_d);
   else
      cmplx_medium_invert_upper<<<1,dim>>>(dim,Ure_d,Uim_d,invUre_d,invUim_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUre_h,invUre_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUim_h,invUim_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUre[i][j] = invUre_h[ix];
         invUim[i][j] = invUim_h[ix++];
      }

   free(Ure_h); free(invUre_h);
   free(Uim_h); free(invUim_h);
}

void GPU_dbl_upper_tiled_solver
 ( int dim, int szt, int nbt, double **U, double *b, double *x,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *D_h = new double[nbr];    // the diagonal tiles on the host
   double *D_d;                      // diagonal tiles on the device
   double *invD_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invD_d;                   // invD_d is invD_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++) D_h[ix++] = U[offset+i][offset+j];
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&D_d,sznum);
   cudaMalloc((void**)&invD_d,sznum);
   cudaMemcpy(D_d,D_h,sznum,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *invlapms = 0.0;
   *mullapms = 0.0;
   *sublapms = 0.0;
   *totlapms = 0.0;
   *addcnt = 0; *mulcnt = 0; *divcnt = 0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   dbl_invert_tiles<<<nbt,szt>>>(szt,D_d,invD_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhs_d;                    // right hand side on device
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhs_d,szrhs);
   cudaMemcpy(rhs_d,b,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl_multiply_inverse<<<1,szt>>>(szt,nbt-1,invD_d,rhs_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;         // #doubles in column of U
   double *Ucol_h = new double[nbrUcol];  // column of U on host
   double *Ucol_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucol_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++) Ucol_h[ix++] = U[rowoff+i][coloff+j];
      }
      cudaMemcpy(Ucol_d,Ucol_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      dbl_back_substitute<<<k,szt>>>(szt,k,Ucol_d,rhs_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_dbl_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      dbl_multiply_inverse<<<1,szt>>>(szt,k-1,invD_d,rhs_d);
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

   cudaMemcpy(x,rhs_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invD_h,invD_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++) U[offset+i][offset+j] = invD_h[ix++];
   }
   free(D_h); free(invD_h); free(Ucol_h);
}

void GPU_cmplx_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Ure, double **Uim,
   double *bre, double *bim, double *xre, double *xim,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *Dre_h = new double[nbr];    // the diagonal tiles on the host
   double *Dim_h = new double[nbr];    // imaginary parts of diagonal tiles
   double *Dre_d;                      // diagonal tiles on the device
   double *Dim_d;                      // imaginary parts of diagonal tiles
   double *invDre_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invDim_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invDre_d;                   // invDre_d is invDre_h on device
   double *invDim_d;                   // invDim_d is invDim_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Dre_h[ix]   = Ure[offset+i][offset+j];
            Dim_h[ix++] = Uim[offset+i][offset+j];
         }
   }
/*
   cout << "after defining the diagonal tiles ..." << endl;
   for(int i=0; i<nbr; i++)
      cout << "D[" << i << "] : " << Dre_h[i] << "  " << Dim_h[i] << endl;
 */
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Dre_d,sznum);
   cudaMalloc((void**)&Dim_d,sznum);
   cudaMalloc((void**)&invDre_d,sznum);
   cudaMalloc((void**)&invDim_d,sznum);
   cudaMemcpy(Dre_d,Dre_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dim_d,Dim_h,sznum,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *invlapms = 0.0;
   *mullapms = 0.0;
   *sublapms = 0.0;
   *totlapms = 0.0;
   *addcnt = 0; *divcnt = 0; *mulcnt = 0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   cmplx_invert_tiles<<<nbt,szt>>>(szt,Dre_d,Dim_d,invDre_d,invDim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   // testing, add another cudaMemcpy ...
/*
   cudaMemcpy(invDre_h,invDre_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDim_h,invDim_d,sznum,cudaMemcpyDeviceToHost);
   cout << "after inverting the diagonal tiles ..." << endl;
   for(int i=0; i<nbr; i++)
      cout << "invD[" << i << "] : "
           << invDre_h[i] << "  " << invDim_h[i] << endl;

   cout << "before copying b to rhs ... " << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : " << bre[i] << "  " << bim[i] << endl;
 */
   double *rhsre_d;                    // right hand side on device
   double *rhsim_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhsre_d,szrhs);
   cudaMalloc((void**)&rhsim_d,szrhs);
   cudaMemcpy(rhsre_d,bre,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsim_d,bim,szrhs,cudaMemcpyHostToDevice);
/*
   cout << "cmplx_multiply_inverse szt = " << szt
        << ", nbt-1 = " << nbt-1 << endl;
 */
   cudaEventRecord(start);
   cmplx_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDre_d,invDim_d,rhsre_d,rhsim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;           // #doubles in column of U
   // int nbrUcol = dim*szt*szt; looks weird if nbt = 1 (one tile) ...
   double *Ucolre_h = new double[nbrUcol];  // column of U on host
   double *Ucolim_h = new double[nbrUcol];  // imaginary parts of the column
   double *Ucolre_d;
   double *Ucolim_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolre_d,szUcol);
   cudaMalloc((void**)&Ucolim_d,szUcol);

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
               Ucolre_h[ix]   = Ure[rowoff+i][coloff+j];
               Ucolim_h[ix++] = Uim[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolre_d,Ucolre_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolim_d,Ucolim_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      cmplx_back_substitute<<<k,szt>>>
         (szt,k,Ucolre_d,Ucolim_d,rhsre_d,rhsim_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_cmplx_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      cmplx_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDre_d,invDim_d,rhsre_d,rhsim_d);
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

   cudaMemcpy(xre,rhsre_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xim,rhsim_d,szrhs,cudaMemcpyDeviceToHost);
/*
   cout << "after copying rhs to x ..." << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : " << xre[i] << "  " << xim[i] << endl;
 */
   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDre_h,invDre_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDim_h,invDim_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Ure[offset+i][offset+j] = invDre_h[ix];
            Uim[offset+i][offset+j] = invDim_h[ix++];
         }
   }
   free(Dre_h); free(invDre_h); free(Ucolre_h);
   free(Dim_h); free(invDim_h); free(Ucolim_h);
}
