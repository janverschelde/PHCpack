// Using a GPU to compute the square root of n complex numbers.

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <gdd_basic.cu>
#include <gqd_basic.cu>
#include "complexD.cu"
#include "complexH.h"

template <class realH, class realD>
int test ( int dim, int blk, int frq, int prc, int tst );
/*
 * Computes the square root of the first dim numbers, starting at 1.
 * Values for the four parameters are passed via the command line.
 * The type realH is the real type for the host,
 * while realD is the real type for the device.
 *
 * ON ENTRY :
 *   dim    dimension, number of numbers tests;
 *   blk    block size, should be less than or equal to dim;
 *   frq    frequency of the tests, for timing purposes;
 *   prc    working precision, 0 for double, 1 for double double,
 *          and 2 for quad double;
 *   tst    0 for silent or 1 for extra test with output. */

using namespace std;

template <class real>
__global__ void squareRoot
 ( int n, complexD<real> *x, complexD<real> *y )
// Applies Newton's method to compute the square root 
// of the n numbers in x and places the results in y.
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;

   complexD<real> inc;
   complexD<real> c = x[i];
   complexD<real> r = c;

   for(int j=0; j<10; j++)
   {
      inc = r + r;
      inc = (r*r - c)/inc;
      r = r - inc;
   }
   y[i] = r;
}

template <class real>
complexH<real> Square ( const complexH<real> z )
// Returns z*z for testing purposes.
{
   complexH<real> result;
   result.real = z.real*z.real - z.imag*z.imag;
   result.imag = 2*z.real*z.imag;
   return result;
}

template <class realH, class realD>
int test ( int dim, int blk, int frq, int prc, int tst )
{

   // we generate dim complex numbers on the host
   complexH<realH> *xhost = new complexH<realH>[dim];
   for(int i=0; i<dim; i++) 
   {
      xhost[i].real = realH((double) i+1);
      xhost[i].imag = realH(0.0);
   }
   // we copy the dim complex numbers to the device
   size_t s = dim*sizeof(complexD<realD>);
   complexD<realD> *xdevice;
   cudaMalloc((void**)&xdevice,s);
   cudaMemcpy(xdevice,xhost,s,cudaMemcpyHostToDevice);
   // allocate memory for the result
   complexD<realD> *ydevice;
   cudaMalloc((void**)&ydevice,s);
   // invoke the kernel with dim/blk blocks per grid
   // and blk threads per block
   for(int i=0; i<frq; i++)
      squareRoot<<<dim/blk,blk>>>(dim,xdevice,ydevice);
   // copy results from device to host
   complexH<realH> *yhost = new complexH<realH>[dim];
   cudaMemcpy(yhost,ydevice,s,cudaMemcpyDeviceToHost);
   if(tst == 1) // test the result
   {
      for(int k=0; k<dim; k++)
      {
         cout << "testing number " <<  k << endl;
         if(prc == 0) cout << setprecision(16);
         if(prc == 1) cout << setprecision(32);
         if(prc == 2) cout << setprecision(64);
         cout << "        x = " << xhost[k].real << endl;
         cout << "  sqrt(x) = " << yhost[k].real << endl;
         complexH<realH> z = Square<realH>(yhost[k]);
         cout << "sqrt(x)^2 = " << z.real << endl;
      }
   }
   return 0;
}

int main ( int argc, char*argv[] )
{
   if(argc < 6)
   {
      cout << "call with 5 arguments : " << endl;
      cout << "dimension, block size, frequency, ";
      cout << "precision (0, 1, or 2), and check (0 or 1)" << endl;

      return 0;
   }
   else
   {
      int n = atoi(argv[1]); // dimension
      int w = atoi(argv[2]); // block size
      int f = atoi(argv[3]); // frequency
      int p = atoi(argv[4]); // precision
      int t = atoi(argv[5]); // test or not

      if(p == 0) return test<double,double>(n,w,f,p,t);
      if(p == 1) return test<dd_real,gdd_real>(n,w,f,p,t);
      if(p == 2) return test<qd_real,gqd_real>(n,w,f,p,t);
      
      cout << "invalid value for the precision, use 0, 1, or 2" << endl;
   }
}
