// Using a GPU to compute the square root of n complex numbers.

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include "complexH.h"
#include "complex_gd.cu"

using namespace std;

__global__ void squareRoot
 ( int n, gd_complex *x, gd_complex *y )
// Applies Newton's method to compute the square root 
// of the n numbers in x and places the results in y.
{
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   gd_complex inc;
   gd_complex c = x[i];
   gd_complex r = c;
   for(int j=0; j<10; j++)
   {
      inc = r + r;
      inc = (r*r - c)/inc;
      r = r - inc;
   }
   y[i] = r;
}

complexH<double> Square ( const complexH<double> z )
// Returns z*z for testing purposes.
{
   complexH<double> result;
   result.real = z.real*z.real - z.imag*z.imag;
   result.imag = 2*z.real*z.imag;
   return result;
}

int main ( int argc, char*argv[] )
{
   if(argc < 5)
   {
      cout << "call with 4 arguments : " << endl;
      cout << "dimension, block size, frequency, and check (0 or 1)" << endl;
   }
   else
   {
      int n = atoi(argv[1]); // dimension
      int w = atoi(argv[2]); // block size
      int f = atoi(argv[3]); // frequency
      int t = atoi(argv[4]); // test or not
      // we generate n complex numbers on the host
      complexH<double> *xhost = new complexH<double>[n];
      for(int i=0; i<n; i++) 
      {
         xhost[i].real = (double) i+1;
         xhost[i].imag = 0.0;
      }
      // we copy the n complex numbers to the device
      size_t s = n*sizeof(gd_complex);
      gd_complex *xdevice;
      cudaMalloc((void**)&xdevice,s);
      cudaMemcpy(xdevice,xhost,s,cudaMemcpyHostToDevice);
      // allocate memory for the result
      gd_complex *ydevice;
      cudaMalloc((void**)&ydevice,s);
      // invoke the kernel with n/w blocks per grid
      // and w threads per block
      for(int i=0; i<f; i++)
         squareRoot<<<n/w,w>>>(n,xdevice,ydevice);
      // copy results from device to host
      complexH<double> *yhost = new complexH<double>[n];
      cudaMemcpy(yhost,ydevice,s,cudaMemcpyDeviceToHost);
      if(t == 1) // test the result
      {
         for(int k=0; k < n; k++)
         {
            cout << "testing number " <<  k << endl;
            cout << setprecision(16);
            cout << "        x = " << xhost[k].real << endl;
            cout << "  sqrt(x) = " << yhost[k].real << endl;
            complexH<double> z = Square(yhost[k]);
            cout << "sqrt(x)^2 = " << z.real << endl;
         }
      }
   }
   return 0;
}
