/* Defines the function with prototype in write_gpu_timings.h. */

#include <iostream>
#include <iomanip>
#include "write_gpu_timings.h"

void write_GPU_timings
 ( double cnvlapms, double addlapms, double elapsedms, double walltimesec )
{
   using namespace std;

   cout << fixed << setprecision(2);
   cout << "Time spent by convolution kernels : ";
   cout << cnvlapms << " milliseconds." << endl;
   cout << "Time spent by addition kernels    : ";
   cout << addlapms << " milliseconds." << endl;
   cout << "Time spent by all kernels         : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "Total wall clock computation time : ";
   cout << fixed << setprecision(3) << walltimesec << " seconds." << endl;
   cout << scientific << setprecision(16);
}
