// include <sys/time.h> needed for gettimeofday()
#include <sys/time.h>

#include "mgs_small.cu"
// #include "mgs_large_row.cu"
#include "mgs_large_block.cu"
#include "cuda_set.cu"

template <class ComplexType, class RealType>
int GPU_MGS
 ( int neq, int nvr, int BS_QR=256 ) // default for double precision
{
   cuda_set();

   ComplexType* V;
   ComplexType* R;
   ComplexType* P;
   ComplexType* sol;

   struct timeval start, end;
   long seconds, useconds;
   gettimeofday(&start, NULL);

   if(nvr <= BS_QR)
   {
      int small_mgs_size;

      mgs_small_dynamic<ComplexType, RealType>
         (V, R, sol, neq, nvr+1, small_mgs_size);
   }
   else
   {
      mgs_large_block<ComplexType, RealType>(V, R, P, sol, neq, nvr+1);
   }

   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   double timeMS_MGS_GPU = seconds*1000 + useconds/1000.0;
   double timeSec_MGS_GPU = timeMS_MGS_GPU/1000;
   std::cout << "GPU Time MS  = " << timeMS_MGS_GPU << std::endl;
   std::cout << "GPU Time Sec = " << timeSec_MGS_GPU << std::endl;

   cudaDeviceReset();

   return 0;
}
