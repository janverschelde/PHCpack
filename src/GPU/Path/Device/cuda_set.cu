// sets the device with error check

#ifndef __CUDA_SET_CU__
#define __CUDA_SET_CU__

#include "cuda_set.h"

void cuda_set()
{
   cudaSetDevice(0);
   if(cudaSuccess
      != cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte))
   {
      std::cout << "Error setting CUDA device!\n" << std::endl;
   }
/*        report only the error when setting CUDA device
   else
   {
      std::cout << "Successfully set CUDA device!\n" << std::endl;
   }
*/
   cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
}

#endif /* __CUDA_SET_CU__ */
