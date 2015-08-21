#ifndef CUDA_SET_CU_
#define CUDA_SET_CU_

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

#endif /* CUDA_SET_CU_ */
