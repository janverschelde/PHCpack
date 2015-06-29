#ifndef CUDA_SET_CU_
#define CUDA_SET_CU_

#include "cuda_set.h"

void cuda_set(){
	cudaSetDevice(0);
	if (cudaSuccess
			!= cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte))
		std::cout << "Error!\n" << std::endl;
	else {
		std::cout << "Success!\n" << std::endl;
	}
	cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
}

#endif /* CUDA_SET_CU_ */
