/* Defines the test on the tensor core matrix multiplication
 * with a simple kernel. */

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "smDMMA_dims.h"
#include "smDMMA_host.h"
#include "smDMMA_kernels.h"

using namespace std;

int main ( int argc, char **argv )
{
   int compare_host = 0;
   int number_range = 40;

   if(argc > 1) 
   {
      compare_host = 1;
      number_range = atoi(argv[1]);
      cout << "number range : " << number_range << endl;
   }
   cout << "Initializing ..." << endl;

   int dev = 0; // findCudaDevice(argc, (const char **)argv);

   cudaDeviceProp deviceProp;
   // replaced checkCudaErrors by assert()
   assert((cudaGetDeviceProperties(&deviceProp, dev)) == 0);

   // Double precision Tensor cores require 
   // a GPU of Ampere (SM8X) architecture or higher.
   if(deviceProp.major < 8)
   {
      printf("test_smDMMA requires SM 8.0 or higher.  Exiting ...\n");
      return -1;
   }
   else
      printf("... running on a SM 8.0 or higher ...\n");

   // memory allocations on host and device
   printf("M: %d (%d x %d)\n", M_GLOBAL, M, M_TILES);
   printf("N: %d (%d x %d)\n", N_GLOBAL, N, N_TILES);
   printf("K: %d (%d x %d)\n", K_GLOBAL, K, K_TILES);

   double *A_h = NULL;
   double *B_h = NULL;
   double *C_h = NULL;

   double *result_hD = NULL;
   double *result_host = NULL;

   A_h = (double*) malloc(sizeof(double) * M_GLOBAL * K_GLOBAL);
   B_h = (double*) malloc(sizeof(double) * K_GLOBAL * N_GLOBAL);
   C_h = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);

   if(compare_host > 0)
   {
      result_hD   = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);
      result_host = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);
   }
   double *A = NULL;
   double *B = NULL;
   double *C = NULL;
   double *D = NULL;
   // replaced checkCudaErrors by assert()
   assert((cudaMalloc((void**)&A,sizeof(double)*M_GLOBAL*K_GLOBAL)) == 0);
   assert((cudaMalloc((void**)&B,sizeof(double)*N_GLOBAL*K_GLOBAL)) == 0);
   assert((cudaMalloc((void**)&C,sizeof(double)*M_GLOBAL*N_GLOBAL)) == 0);
   assert((cudaMalloc((void**)&D,sizeof(double)*M_GLOBAL*N_GLOBAL)) == 0);

   assert(((unsigned long long)A) % 128 == 0);
   assert(((unsigned long long)B) % 128 == 0);
   assert(((unsigned long long)C) % 128 == 0);
   assert(((unsigned long long)D) % 128 == 0);

   cout << "Initializing matrices on host ..." << endl;
   init_host_matrices(A_h, B_h, C_h, number_range);

   cout << "Preparing data for GPU ..." << endl;

   assert((cudaMemcpy(A, A_h, sizeof(double) * M_GLOBAL * K_GLOBAL,
                      cudaMemcpyHostToDevice)) == 0);
   assert((cudaMemcpy(B, B_h, sizeof(double) * N_GLOBAL * K_GLOBAL,
                      cudaMemcpyHostToDevice)) == 0);
   assert((cudaMemcpy(C, C_h, sizeof(double) * M_GLOBAL * N_GLOBAL,
                      cudaMemcpyHostToDevice)) == 0);
   assert((cudaMemset(D, 0, sizeof(double) * M_GLOBAL * N_GLOBAL)) == 0);

   const double alpha = 1.0f; // 1.1f;
   const double beta = 1.0f; // 1.2f;

   cudaEvent_t start, stop;

   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   cudaEventRecord(start);

   GPU_DMMA(A, B, C, D, alpha, beta);

   if(compare_host > 0)
   {
      assert(cudaMemcpy(result_hD,D,sizeof(double)*M_GLOBAL*N_GLOBAL,
                        cudaMemcpyDeviceToHost) == 0);
   }
   assert((cudaEventRecord(stop)) == 0);
   assert((cudaEventSynchronize(stop)) == 0);

   if(compare_host > 0)
   {
      cout << "Verifying correctness of the computations ..." << endl;

      memcpy(result_host, C_h, sizeof(double) * M_GLOBAL * N_GLOBAL);

      cout << "Computing matrix matrix multiplication on host ..." << endl;
      matMultiplyOnHost(A_h,B_h,result_host,alpha,beta,M_GLOBAL,K_GLOBAL,
                        K_GLOBAL,N_GLOBAL,M_GLOBAL,N_GLOBAL);

      cout << scientific << setprecision(16);
      cout << "  host : " << result_host[N_GLOBAL-1] << endl;
      cout << "device : " << result_hD[N_GLOBAL-1] << endl;

      size_t number_of_matches = 0;
      for(int i=0; i<N_GLOBAL*M_GLOBAL; i++)
      {
         if(fabs(result_hD[i] - result_host[i]) > 0.1f)
         {
            cout << "mismatch i= " << result_hD[i]
		 << " result_host= " << result_host[i] << endl;
            break;
         }
         else
         {
            number_of_matches++;
         }
      }
      cout << "number_of_matches = " << number_of_matches
	   << " out of = " << N_GLOBAL*M_GLOBAL << endl;
      free(result_hD);
      free(result_host);
   }
   float milliseconds = 0;

   cudaEventElapsedTime(&milliseconds, start, stop);

   cout << "Time: " << milliseconds << " ms" << endl;

   cout << fixed << setprecision(2);

   cout << "FP64 TFLOPS: "
        << (((double)M_GLOBAL*N_GLOBAL*K_GLOBAL*2)/(milliseconds/1000.))/1e12
	<< endl;

   free(A_h);
   free(B_h);
   free(C_h);
   assert((cudaFree((void*)A)) == 0);
   assert((cudaFree((void*)B)) == 0);
   assert((cudaFree((void*)C)) == 0);
   assert((cudaFree((void*)D)) == 0);

   return 0;
}
