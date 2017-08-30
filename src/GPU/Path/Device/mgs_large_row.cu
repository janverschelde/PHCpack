// row reduction kernels for large dimensional problems

#include "mgs_large_backsub.cu"
#include "mgs_large_norm.cu"
#include "mgs_large_row_reduce.cu"
#include "log2ceil.h"

template <class ComplexType, class RealType>
void mgs_large_reduce
 ( ComplexType* V, ComplexType* R, ComplexType* sol, int rows, int cols );
//   int maxrounds=128, int shmemsize=512, int BS_QR=256 );
// default values are for double precision

template <class ComplexType, class RealType>
void mgs_large
 ( ComplexType* V, ComplexType* R, ComplexType* sol, int rows, int cols,
   int maxrounds=128, int shmemsize=512, int BS_QR=256 )
// default values are for double precision
{
   mgs_large_reduce<ComplexType,RealType>
      (V, R, sol, rows, cols, maxrounds, BS_QR);
   mgs_large_backsubstitution<ComplexType,RealType>
      (R, sol, rows, cols, BS_QR);
}

template <class ComplexType, class RealType>
void mgs_large_reduce
 ( ComplexType* V, ComplexType* R, ComplexType* sol, int rows, int cols,
   int maxrounds=128, int shmemsize=512, int BS_QR=256 )
// default values are for double precision
{
   int BS = min(BS_QR,rows);
   // int BS = 32;

   int rowsLog2 = log2ceil(rows); // ceil for sum reduction
   // int dimR = cols*(cols+1)/2;

   RealType* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(RealType));

   RealType* sums_global; // norm of the pivot column
   cudaMalloc((void**)&sums_global,maxrounds*sizeof(RealType));

   int rf = ceil(((double) rows)/BS);
   int rfLog2 = log2ceil(rf);
   int BSLog2 = log2ceil(BS);
   int lastBSLog2 = log2ceil(rows-BS*(rf-1));

   /*std::cout << "BS     = " << BS << std::endl;
     std::cout << "rf     = " << rf << std::endl;
     std::cout << "rfLog2 = " << rfLog2 << std::endl;
     std::cout << "BSLog2 = " << BSLog2 << std::endl;
     std::cout << "lastBSLog2 = " << lastBSLog2 << std::endl;*/

   for(int piv=0; piv<cols-1; piv++)
   {
      mgs_large_normalize_kernel1<ComplexType,RealType><<<rf,BS>>>
         (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,
          pivnrm,lastBSLog2,sums_global,maxrounds,shmemsize);
      mgs_large_normalize_kernel2<ComplexType,RealType><<<rf,BS>>>
         (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,
          pivnrm,lastBSLog2,sums_global,maxrounds,shmemsize);
      // XXX BS should be greater than maxround
      mgs_large_row_reduce_kernel<ComplexType,RealType><<<cols-piv-1,BS>>>
         (V,R,cols,rows,rowsLog2, piv,rf,rfLog2,BS,BSLog2,pivnrm,lastBSLog2,
          maxrounds,BS_QR);
   }
}
