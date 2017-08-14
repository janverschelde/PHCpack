#include "mgs_large_backsub.cu"
#include "mgs_large_norm.cu"
#include "mgs_large_row_reduce.cu"

void mgs_large_reduce ( GT* V, GT* R, GT* sol, int rows, int cols );

void mgs_large ( GT* V, GT* R, GT* sol, int rows, int cols ) 
{
   mgs_large_reduce(V, R, sol, rows, cols);
   mgs_large_backsubstitution(R, sol, rows, cols);
}

void mgs_large_reduce ( GT* V, GT* R, GT* sol, int rows, int cols )
{
   int BS = min(BS_QR,rows);
   // int BS = 32;

   int rowsLog2 = log2ceil(rows); // ceil for sum reduction
   // int dimR = cols*(cols+1)/2;

   T* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(T));

   T* sums_global; // norm of the pivot column
   cudaMalloc((void**)&sums_global,maxrounds*sizeof(T));

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
      mgs_large_normalize_kernel1<<<rf,BS>>>
         (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,
          pivnrm,lastBSLog2,sums_global);
      mgs_large_normalize_kernel2<<<rf,BS>>>
         (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,
          pivnrm,lastBSLog2,sums_global);
      // XXX BS should be greater than maxround
      mgs_large_row_reduce_kernel<<<cols-piv-1,BS>>>
         (V,R,cols,rows,rowsLog2, piv,rf,rfLog2,BS,BSLog2,pivnrm,lastBSLog2);
   }
}
