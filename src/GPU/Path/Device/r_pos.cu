// auxiliary function for modified Gram-Schmidt method

#ifndef __PATH_GPU_MGS_R_POS_CU__
#define __PATH_GPU_MGS_R_POS_CU__

__device__ inline int r_pos ( int x, int y, int cols )
{
   return cols*(cols+1)/2 -1 - (y*(y+1)/2 -(x-y));
}

#endif
