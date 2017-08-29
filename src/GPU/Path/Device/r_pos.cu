__device__ inline int r_pos ( int x, int y, int cols )
{
   return cols*(cols+1)/2 -1 - (y*(y+1)/2 -(x-y));
}
