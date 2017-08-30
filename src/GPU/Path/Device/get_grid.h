// auxiliary functions needed in mgs_large_block_reduce.cu

#ifndef __PATH_GPU_MGS_GET_GRID_H__
#define __PATH_GPU_MGS_GET_GRID_H__

dim3 get_grid ( int NB, int n_path=1 );

dim3 get_grid ( int n_job, int BS, int n_path, int n_thread_per_job = 1 );

#endif
