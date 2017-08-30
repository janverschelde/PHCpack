// definition of auxiliary functions, with prototype in get_grid.h

#ifndef __PATH_GPU_MGS_GET_GRID_CU__
#define __PATH_GPU_MGS_GET_GRID_CU__

dim3 get_grid ( int NB, int n_path )
{
   int NBy = 1;
   int NBx = NB;
   while(NBx > 65535)
   {
      NBy++;
      NBx = (NB-1)/NBy + 1;
   }
   return dim3(NBx,NBy,n_path);
}

dim3 get_grid ( int n_job, int BS, int n_path, int n_thread_per_job )
{
   int NB = get_NB(n_job, BS, n_thread_per_job);
   return get_grid(NB, n_path);
}

#endif
