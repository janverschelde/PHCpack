__device__ inline int r_pos ( int x, int y, int cols )
{
   return cols*(cols+1)/2 -1 - (y*(y+1)/2 -(x-y));
}

#include "mgs_small.cu"
#include "mgs_large.cu"

int GPU_MGS
 ( const CPUInstHom& hom, CT*& sol_gpu, CT*& matrix_gpu_q,
   CT*& matrix_gpu_r, CT* V)
{
   cout << "GPU MGS" << endl;

   // CUDA configuration
   cuda_set();

   GPUWorkspace workspace(0, 0, 0, hom.n_eq, hom.dim, 0, 0, 0, 1);

   workspace.init_V_value(V);

   std::cout << "workspace.n_path = " << workspace.n_path << std::endl;

   struct timeval start, end;
   long seconds, useconds;
   gettimeofday(&start, NULL);

   if(hom.dim <= BS_QR)
   {
      mgs_small_dynamic
         (workspace.V, workspace.R, workspace.sol, hom.n_eq, hom.dim+1,
          workspace.small_mgs_size, workspace.n_matrix, workspace.n_matrix_R,
          workspace.n_path);
      // mgs_small_template
      //    (workspace.V, workspace.R, workspace.sol, hom.n_eq, hom.dim+1,
      //     workspace.n_matrix, workspace.n_matrix_R, workspace.n_path);
   }
   else
   {
      mgs_large_block(workspace.V, workspace.R, workspace.P, workspace.sol,
                      hom.n_eq, hom.dim+1);
      // mgs_large_orig(workspace.V, workspace.R, workspace.sol,
      //                hom.n_eq, hom.dim+1);
      // mgs_large_old(workspace.V, workspace.R, workspace.sol,
      //               hom.n_eq, hom.dim+1);
   }

   sol_gpu = workspace.get_sol();

   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   double timeMS_MGS_GPU = seconds*1000 + useconds/1000.0;
   double timeSec_MGS_GPU = timeMS_MGS_GPU/1000;
   std::cout << "GPU Time MS  = " << timeMS_MGS_GPU << std::endl;
   std::cout << "GPU Time Sec = " << timeSec_MGS_GPU << std::endl;

   matrix_gpu_q = workspace.get_matrix();
   matrix_gpu_r = workspace.get_matrix_r();

   cudaDeviceReset();
   return 0;
}

int GPU_MGS_Mult
 ( const CPUInstHom& hom, CT**& sol_gpu, CT**& matrix_gpu_q,
   CT**& matrix_gpu_r, CT* V, int n_path ) 
{
   cout << "GPU MGS Mult" << endl;

   // CUDA configuration
   cuda_set();

   GPUWorkspace workspace(0, 0, 0, hom.n_eq, hom.dim, 0, 0, 0, n_path);

   std::cout << "n_path = " << n_path << std::endl;

   workspace.init_V_value(V);

   std::cout << "workspace.n_path = " << workspace.n_path << std::endl;

   struct timeval start, end;
   long seconds, useconds;
   gettimeofday(&start, NULL);

   if(hom.dim <= 32)
   {
      // Up to 32 for multiple
      mgs_small_dynamic
         (workspace.V, workspace.R, workspace.sol, hom.n_eq, hom.dim+1,
          workspace.small_mgs_size, workspace.n_matrix, workspace.n_matrix_R,
          workspace.n_path);
      // mgs_small_template
      //   (workspace.V,workspace.R,workspace.sol,hom.n_eq,hom.dim+1,
      //    workspace.n_matrix, workspace.n_matrix_R, workspace.n_path);
      // mgs_small(workspace.V,workspace.R,workspace.sol,hom.n_eq,hom.dim+1,
      //           workspace.n_matrix,workspace.n_path);
      // mgs_small1(workspace.V, workspace.R,workspace.sol,hom.n_eq,hom.dim+1,
      //            workspace.n_matrix,workspace.n_matrix_R,workspace.n_path);
   }
   else
   {
      // Only work for single
      mgs_large_block(workspace.V, workspace.R, workspace.P,
                      workspace.sol, hom.n_eq, hom.dim+1);
   }
   sol_gpu = workspace.get_sol_mult();

   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   double timeMS_MGS_GPU = seconds*1000 + useconds/1000.0;
   double timeSec_MGS_GPU = timeMS_MGS_GPU/1000;
   std::cout << "GPU Time MS  = " << timeMS_MGS_GPU << std::endl;
   std::cout << "GPU Time Sec = " << timeSec_MGS_GPU << std::endl;

   cudaDeviceReset();
   return 0;
}
