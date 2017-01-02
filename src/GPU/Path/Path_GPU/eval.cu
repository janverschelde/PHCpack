#if(path_precision == 0)
#include "eval_mon_d.cu"
#elif(path_precision == 1)
#include "eval_mon_dd.cu"
#else
#include "eval_mon_qd.cu"
#endif

#include "eval_sum.cu"
#include "eval_base.cu"
#include "eval_coef.cu"
#include "eval_mult.cu"
#include "path_init.cu"

void eval ( GPUWorkspace& workspace, const GPUInst& inst )
{
   // std::cout << "inst.dim = " << inst.dim << std::endl;

   if(inst.PED_hom == true)
   {
      // Evaluate homotopy coefficient
      eval_coef_kernel<<<inst.coef_grid, inst.coef_BS>>>
         (workspace.coef, inst.coef, inst.n_coef, workspace.t,
          workspace.one_minor_t, workspace.workspace_size);
   }
   else
   {
      // Copy homotopy coefficient, should be combined to base kernel
      eval_coef_copy_kernel<<<inst.coef_grid, inst.coef_BS>>>
         (workspace.coef, inst.coef, inst.n_coef);
   }

   if(inst.max_deg_base!= NULL)
   {
      eval_base_table_kernel<<<1, inst.dim>>>
         (inst.dim, inst.max_deg_base, inst.base_table_start,
          workspace.x, workspace.deg_table);

      std::cout << "n_constant = " << inst.n_constant
                << " n_mon_base_start = " << inst.n_mon_base_start
                << std::endl;

      eval_base_kernel<<<1, inst.n_mon_base>>>
         (inst.n_mon_base, inst.mon_pos, inst.mon_exp,
          inst.mon_pos_start + inst.n_mon_base_start,
          inst.base_table_start, workspace.deg_table,
          workspace.coef + inst.n_mon_base_start);
   }

   // CT** gpu_workspace_all = new CT*[workspace.n_path];
   // CT** gpu_matrix_all = new CT*[workspace.n_path];
   // for(int path_idx=0; path_idx<workspace.n_path; path_idx++)
   // {
   //    gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
   //    gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
   //    std::cout << "path_idx = " << path_idx << std::endl;
   //    for(int coef_idx=0; coef_idx<10; coef_idx++)
   //    {
   //       std::cout << coef_idx << " "
   //                 <<  gpu_workspace_all[path_idx][coef_idx];
   //    }
   // }

   eval_mon(workspace, inst);
   eval_sum(workspace, inst);
}

int GPU_Eval
 ( const CPUInstHom& hom, CT* cpu_sol0, CT* cpu_t, CT**& gpu_workspace_all,
   CT**& gpu_matrix_all, int n_path, int* x_t_idx, int n_predictor )
{

   cout << "GPU Eval" << endl;
   // CUDA configuration
   cuda_set();

   GPUInst inst(hom, n_path);
   GPUWorkspace workspace
     (inst.mon_pos_size, inst.n_coef, inst.n_constant, inst.n_eq,
      inst.dim, n_predictor, inst.alpha, inst.base_table_size, n_path);

   dim3 update_t_grid = get_grid(n_path,inst.predict_BS,1);
   std::cout << "n_path = " << n_path << std::endl;
   std::cout << "inst.n_sum_zero = " << inst.n_sum_zero << std::endl;
   path_mult_init_kernel<<<update_t_grid, inst.predict_BS>>>
     (workspace.t_mult, workspace.t_last_mult, workspace.delta_t_mult,
      workspace.t_array, workspace.one_minor_t, workspace.alpha_gpu,
      workspace.path_success, workspace.newton_success, workspace.end_range,
      workspace.n_success, workspace.max_f_val_last_gpu, 
      workspace.n_point_mult, workspace.x_t_idx_mult,
      workspace.workspace_size, workspace.path_idx, n_path,
      workspace.matrix, inst.n_sum_zero, inst.sum_zeros,
      workspace.n_predictor);

   // x_t_idx, x_array, t_array
   if(x_t_idx!=NULL)
   {
      workspace.update_x_t_idx_all(x_t_idx);
      workspace.update_x_t_value_array(cpu_sol0, cpu_t, x_t_idx);
   }
   // x_mult, t_mult
   if(n_path>1)
   {
      workspace.update_x_mult_vertical(cpu_sol0, x_t_idx);
      workspace.update_t_value_mult2(cpu_t, x_t_idx);
   }
   // workspace.init_workspace_eq(inst.n_pos_total_eq, workspace.n_path);

   bool mult = false;
   if(n_path>1)
   {
      mult = true;
   }

   struct timeval start, end;
   long seconds, useconds;
   double timeMS_gpu;
   gettimeofday(&start, NULL);

   if(mult)
   {
      eval_mult(workspace, inst);
   }
   else
   {
      eval(workspace, inst);
   }

   // Get result
   gpu_matrix_all = workspace.get_matrix_mult();
   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   timeMS_gpu = seconds*1000 + useconds/1000.0;
   std::cout.precision(8);
   cout << "Path GPU Eval    Time: "<< timeMS_gpu << endl;

   // Get workspace for test
   gpu_workspace_all = workspace.get_workspace_mult();

   cudaDeviceReset();
   return 0;
}
