#ifndef PATH_GPU_CU_
#define PATH_GPU_CU_

#include "cuda_set.cu"

#include "complex.cu"
#include "parameter.h"

#include "path_gpu_data.cu"

#include "predict.cu"
#include "newton.cu"
#include "path_gpu_mult.cu"

bool path_single
 ( GPUWorkspace& workspace, GPUInst& inst, Parameter path_parameter,
   CT cpu_t, int n_path, int inverse = 0, int verbose = 0 )
{
   bool debug = false; // debug = true;

   /* bool Record = false; // Record = true;

     if(Record)
     {
        workspace.path_record.init_record(path_parameter.max_step);
     }
    */

   int n_point = 1;
   int n_step = 0;

   // Parameters
   CT delta_t = CT(path_parameter.max_delta_t,0);

   CT* tmp_t = (CT *)malloc(sizeof(CT));
   CT* tmp_t_last = (CT *)malloc(sizeof(CT));
   *tmp_t_last = cpu_t;

   int n_success = 0;

   while(tmp_t_last->real < T1(1)) 
   {
      if(debug)
      {
         std::cout << "n_point = " << n_point << ", n_step = "
                   << n_step << std::endl;
      }
      if(delta_t.real + tmp_t_last->real < 1) 
      {
         *tmp_t = *tmp_t_last + delta_t;
      }
      else
      {
         *tmp_t = CT(1,0);
      }
      if(debug)
      {
         std::cout << "delta_t = " << delta_t;
         std::cout << "tmp_t   = " << *tmp_t;
      }
      bool end_range = false;

      if(tmp_t->real>0.9)
      {
         end_range = true;
      }
      if(inverse == 0)
      {
         workspace.update_t_value(*tmp_t);
      }
      else
      {
         workspace.update_t_value_inverse(*tmp_t);
      }
      int n_predictor = min(workspace.n_predictor, n_point);

      if(debug)
      {
         std::cout << "n_predictor   = " << n_predictor << std::endl;
      }

      // int BS_pred = 32;
      // int nBS_pred = (inst.dim-1)/BS_pred+1;
      // std::cout << "workspace.x_t_idx = " << workspace.x_t_idx << std::endl;

      predict_newton_kernel<<<inst.predict_grid, inst.predict_BS>>>
         (workspace.x_array, workspace.t_array, n_predictor, inst.dim,
          workspace.x_t_idx_mult, workspace.workspace_size);

      if(debug)
      {
         std::cout << "Predict X:" << std::endl;
         workspace.print_x();

         /*std::cout << "X Array:" << std::endl;
	   workspace.print_x_array();*/
      }
      /*
       if(Record)
       {
          hom.path_data_gpu.add_step_empty();
          hom.path_data_gpu.update_step_t(delta_t, *tmp_t);
       }
       */
      bool newton_success = newton_single
         (workspace, inst, path_parameter, end_range);

      if(newton_success == 1)
      {
         if(debug)
         {
            std::cout << "---------- success -----------"<< std::endl;
	 }
         n_point++;
         workspace.update_x_t_idx();
         *tmp_t_last = *tmp_t;
         n_success++;
      }
      else
      {
         delta_t.real = delta_t.real*path_parameter.step_decrease;
         if(debug)
         {
            std::cout << "Decrease delta_t = " << delta_t << std::endl;
         }
         // std::cout << "      tmp_t_last = " << *tmp_t_last << std::endl;
         if(delta_t.real < path_parameter.min_delta_t)
         {
            break;
         }
         n_success = 0;
      }
      if(n_success > 1)
      {
         delta_t.real = delta_t.real*path_parameter.step_increase;;
         if(delta_t.real > path_parameter.max_delta_t)
         {
            delta_t.real = path_parameter.max_delta_t;
         }
         if(debug)
         {
            std::cout << "Increase delta_t = " << delta_t << std::endl;
         }
      }
      T1 max_delta_t_real;
      // std::cout << "tmp_t->real = " << tmp_t_last->real << std::endl;
      if(tmp_t_last->real > 0.9)
      {
         max_delta_t_real = 1E-2;
      }
      else
      {
         max_delta_t_real = path_parameter.max_delta_t;
      }
      if(delta_t.real > max_delta_t_real)
      {
         delta_t.real = max_delta_t_real;
      }
      n_step++;
      if(n_step >= path_parameter.max_step) 
      {
         break;
      }
   }
   bool success = 0;
   if(verbose > 0)
   {
      std::cout << "----------- Path Tracking Report ----------" << std::endl;
   }
   if(tmp_t_last->real == 1)
   {
      success = 1;
      if(verbose > 0) std::cout << "Success" << std::endl;
   }
   else
   {
      if(verbose > 0) std::cout << "Fail" << std::endl;
   }
   inst.n_step_GPU = n_step;
   inst.n_point_GPU = n_point;

   return success;
}

bool GPU_Path
 ( CPUInstHom& hom, Parameter path_parameter, CT* cpu_sol0, CT cpu_t, 
   CT*& x_gpu, int n_path, int inverse, int verbose )
{
	cuda_set();
   if(verbose > 0)
   {
      std::cout << "n_path = " << n_path << std::endl;
   }
   GPUInst inst(hom,n_path);
   // std::cout << "mon_pos_block_size = "
   //           << hom.CPU_inst_hom_block.mon_pos_block_size
   //           << " inst.n_mon_level[0] = " << inst.n_mon_level[0]
   //           << " mon_pos_size = " << mon_pos_size << std::endl;
   GPUWorkspace workspace(inst.mon_pos_size, inst.n_coef, inst.n_constant, 
      inst.n_eq, inst.dim, path_parameter.n_predictor, inst.alpha);
   workspace.update_x_t(cpu_sol0, cpu_t);

   struct timeval start, end;
   long seconds, useconds;
   gettimeofday(&start, NULL);

   bool success = path_single(workspace,inst,path_parameter,cpu_t,n_path,
                              inverse,verbose);
   x_gpu = workspace.get_x_last();

   gettimeofday(&end, NULL);

   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   double timeMS_Path_GPU = seconds*1000 + useconds/1000.0;
   double timeSec_Path_GPU = timeMS_Path_GPU/1000;
   if(verbose > 0)
   {
      cout << "Path GPU Test MS   Time: "<< timeMS_Path_GPU << endl;
      cout << "Path GPU Test      Time: "<< timeSec_Path_GPU << endl;
      cout << "Path GPU Step     Count: "<< inst.n_step_GPU << endl;
      cout << "Path GPU Point    Count: "<< inst.n_point_GPU << endl;
      cout << "Path GPU Eval     Count: "<< inst.n_eval_GPU << endl;
      cout << "Path GPU MGS      Count: "<< inst.n_mgs_GPU << endl;
   }
   hom.timeSec_Path_GPU = timeSec_Path_GPU;
   hom.n_step_GPU = inst.n_step_GPU;
   hom.n_point_GPU = inst.n_point_GPU;
   hom.n_eval_GPU = inst.n_eval_GPU;
   hom.n_mgs_GPU = inst.n_mgs_GPU;

   return success;
}

#endif /* PATH_GPU_CU_ */
