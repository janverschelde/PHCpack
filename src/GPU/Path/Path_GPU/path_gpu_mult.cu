__global__ void update_t_kernel
 ( GT* t_mult, GT* t_last_mult, GT* delta_t_mult, GT* t_array_mult,
   int* path_success, int* newton_success, int* end_range, int* n_success,
   int* n_point_mult, int* x_t_idx_mult, int workspace_size, int n_path,
   double step_increase, double step_decrease, int n_predictor )
{
   int path_idx = threadIdx.x+blockIdx.x*blockDim.x;
   if(path_idx<n_path)
   {
      if(path_success[path_idx] == 0)
      {
         // Check delta_t
         if(newton_success[path_idx] == 1)
         {
            if(t_mult[path_idx].real >= 1.0)
            {
               path_success[path_idx] = 1;
            }
            else
            {
               // n_point_mult[path_idx]++;
               t_last_mult[path_idx] = t_mult[path_idx];
               n_success[path_idx]++;
               if(n_success[path_idx] > 1) 
               {
                  // step increase
                  delta_t_mult[path_idx].real
                     = delta_t_mult[path_idx].real*step_increase;
               }
               T max_delta_t;
               if(t_last_mult[path_idx].real > 0.9)
               {
                  max_delta_t = GT(1E-2,0.0).real;
               }
               else
               {
                  max_delta_t = GT(1E-1,0.0).real;;
               }
               if(delta_t_mult[path_idx].real > max_delta_t) 
               {
                  delta_t_mult[path_idx].real = max_delta_t;
               }
            }
         }
         else
         {
            // std::cout << " fail" << std::endl;
            delta_t_mult[path_idx].real
               = delta_t_mult[path_idx].real*step_decrease;
            if(delta_t_mult[path_idx].real < 1E-7)
            {
               path_success[path_idx] = -1;
            }
            n_success[path_idx] = 0;
         }
      }
      if(path_success[path_idx] == 0)
      {
         // Check t
         t_mult[path_idx] = t_last_mult[path_idx] + delta_t_mult[path_idx];
         if(t_mult[path_idx].real > 1.0)
         {
            t_mult[path_idx].real = GT(1.0,0.0).real;
         }
         if(t_mult[path_idx].real > 0.9)
         {
            end_range[path_idx] = 1;
         }
         else
         {
            end_range[path_idx] = 0;
         }
         // GT* t = t_array + workspace_size*path_idx + x_t_idx_mult[path_idx];
         // *t = t_mult[path_idx];
         GT* t = t_array_mult + (n_predictor+1)*path_idx
               + x_t_idx_mult[path_idx];
         *t = t_mult[path_idx];
         // GT* one_minor_t = one_minor_t_mult + workspace_size*path_idx;
         // *one_minor_t = (*alpha)*(GT(1.0,0) - t_mult[path_idx]);
      }
   }
}

int path_mult
 ( GPUWorkspace& workspace, GPUInst& inst, Parameter path_parameter,
   CT* cpu_t, int inverse = 0, int verbose = 0 )
{
   bool debug = false; //debug = true;
   int path_idx_test = 0;

   int n_point = 1;
   int n_step = 0;
   int n_path = workspace.n_path;

   /*
    for(int path_idx=0; path_idx<n_path; path_idx++)
    {
       workspace.end_range_host[path_idx] = 0;
    }
    cudaMemcpy(workspace.end_range,workspace.end_range_host,
       n_path*sizeof(int),cudaMemcpyHostToDevice);
    // Parameters
    CT* delta_t = (CT *)malloc(n_path*sizeof(CT));
    CT* tmp_t = (CT *)malloc(n_path*sizeof(CT));
    CT* tmp_t_last = (CT *)malloc(n_path*sizeof(CT));
    for(int path_idx=0; path_idx<n_path; path_idx++)
    {
       delta_t[path_idx] = CT(path_parameter.max_delta_t,0);
       tmp_t_last[path_idx] = cpu_t[path_idx];
       cpu_t[path_idx] += delta_t[path_idx];
    }
    for(int path_idx=0; path_idx<n_path; path_idx++)
    {
       comp1_qd2gqd(tmp_t_last+path_idx, workspace.t_last_mult_host+path_idx);
       comp1_qd2gqd(cpu_t+path_idx, workspace.t_mult_host+path_idx);
       comp1_qd2gqd(delta_t+path_idx, workspace.delta_t_mult_host+path_idx);
    }
    cudaMemcpy(workspace.t_last_mult,workspace.t_last_mult_host,
       n_path*sizeof(GT),cudaMemcpyHostToDevice);
    cudaMemcpy(workspace.t_mult, workspace.t_mult_host, n_path*sizeof(GT),
       cudaMemcpyHostToDevice);
    cudaMemcpy(workspace.delta_t_mult,workspace.delta_t_mult_host,
       n_path*sizeof(GT),cudaMemcpyHostToDevice);
    for(int path_idx=0; path_idx<n_path; path_idx++)
    {
       workspace.path_success_host[path_idx] = 0;
       workspace.n_success_host[path_idx] = 0;
    }
    cudaMemcpy(workspace.path_success,workspace.path_success_host,
               n_path*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(workspace.n_success,workspace.n_success_host,
               n_path*sizeof(int),cudaMemcpyHostToDevice);
   */
   dim3 update_t_grid = get_grid(n_path,inst.predict_BS,1);

   path_mult_init_kernel
      <<<update_t_grid,inst.predict_BS>>>
      (workspace.t_mult,workspace.t_last_mult,workspace.delta_t_mult,
       workspace.t_array,workspace.one_minor_t,workspace.alpha_gpu,
       workspace.path_success,workspace.newton_success,workspace.end_range,
       workspace.n_success,workspace.max_f_val_last_gpu,
       workspace.n_point_mult,workspace.x_t_idx_mult,workspace.workspace_size,
       workspace.path_idx, n_path,workspace.matrix,inst.n_sum_zero,
       inst.sum_zeros, workspace.n_predictor);

   path_mult_x_init_kernel
      <<<n_path, inst.dim>>>
      (workspace.x_mult, workspace.x_array,workspace.workspace_size,
       inst.dim, workspace.n_predictor);

   cudaMemcpy(workspace.path_success_host,workspace.path_success,
              n_path*sizeof(int),cudaMemcpyDeviceToHost);

   int n_path_continuous = n_path;
   int n_path_fail = 0;

   while(true)
   {
      if(debug)
      {
         std::cout << "n_point = " << n_point
                   << ", n_step = " << n_step
                   << " n_path_continuous = " << n_path_continuous
                   << std::endl;
      }
      // std::cout << n_step << " " 
      //  << "workspace.path_success_host[path_idx_test]" 
      // << workspace.path_success_host[path_idx_test] << std::endl;
      // int n_predictor = min(workspace.n_predictor, n_point);
      // std::cout << "n_predictor   = " << n_predictor << std::endl;
      // int BS_pred = 32;
      // int nBS_pred = (inst.dim-1)/BS_pred+1;
      // std::cout << "workspace.x_t_idx = " << workspace.x_t_idx << std::endl;
      dim3 predict_grid = get_grid(inst.dim,inst.predict_BS,n_path_continuous);
      predict_newton_kernel
         <<<predict_grid, inst.predict_BS>>>
         (workspace.x_array,workspace.t_array,workspace.n_predictor,
          inst.dim, workspace.x_t_idx_mult, workspace.workspace_size,
          workspace.n_point_mult, workspace.path_idx);

      if(debug)
      {
         std::cout << "t:" << std::endl;
         workspace.print_t_mult(path_idx_test);
         std::cout << "delta:" << std::endl;
         workspace.print_delta_t_mult(path_idx_test);
         // std::cout << "Predict X:" << std::endl;
         // workspace.print_x_mult(path_idx_test);
         /*
          std::cout << "T Array:" << std::endl;
          workspace.print_t_array();
          std::cout << "X Array:" << std::endl;
          workspace.print_x_array();
         */
      }
      newton_align(workspace, inst, path_parameter);

      update_t_kernel
         <<<update_t_grid, inst.predict_BS>>>
         (workspace.t_mult,workspace.t_last_mult,workspace.delta_t_mult,
          workspace.t_array,workspace.path_success,workspace.newton_success,
          workspace.end_range, workspace.n_success,workspace.n_point_mult,
          workspace.x_t_idx_mult, workspace.workspace_size, n_path,
          path_parameter.step_increase,path_parameter.step_decrease,
          workspace.n_predictor);

      // update kernel
      /*
       for(int path_idx=0; path_idx<n_path; path_idx++){ }
       cudaMemcpy(workspace.end_range,workspace.end_range_host,
                  n_path*sizeof(int),cudaMemcpyHostToDevice);
      */
      /*
       if(debug && path_idx==path_idx_test)
       {
          // std::cout << "Decrease delta_t = " << delta_t[path_idx]
          //           << std::endl;
          // std::cout << "      tmp_t_last = " << *tmp_t_last << std::endl;
       }
       if(debug && path_idx==path_idx_test)
       {
          std::cout << path_idx << " " << delta_t[path_idx].real 
                    << " " << tmp_t[path_idx].real << std::endl;
       }
      */
      /*
       if(inverse == 0)
       {
          workspace.update_t_value_mult(tmp_t);
       }
       else
       {
          // workspace.update_t_value_inverse(*tmp_t);
       }
      */
      cudaMemcpy(workspace.path_success_host,workspace.path_success,
                 n_path*sizeof(int),cudaMemcpyDeviceToHost);

      int n_path_success = 0;
      n_path_fail = 0;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         if(workspace.path_success_host[path_idx] == 1)
         {
            n_path_success++;
         }
         else if(workspace.path_success_host[path_idx] == -1)
         {
            n_path_fail++;
         }
      }
      if(n_path_success + n_path_fail == workspace.n_path)
      {
         if(verbose > 0)
         {
            std::cout << "ALL PATH Finished!!!" << std::endl;
            std::cout << "Success " << n_path_success << std::endl;
            std::cout << "Fail    " << n_path_fail << std::endl;
         }
         break;
      }
      cudaMemcpy(workspace.path_success,workspace.path_success_host,
                 n_path*sizeof(int),cudaMemcpyHostToDevice);
      n_step++;
      if(n_step >= path_parameter.max_step) 
      {
         break;
      }
      n_path_continuous = 0;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // std::cout << path_idx << " "
         //           << workspace.newton_success_host[path_idx] << std::endl;
         if(workspace.path_success_host[path_idx] == 0)
         {
            workspace.path_idx_host[n_path_continuous] = path_idx;
            n_path_continuous += 1;
         }
      }
      workspace.n_path_continuous = n_path_continuous;
      cudaMemcpy(workspace.path_idx,workspace.path_idx_host,
                 n_path*sizeof(int),cudaMemcpyHostToDevice);
   }
   path_mult_x_finish_kernel
      <<<n_path, inst.dim>>>
         (workspace.x_mult,workspace.x_array, workspace.x_t_idx_mult,
          workspace.workspace_size, inst.dim, workspace.n_predictor);

   int n_path_success = 0;
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      if(workspace.path_success_host[path_idx] == 1)
      {
         n_path_success++;
      }
   }
   if(verbose > 0)
   {
      std::cout << "-------------- Path Tracking Report ---------------"
                << std::endl;
      std::cout << "Success " << n_path_success << std::endl;
      std::cout << "Fail "  << workspace.n_path-n_path_success << std::endl;
   }
   inst.n_step_GPU = n_step;
   inst.n_point_GPU = n_point;

   return n_path_success;
}

bool* GPU_Path_mult
 ( CPUInstHom& hom, Parameter path_parameter, CT* cpu_sol0, CT* cpu_t,
   CT**& x_gpu, int n_path, int inverse, int verbose )
{
   cuda_set();

   if(verbose > 0) std::cout << "n_path = " << n_path << std::endl;

   GPUInst inst(hom, n_path);
   GPUWorkspace workspace(inst.mon_pos_size,inst.n_coef,inst.n_constant,
      inst.n_eq,inst.dim,path_parameter.n_predictor,inst.alpha,
      inst.base_table_size,n_path);

   if(verbose > 0)
   {
      std::cout << "workspace.n_path = " << workspace.n_path << std::endl;
   }

   /*
    int* x_t_idx = (int *)malloc(n_path*sizeof(int));
    for(int sys_idx=0; sys_idx<n_path; sys_idx++)
    {
       x_t_idx[sys_idx] = 1;
    }
    workspace.update_x_t_value_mult(cpu_sol0, cpu_t);
    workspace.update_x_t_idx_all(x_t_idx);
    delete[] x_t_idx;
   */
   workspace.update_x_mult_horizontal(cpu_sol0);
   // workspace.print_x_last_mult();

   struct timeval start, end;
   long seconds, useconds;
   gettimeofday(&start, NULL);

   path_mult(workspace, inst, path_parameter, cpu_t, inverse, verbose);
   x_gpu = workspace.get_mult_x_horizontal();

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

   bool* path_success = (bool *)malloc(n_path*sizeof(bool));
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      if(workspace.path_success_host[path_idx] == 1)
      {
         path_success[path_idx] = true;
      }
      else
      {
         path_success[path_idx] = false;
      }
   }
   return path_success;
}
