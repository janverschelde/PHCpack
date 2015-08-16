/* eval_test.cpp, created on: Feb 8, 2015 by yxc with edits by jv */

#include "eval_test.h"

T1 eval_compare 
 ( const CPUInstHom& cpu_inst_hom, CT* gpu_workspace, CT* gpu_matrix,
   const CT* cpu_workspace, const CT* cpu_matrix )
{
   int n_coef = cpu_inst_hom.CPU_inst_hom_coef.n_coef;
   int n_workspace_size = n_coef + cpu_inst_hom.CPU_inst_hom_mon.mon_pos_size;
   int n_workspace = n_coef+cpu_inst_hom.CPU_inst_hom_mon.mon_pos_size;

   std::cout << "----- Coef Check CPU vs GPU ----"<< n_coef << std::endl;
   err_check_workspace(cpu_workspace, gpu_workspace, n_coef);

   if(MON_EVAL_METHOD != 1)
   {
      std::cout << "----- Workspace Check CPU vs GPU ----" << std::endl;
      err_check_workspace(cpu_workspace, gpu_workspace, n_workspace);
      /*std::cout << "n_workspace = " << n_workspace << std::endl;
        for(int i=10; i<20; i++)
        {
           std::cout << i << std::endl 
                    << cpu_workspace[i] << gpu_workspace[i];
        }
      */
   }
   std::cout << "----- Jacobian and Fun Check CPU vs GPU ----"
             << cpu_inst_hom.n_eq*(cpu_inst_hom.dim+1) << std::endl;
   T1 err = err_check_workspace(cpu_matrix,gpu_matrix,
      cpu_inst_hom.n_eq*(cpu_inst_hom.dim+1),
      cpu_inst_hom.n_eq*(cpu_inst_hom.dim+1));

   return err;
}

T1 eval_test_classic
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, CT* sol0, CT t,
   PolySys& Classic_Sys, int n_path )
{
   struct timeval start, end;
   long seconds, useconds;
   double timeMS_classic;
   double timeMS_cpu;
   double timeMS_gpu;

   int n_eq = cpu_inst_hom.n_eq;
   int dim = cpu_inst_hom.dim;

   if(n_path<=0)
   {
      std::cout << "Default number of path" << std::endl;
      n_path = 1000;
   }
   int n_predictor = workspace_cpu.n_predictor;
   std::cout << "n_path = " << n_path << std::endl;
   CT* sol = new CT[n_path*dim*(n_predictor+1)];
   CT* sol_tmp = sol;
   for(int sol_idx=0; sol_idx<n_path; sol_idx++)
   {
      for(int pred_idx=0; pred_idx<n_predictor+1; pred_idx++)
      {
         for(int x_idx=0; x_idx<dim; x_idx++)
         {
            int r = rand();
            T1 tmp = T1(r);
            // sol_tmp[x_idx] = CT(sin(tmp),cos(tmp));
            sol_tmp[x_idx] = CT(x_idx+1,0.0);
            // sol_tmp[x_idx] = CT(1,0.0);
         }
         sol_tmp += dim;
      }
   }
   CT* t_mult = new CT[n_path*(n_predictor+1)];
   for(int sol_idx=0; sol_idx<n_path*(n_predictor+1); sol_idx++)
   {
      double r = 1.0*rand()/RAND_MAX;
      // t_mult[sol_idx] = CT(r,0.0);
      t_mult[sol_idx] = CT(1,0.0);
   }
   int* x_t_idx = new int[n_path];
   for(int sol_idx=0; sol_idx<n_path; sol_idx++)
   {
      x_t_idx[sol_idx] = rand()%(n_predictor+1);
   }
   std::cout << "----- CPU Evaluation ----" << std::endl;
   Workspace* workspace_cpu_all = new Workspace[n_path];
   for(int sol_idx=0; sol_idx<n_path; sol_idx++)
   {
      cpu_inst_hom.init_workspace(workspace_cpu_all[sol_idx]);
   }
   gettimeofday(&start, NULL);
   for(int sol_idx=0; sol_idx<n_path; sol_idx++)
   {
      CT* tmp_sol = sol+sol_idx*dim*(n_predictor+1)+dim*x_t_idx[sol_idx];
      CT* t_tmp = t_mult+sol_idx*(n_predictor+1)+x_t_idx[sol_idx];
      cpu_inst_hom.eval(workspace_cpu_all[sol_idx], tmp_sol, *t_tmp);
   }
   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   timeMS_cpu = seconds*1000 + useconds/1000.0;

   bool classic_check = false;
   if(classic_check)
   {
      std::cout << "----- Class Evaluation ----" << std::endl;
      CT* workspace_classic = new CT[n_path*n_eq*(dim+1)];
      CT** f_val = new CT*[n_path];
      CT* tmp_workspace = workspace_classic;
      CT*** deri_val = new CT**[n_path];
      CT** deri_space = new CT*[n_path*n_eq];

      for(int sol_idx=0; sol_idx<n_path; sol_idx++)
      {
         f_val[sol_idx] = tmp_workspace;
         tmp_workspace += n_eq;
         deri_val[sol_idx] = deri_space + sol_idx*n_eq;
         for(int i=0; i<n_eq; i++)
         {
            deri_val[sol_idx][i] = tmp_workspace;
            tmp_workspace += dim;
         }
      }
      gettimeofday(&start, NULL);
      for(int sol_idx=0; sol_idx<n_path; sol_idx++)
      {
         CT* tmp_sol = sol+sol_idx*dim*(n_predictor+1)+dim*x_t_idx[sol_idx];
         Classic_Sys.eval(tmp_sol, f_val[sol_idx], deri_val[sol_idx]);
      }
      gettimeofday(&end, NULL);
      seconds  = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      timeMS_classic = seconds*1000 + useconds/1000.0;

      // Check two CPU method
      std::cout << "----- Classic Evaluation Check ----" << std::endl;
      for(int sol_idx=0; sol_idx<n_path; sol_idx++)
      {
         err_check_class_workspace(deri_val[sol_idx],f_val[sol_idx],
            workspace_cpu_all[sol_idx].matrix, n_eq, dim);
      }
      delete[] workspace_classic;
      delete[] f_val;
      delete[] deri_val;
      delete[] deri_space;
   }
   std::cout << "----- GPU Evaluation ----" << std::endl;
   CT** gpu_workspace_all;
   CT** gpu_matrix_all;
   gettimeofday(&start, NULL);
   GPU_Eval(cpu_inst_hom,sol,t_mult,gpu_workspace_all,gpu_matrix_all,n_path,
      x_t_idx, n_predictor);
   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   timeMS_gpu = seconds*1000 + useconds/1000.0;
   std::cout << "----- CPU vs GPU Evaluation Check----" << std::endl;
   T1 err = 0;
   for(int sol_idx=0; sol_idx<n_path; sol_idx++)
   {
      // std::cout << "sol_idx = " << sol_idx << std::endl;
      T1 err_tmp = eval_compare(cpu_inst_hom,gpu_workspace_all[sol_idx],
         gpu_matrix_all[sol_idx],workspace_cpu_all[sol_idx].all,
         workspace_cpu_all[sol_idx].matrix);
      if(err_tmp > err)
      {
         err = err_tmp;
      }
      // std::cout << "err = " << err_tmp << std::endl;
   }
   delete[] x_t_idx;
   delete[] t_mult;
   delete[] sol;

   for(int sol_idx=0; sol_idx<n_path; sol_idx++)
   {
      delete[] gpu_workspace_all[sol_idx];
      delete[] gpu_matrix_all[sol_idx];
   }
   delete[] gpu_workspace_all;
   delete[] gpu_matrix_all;

   std::cout << "err = " << err << std::endl;
   std::cout << "Classic Eval time " << timeMS_classic << std::endl;
   std::cout << "CPU     Eval time " << timeMS_cpu << std::endl;
   std::cout << "GPU     Eval time " << timeMS_gpu << std::endl;

   return err;
}
