/* path_test.cpp, created on Feb 8, 2015 by yxc with edits by jv */

#include "path_test.h"

T1 path_test_reverse
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT t, int device_option)
{
   std::cout << "--------- Path Tracking Reverse Test ----------" << std::endl;

   bool cpu_test = 0;
   bool gpu_test = 0;
   if(device_option == 1)
   {
      cpu_test = 1;
      gpu_test = 0;
   }
   else if(device_option == 2)
   {
      cpu_test = 0;
      gpu_test = 1;
   }
   else
   {
      std::cout << "Device_option Invalid. Choose from the following:"
                << std::endl
                << "  1. CPU" << std::endl
                << "  2. GPU" << std::endl
                << "Your device_option = " << device_option << std::endl;
   }

   double timeSec_Predict = 0;
   double timeSec_Eval = 0;
   double timeSec_MGS = 0;
   double timeSec = 0;

   CT* x_cpu;
   CT* x_cpu_target;

   if(cpu_test == 1)
   {
      clock_t begin = clock();
      bool success = path_tracker
         (workspace_cpu,cpu_inst_hom,path_parameter,
          timeSec_Predict,timeSec_Eval,timeSec_MGS);
      clock_t end = clock();
      timeSec += (end - begin) / static_cast<double>( CLOCKS_PER_SEC );
      x_cpu_target = workspace_cpu.x_last;
   }

   if(gpu_test == 1)
   {
      bool success = GPU_Path(cpu_inst_hom,path_parameter,sol0,
                              t,x_cpu_target);
   }
   // Path Tracking Reverse Test
   if(cpu_test == 1)
   {
      t = CT(0.0,0.0);
      workspace_cpu.init_x_t_idx();
      workspace_cpu.update_x_t(x_cpu_target, t);
      bool success = path_tracker(workspace_cpu,cpu_inst_hom,path_parameter,
                                  timeSec_Predict,timeSec_Eval,timeSec_MGS,1);
      x_cpu = workspace_cpu.x_last;
   }
   if(gpu_test == 1)
   {
      bool success = GPU_Path(cpu_inst_hom,path_parameter,x_cpu_target,
                              t,x_cpu,1);
   }
   // delete[] x_cpu;
   /*for(int i=0; i<cpu_inst_hom.dim; i++) 
     {
        std::cout << i << " " << x_cpu[i];
     }*/

   std::cout << "--------- Start solution vs Reverse solution ----------"
             << std::endl;
   T1 err = err_check_workspace(sol0, x_cpu, cpu_inst_hom.dim);
   cout << "Path CPU Predict   Time: "<< timeSec_Predict << endl;
   cout << "Path CPU Eval      Time: "<< timeSec_Eval << endl;
   cout << "Path CPU MGS       Time: "<< timeSec_MGS << endl;

   return err;
}

bool path_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* sol0, CT*& sol_new, CT t, PolySys& Target_Sys,
   int device_option )
{
   // std::cout << "--------- Path Tracking Test ----------" << std::endl;

   bool Record = false;
   Record = true;

   double timeSec_Path_CPU = 0;
   double timeMS_Path_CPU = 0;
   double timeSec_Predict_CPU = 0;
   double timeSec_Eval_CPU = 0;
   double timeSec_MGS_CPU = 0;
   CT* x_cpu = NULL;
   CT* x_gpu = NULL;

   bool cpu_test = 0;
   bool gpu_test = 0;
   bool cpu_success = 0;
   bool gpu_success = 0;

   T1 err = 0;

   if(device_option == 0)
   {
      cpu_test = 1;
      gpu_test = 1;
   }
   else if(device_option == 1)
   {
      cpu_test = 1;
      gpu_test = 0;
   }
   else if(device_option == 2)
   {
      cpu_test = 0;
      gpu_test = 1;
   }
   else
   {
      std::cout << "Device_option Invalid. Choose from the following:"
                << std::endl
                << "  0. CPU and GPU" << std::endl
                << "  1. CPU" << std::endl
                << "  2. GPU" << std::endl
                << "Your device_option = " << device_option << std::endl;
   }
   if(cpu_test == 1)
   {
      if(Record == true)
      {
         cpu_inst_hom.path_data.clear();
      }
      struct timeval start, end;
      long seconds, useconds;
      gettimeofday(&start, NULL);
      cpu_success = path_tracker(workspace_cpu, cpu_inst_hom, path_parameter,
         timeSec_Predict_CPU,timeSec_Eval_CPU,timeSec_MGS_CPU);
      x_cpu = workspace_cpu.x_last;
      gettimeofday(&end, NULL);
      if(Record == true)
      {
         cpu_inst_hom.path_data.print_phc();
      }
      seconds = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      timeMS_Path_CPU = seconds*1000 + useconds/1000.0;
      cpu_inst_hom.timeSec_Path_CPU = timeMS_Path_CPU/1000;
   }
   if(gpu_test == 1)
   {
      gpu_success = GPU_Path(cpu_inst_hom,path_parameter,sol0,t,x_gpu);
      /*if(gpu_success == 1)
        {
           // cyclic solution 0 pieri solution 2/3
           string file_name = sol_filename(cpu_inst_hom.dim, 0);
           write_complex_array(file_name, x_gpu, cpu_inst_hom.dim);
        }*/
   }
   if(cpu_test == 1)
   {
      // Print Time and Error Report
      cout << "Path CPU Path MS   Time: "<< timeMS_Path_CPU << endl;
      cout << "Path CPU Path      Time: "
           << cpu_inst_hom.timeSec_Path_CPU << endl;
      // cout << "Path CPU Predict   Time: "<< timeSec_Predict_CPU << endl;
      // cout << "Path CPU Eval      Time: "<< timeSec_Eval_CPU << endl;
      // cout << "Path CPU MGS       Time: "<< timeSec_MGS_CPU << endl;
      cout << "Path CPU Step     Count: "<< cpu_inst_hom.n_step_CPU << endl;
      cout << "Path CPU Point    Count: "<< cpu_inst_hom.n_point_CPU << endl;
      cout << "Path CPU Eval     Count: "<< cpu_inst_hom.n_eval_CPU << endl;
      cout << "Path CPU MGS      Count: "<< cpu_inst_hom.n_mgs_CPU << endl;
      cpu_inst_hom.timeSec_Path_CPU = timeSec_Path_CPU;
   }
   if(cpu_test == 1 && gpu_test == 1)
   {
      std::cout << "--------- Path Tracking Error CPU vs GPU ----------"
                << std::endl;
      err = err_check_workspace(x_cpu, x_gpu, cpu_inst_hom.dim);
      std::cout << "err = " << err << std::endl;
   }
   if(cpu_test == 1)
   {
      std::cout << "CPU Solution" << std::endl;
      T1 max_x = 0;
      for(int i=0; i<cpu_inst_hom.dim; i++)
      {
         // std::cout << i << "  " << x_cpu[i];
         if(abs(x_cpu[i].real) > max_x)
         {
            max_x = abs(x_cpu[i].real);
         }
         if(abs(x_cpu[i].imag) > max_x)
         {
            max_x = abs(x_cpu[i].imag);
         }
      }
      std::cout << "Max abs(x): " << max_x << std::endl;
   }
   if(gpu_test == 1)
   {
      std::cout << "GPU Solution Max Abs ";
      T1 max_x = 0;
      for(int i=0; i<cpu_inst_hom.dim; i++)
      {
         // std::cout << i << "  " << x_gpu[i];
         if(abs(x_gpu[i].real) > max_x)
         {
            max_x = abs(x_gpu[i].real);
         }
         if(abs(x_gpu[i].imag) > max_x)
         {
            max_x = abs(x_gpu[i].imag);
         }
      }
      std::cout << max_x << std::endl;
   }
   if(cpu_test == 1)
   {
      std::cout << "CPU Path: ";
      if(cpu_success == 1)
      {
         std::cout << "Success!" << std::endl;
         int res_check = 1;
         if(res_check == 1)
         {
            int n_eq = cpu_inst_hom.n_eq;
            int dim = cpu_inst_hom.dim;
            CT* workspace_classic = new CT[n_eq*(dim+1)];
            CT* f_val = workspace_classic;
            CT** deri_val = new CT*[n_eq];
            CT* tmp_workspace = workspace_classic + n_eq;
            for(int i=0; i<n_eq; i++)
            {
               deri_val[i] = tmp_workspace;
               tmp_workspace += dim;
            }
            Target_Sys.eval(x_cpu,f_val,deri_val);
            T1 max_residual = 0;
            for(int i=0; i<n_eq; i++)
            {
               // std::cout << i << " " << f_val[i];
               T1 tmp_residual = f_val[i].real*f_val[i].real 
                               + f_val[i].imag*f_val[i].imag;
               if(tmp_residual > max_residual)
               {
                  max_residual = tmp_residual;
               }
            }
            max_residual = sqrt(max_residual);
            cpu_inst_hom.max_residual = max_residual;
            std::cout << "CPU Residual Check: "
                      << max_residual << std::endl;
            delete[] deri_val;
            delete[] workspace_classic;
         }
      }
      else
      {
         std::cout << "Fail!" << std::endl;
      }
      sol_new = x_cpu;
   }
   if(gpu_test == 1)
   {
      std::cout << "GPU Path: ";
      if(gpu_success == 1)
      {
         std::cout << "Success!" << std::endl;
         int res_check = 0;
         if(res_check == 1)
         {
            int n_eq = cpu_inst_hom.n_eq;
            int dim = cpu_inst_hom.dim;
            CT* workspace_classic = new CT[n_eq*(dim+1)];
            CT* f_val = workspace_classic;
            CT** deri_val = new CT*[n_eq];
            CT* tmp_workspace = workspace_classic + n_eq;
            for(int i=0; i<n_eq; i++)
            {
               deri_val[i] = tmp_workspace;
               tmp_workspace += dim;
            }
            Target_Sys.eval(x_gpu,f_val,deri_val);
            T1 max_residual = 0;
            for(int i=0; i<n_eq; i++)
            {
               T1 tmp_residual = f_val[i].real*f_val[i].real
                               + f_val[i].imag*f_val[i].imag;
               if(tmp_residual > max_residual)
               {
                  max_residual = tmp_residual;
               }
            }
            max_residual = sqrt(max_residual);
            std::cout << "GPU Residual Check: "
                      << max_residual << std::endl;
            delete[] deri_val;
            delete[] workspace_classic;
         }
      }
      else
      {
         std::cout << "Fail!" << std::endl;
      }
      sol_new = x_gpu;
   }
   // if cpu or gpu tested, it has to be a success
   // no test or success
   cpu_inst_hom.success_CPU = cpu_success;
   cpu_inst_hom.success_GPU = gpu_success;
   bool success = (!cpu_test||cpu_success) && (!gpu_test||gpu_success);
   return success;
}
