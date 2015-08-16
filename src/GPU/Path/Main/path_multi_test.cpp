/* multi_path_test.cpp, created on Feb 12, 2015 by yxc with edits by jv */

#include "path_multi_test.h"
#include "path_data.h"

class path_report
{
   public:

      bool success;
      CT t;
      double path_time;
      int n_step;
      int n_point;
      int n_eval;
      int n_mgs;
      T1 max_x;

      void init 
       ( bool success, CT t, double path_time, int n_step, int n_point,
         int n_eval, int n_mgs, int dim, CT* x_last )
      {
         this->success = success;
         this->t = t;
         this->path_time = path_time;
         this->n_step = n_step;
         this->n_point = n_point;
         this->n_eval = n_eval;
         this->n_mgs = n_mgs;
         max_x = 0.0;
         for(int i=0; i<dim; i++)
         {
            T1 x_tmp = x_last[i].real*x_last[i].real
                     + x_last[i].imag*x_last[i].imag;
            if(x_tmp > max_x)
            {
               max_x = x_tmp;
            }
         }
      }

      path_report()
      {
         success = 0;
         t = CT(0.0,0.0);
         path_time = 0;
         n_step = 0;
         n_point = 0;
         n_eval = 0;
         n_mgs = 0;
         max_x = 0.0;
      };

      path_report
       ( bool success, CT t, double path_time, int n_step, int n_point,
         int n_eval, int n_mgs, int dim, CT* x_last )
      {
         init(success,t,path_time,n_step,n_point,n_eval,n_mgs,dim,x_last);
      }

      void print()
      {
         std::cout << "success   = " << success << std::endl;
         std::cout << "t         = " << t.real << std::endl;
         std::cout << "path_time = " << path_time << std::endl;
         std::cout << "n_step    = " << n_step << std::endl;
         std::cout << "n_point   = " << n_point << std::endl;
         std::cout << "n_eval    = " << n_eval << std::endl;
         std::cout << "n_mgs     = " << n_mgs << std::endl;
      }

      void print_line()
      {
         std::cout << std::fixed
                   << std::setw(2) << success
                   << std::setw(10) << std::setprecision(5) << path_time
                   << std::setw(6) << n_step
                   << std::setw(5) << n_point
                   << std::setw(8) << n_eval
                   << std::setw(8) << n_mgs
                   << std::setw(25) << std::setprecision(15) << t.real
                   << std::setw(35) << std::setprecision(10) << max_x
                   << std::endl;
      }
};

class path_report_set
{
   public:

     int n_path;
     path_report* reports;

     path_report_set(int n_path)
     {
        this->n_path = n_path;
        reports = new path_report[n_path];
     }

     ~path_report_set()
     {
        delete[] reports;
     }

     void load_report
      ( int i, CPUInstHom& cpu_inst_hom, Workspace& workspace_cpu )
     {
        reports[i].init(cpu_inst_hom.success_CPU,cpu_inst_hom.t_CPU,
           cpu_inst_hom.timeSec_Path_CPU,cpu_inst_hom.n_step_CPU,
           cpu_inst_hom.n_point_CPU,cpu_inst_hom.n_eval_CPU,
           cpu_inst_hom.n_mgs_CPU,cpu_inst_hom.dim, workspace_cpu.x_last);
     }

     void path_number_report()
     {
        int n_path_success = 0;
        int n_path_fail = 0;
        for(int i=0; i<n_path; i++)
        {
           if(reports[i].success==true)
           {
              n_path_success++;
           }
           else
           {
              n_path_fail++;
           }
        }
        std::cout << std::fixed
                  << "Number of Path Success : " << setw(8) << n_path_success
                  << std::endl
                  << "Number of Path Fail    : " << setw(8) << n_path_fail
                  << std::endl;
     }

     void print_fail_path_reports()
     {
        for(int i=0; i<n_path; i++)
        {
           if(reports[i].success == false)
           {
              reports[i].print_line();
           }
        }
     }
};

void compare_sol_with_phc(int dim, PolySolSet& sols)
{
   std::ostringstream cyclic_filename_sol;
   // cyclic_filename_sol << "../Problems/MultiPath/cyclic"
   //                     << dim << ".path_sol";
   cyclic_filename_sol << "../Problems/MultiPath/game8two.path_sol";

   string filename_sol = cyclic_filename_sol.str();
   std::cout << filename_sol << std::endl;
   ifstream myfile_sol(filename_sol.c_str());

   if(myfile_sol.is_open() == false)
   {
      std::cout << "Can't open sol file..." << std::endl;
      return;
   }

   PolySolSet sols_phc(myfile_sol);

   // sols_phc.print_short();
   // sols.print_short();

   int n_sol_diff = 0;
   int start_idx = 0;
   // start_idx = 80;
   for(int i=0; i<sols.n_sol; i++)
   {
      // (sols_phc.sols[start_idx+i])->print();
      // (sols.sols[i])->print();
      if(sols_phc.sols[start_idx+i]->info != "success" 
         or sols.sols[i]->info != "success")
      {
         if(sols_phc.sols[start_idx+i]->info != sols.sols[i]->info)
         {
            std::cout << i << " " << sols_phc.sols[start_idx+i]->info
                      << " " << sols.sols[i]->info << std::endl;
            n_sol_diff++;
         }
      }
      else
      {
         bool check_sol
          = ((*(sols_phc.sols[start_idx+i])) == (*(sols.sols[i])));
         if(check_sol == false)
         {
            std::cout << i << " Different Result" 
            // << std::endl;
            // std::cout << sols_phc.sols[i]->sol[0]
            << sols.sols[i]->sol[0] << std::endl;
            n_sol_diff++;
         }
      }
   }
   std::cout << "Different Solution " << n_sol_diff << std::endl;

   // sols_phc.compare(sols);
}

void compare_path_with_phc
 ( Path& path_data )
{
   string file_name = "../Problems/test/cyclic10.path";
   ifstream path_file(file_name.c_str());
   if(path_file.is_open()==false)
   {
      std::cout << "Cannot open file" << file_name << std::endl;
   }
   int dim = 10;
   Path path_data_phc(dim,path_file);
   // path_data_phc.print_phc();
   path_data.print_phc();
   path_data_phc.compare(path_data);
}

bool path_multiple_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom, 
   Parameter path_parameter, PolySolSet& start_sol_set, PolySys& Target_Sys,
   int device_option, int n_path_test )
{
   bool cpu_test = 0;
   bool gpu_test = 0;
   if(device_option == 0)
   {
      cpu_test = 1;
      gpu_test = 1;
      std::cout << "CPU + GPU Testing..." << std::endl;
   }
   else if(device_option == 1)
   {
      cpu_test = 1;
      gpu_test = 0;
      std::cout << "CPU Testing..." << std::endl;
   }
   else if(device_option == 2)
   {
      cpu_test = 0;
      gpu_test = 1;
      std::cout << "GPU Testing..." << std::endl;
   }
   int n_path_success = 0;
   int n_path_fail = 0;
   int n_path = start_sol_set.n_sol;
   int start_idx = 0;
   // start_idx = 80;
   // n_path = min(start_idx+1,n_path);
   if(n_path_test > 0)
   {
      n_path = min(n_path_test,n_path);
   }
   path_report_set reports(n_path);
   PolySolSet sols(cpu_inst_hom.dim);
   // CT gamma = read_gamma(cpu_inst_hom.dim);
   // cpu_inst_hom.update_alpha(gamma);
   std::cout << "gamma = " << cpu_inst_hom.CPU_inst_hom_coef.alpha
             << std::endl;
   CT** cpu_sol = new CT*[n_path];
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      cpu_sol[path_idx] = new CT[cpu_inst_hom.dim];
   }
   bool* path_success_cpu = new bool[n_path];
   bool* path_success_gpu = NULL;
   CT** x_gpu = NULL;
   int dim = cpu_inst_hom.dim;
   double timeMS_Path_CPU;

   if(cpu_test==true)
   {
      struct timeval start, end;
      long seconds, useconds;

      gettimeofday(&start, NULL);
      for(int path_idx=start_idx; path_idx<n_path; path_idx++)
      {
         std::cout << "Sol " << path_idx << std::endl;
         CT* sol0 = start_sol_set.get_sol(path_idx);
         /*for(int var_idx=0; var_idx<cpu_inst_hom.dim; var_idx++)
           {
              std::cout << var_idx << " " << sol0[var_idx] << std::endl;
           }*/
         CT* sol_new = NULL;
         CT t(0.0,0.0);
         workspace_cpu.init_x_t_idx();
         workspace_cpu.update_x_t(sol0,t);
         workspace_cpu.path_idx = path_idx;
         path_success_cpu[path_idx] = path_test(workspace_cpu,cpu_inst_hom,
            path_parameter,sol0,sol_new,t,Target_Sys,1);

         /*string sol_info;
           if(path_success == 0)
           {
              sol_info = "failure";
              for(int sol_idx=0; sol_idx<cpu_inst_hom.dim; sol_idx++)
              {
                 if((abs(workspace_cpu.x_last[sol_idx].real) > 1E6)
                    or (abs(workspace_cpu.x_last[sol_idx].imag) > 1E6))
                 {
                    sol_info = "infinity";
                    break;
                 }
              }
           }
           else
           {
              sol_info = "success";
           }
           sols.add_sol(workspace_cpu.x_last,cpu_inst_hom.max_residual,
              cpu_inst_hom.max_delta_x, path_idx, sol_info);
         */
         for(int x_idx=0; x_idx<cpu_inst_hom.dim; x_idx++)
         {
            cpu_sol[path_idx][x_idx] = workspace_cpu.x_last[x_idx];
         }
         delete[] sol0;
         // cpu_inst_hom.path_data.print();
         // compare_path_with_phc(cpu_inst_hom.path_data);
         // reports.load_report(path_idx, cpu_inst_hom, workspace_cpu);
      }
      // std::cout << "n_path_success_cpu = " << n_path_success_cpu
      //           << std::endl;

      gettimeofday(&end, NULL);
      seconds  = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      timeMS_Path_CPU = seconds*1000 + useconds/1000.0;
      // cpu_inst_hom.timeSec_Path_CPU = timeMS_Path_CPU/1000;
   }
   if(gpu_test==true)
   {
      CT* sol0 = new CT[n_path*cpu_inst_hom.dim];
      CT* sol_tmp_mult = sol0;
      for(int sys_idx=0; sys_idx<n_path; sys_idx++)
      {
         CT* sol_tmp = start_sol_set.get_sol(sys_idx);
         for(int i=0; i<cpu_inst_hom.dim; i++)
         {
            *sol_tmp_mult++ = sol_tmp[i];
            // std::cout << i << " " << sol_tmp[i] << std::endl;
         }
      }
      CT* t = new CT[n_path];
      for(int sys_idx=0; sys_idx<n_path; sys_idx++)
      {
         t[sys_idx] = CT(0.0,0.0);
      }
      path_success_gpu = GPU_Path_mult(cpu_inst_hom,path_parameter,
         sol0,t,x_gpu,n_path);
      delete[] t;
   }
   if(cpu_test==true && gpu_test==true)
   {
      int n_path_success_cpu = 0;
      int n_path_success_gpu = 0;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         if(path_success_cpu[path_idx] != path_success_gpu[path_idx])
         {
            std::cout << path_idx << " ";
         }
         if(path_success_cpu[path_idx]==true)
         {
            n_path_success_cpu++;
         }
         if(path_success_gpu[path_idx]==true)
         {
            n_path_success_gpu++;
         }
      }
      std::cout << std::endl;
      std::cout << "n_path_success_cpu = " << n_path_success_cpu << std::endl;
      std::cout << "n_path_success_gpu = " << n_path_success_gpu << std::endl;

      T1 err = 0;
      int n_path_err = 0;
      for(int path_idx=start_idx; path_idx<n_path; path_idx++)
      {
         bool path_err = 0;
         for(int x_idx=0; x_idx<dim; x_idx++)
         {
            if((path_success_cpu[path_idx]==true))
            {
               CT err_tmp = x_gpu[path_idx][x_idx]-cpu_sol[path_idx][x_idx];
               if(abs(err_tmp.real) > err)
               {
                  err = err_tmp.real;
               }
               if(abs(err_tmp.imag) > err)
               {
                  err = err_tmp.imag;
               }
               if((abs(err_tmp.real) > 1E-5 || abs(err_tmp.imag) > 1E-5))
               {
                  path_err = true;
                  std::cout << path_idx << " " << x_idx << " " << err_tmp;
               }
            }
         }
         if(path_err)
         {
            n_path_err++;
            std::cout << path_idx << " ";
         }
      }
      std::cout << std::endl;
      std::cout << "n_path_err = " << n_path_err << std::endl;
      std::cout << "sol err: " << err << std::endl;
      std::cout << "x_gpu[0][0]   = " << x_gpu[0][0];
      std::cout << "cpu_sol[0][0] = " << cpu_sol[0][0];
   }
   if(cpu_test==true)
   {
      std::cout << "timeMS_Path_CPU = " << timeMS_Path_CPU << std::endl;
   }
   // reports.print_fail_path_reports();
   // reports.path_number_report();
   bool compare_with_phc = false;
   if(compare_with_phc)
   {
      compare_sol_with_phc(cpu_inst_hom.dim, sols);
   }
   // sols.print_short();
   return true;
}
