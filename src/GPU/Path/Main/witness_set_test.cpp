/* The file "witness_set_test.cpp" was written by Xiangcheng Yu
 * on 8 February 2015. */

#include "witness_set_test.h"

int witness_set_test
 ( Workspace& workspace_cpu, CPUInstHom& cpu_inst_hom,
   Parameter path_parameter, CT* x_cpu, CT t, int device_option )
{
   std::cout << "--------- Witness Set Test ----------" << std::endl;

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
   CT* x_new = NULL;

   PolySolSet witness_set_start(cpu_inst_hom.dim);
   PolySolSet witness_set_target(cpu_inst_hom.dim);

   witness_set_start.add_diff_sol(x_cpu);
   std::cout << "Witness Set1" << std::endl;
   witness_set_start.print_short();

   double timeSec_Predict = 0;
   double timeSec_Eval = 0;
   double timeSec_MGS = 0;
   double timeSec = 0;

   for(int i=0; i<10; i++)
   {
      clock_t begin = clock();

      bool success = 0;
      int alpha_round = 0;
      while(success==0)
      {
         cpu_inst_hom.path_data.clear();
         t = CT(0.0,0.0);
         workspace_cpu.init_x_t_idx();
         workspace_cpu.update_x_t(x_cpu, t);
         cpu_inst_hom.update_alpha();
         if(cpu_test == 1)
         {
            success = path_tracker(workspace_cpu,cpu_inst_hom,path_parameter,
                                   timeSec_Predict, timeSec_Eval, timeSec_MGS);
         }
         if(gpu_test == 1)
         {
            success = GPU_Path(cpu_inst_hom,path_parameter,x_cpu,t,x_new);
         }
         alpha_round++;
         if(success == 0 && alpha_round > 5)
         {
            break;
         }
      }
      if(success == 0 && alpha_round > 5)
      {
         std::cout << "Path Fail! Start -> Target" << std::endl \
                   << "alpha round = " << alpha_round << std::endl \
                   << "path  round = " << i << std::endl;
         break;
      }
      cpu_inst_hom.path_data.print_phc();
      clock_t end = clock();
      timeSec += (end - begin) / static_cast<double>( CLOCKS_PER_SEC );

      if(cpu_test == 1)
      {
         workspace_cpu.copy_x_last(x_cpu);
      }
      if(gpu_test == 1)
      {
         delete[] x_cpu;
         x_cpu = x_new;
      }
      witness_set_target.add_diff_sol(x_cpu);

      std::cout << "Witness Set Target" << std::endl;
      witness_set_target.print_short();
      if(witness_set_target.n_sol == 4 and witness_set_start.n_sol == 4)
      {
         break;
      }
      // Path Tracking Reverse Test
      success = 0;
      alpha_round = 0;
      while(success==0)
      {
         cpu_inst_hom.path_data.clear();
         t = CT(0.0,0.0);
         workspace_cpu.init_x_t_idx();
         workspace_cpu.update_x_t(x_cpu, t);
         cpu_inst_hom.update_alpha();
         // success = path_tracker(workspace_cpu,cpu_inst_hom,path_parameter,
         //                        timeSec_Predict,timeSec_Eval,timeSec_MGS,1);
         if(cpu_test == 1)
         {
            success = path_tracker(workspace_cpu,cpu_inst_hom,path_parameter,
                                   timeSec_Predict,timeSec_Eval,timeSec_MGS,1);
         }
         if(gpu_test == 1)
         {
            success = GPU_Path(cpu_inst_hom,path_parameter,x_cpu,t,x_new,1);
         }
         alpha_round++;
         if(success == 0 && alpha_round > 5)
         {
            break;
         }
      }
      cpu_inst_hom.path_data.print_phc();

      if(success == 0 && alpha_round > 5)
      {
         std::cout << "Path Fail! Target -> Start" << std::endl
                   << "alpha round = " << alpha_round << std::endl
                   << "path  round = " << i << std::endl;
         break;
      }
      if(cpu_test == 1)
      {
         workspace_cpu.copy_x_last(x_cpu);
      }
      if(gpu_test == 1)
      {
         delete[] x_cpu;
         x_cpu = x_new;
      }
      witness_set_start.add_diff_sol(x_cpu);

      std::cout << "Witness Set Start" << std::endl;
      witness_set_start.print_short();
      if(witness_set_target.n_sol == 4 and witness_set_start.n_sol == 4)
      {
         break;
      }
   }
   std::cout << "Witness Set Start" << std::endl;
   witness_set_start.print_short();
   std::cout << "Witness Set Target" << std::endl;
   witness_set_target.print_short();

   /*for(int i=0; i<cpu_inst_hom.dim; i++)
     {
        std::cout << i << " " << x_cpu[i];
     }*/

    /*CT* x_gpu = GPU_Path(cpu_inst_hom, sol0, t, n_predictor, max_it,
                           err_max_delta_x, max_step);
      cout << "Path CPU Test      Time: "<< timeSec << endl;
      std::cout << "--------- Path Tracking Error CPU vs GPU ----------"
                << std::endl;
      T1 err = err_check_workspace(x_cpu, x_gpu, cpu_inst_hom.dim);
      std::cout << " x_cpu[0] = " << x_cpu[0];
      std::cout << " x_gpu[0] = " << x_gpu[0];
      free(x_gpu);
     */

   // T1 err = err_check_workspace(sol0, x_cpu, cpu_inst_hom.dim);

   cout << "Path CPU Predict   Time: "<< timeSec_Predict << endl;
   cout << "Path CPU Eval      Time: "<< timeSec_Eval << endl;
   cout << "Path CPU MGS       Time: "<< timeSec_MGS << endl;

   return witness_set_start.n_sol;
}
