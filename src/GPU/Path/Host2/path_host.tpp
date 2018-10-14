// The file path_host.tpp contains the definition of the function
// with prototype in the file path_host.h to track one single path.

template <class ComplexType, class RealType>
bool path_tracker
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter,
   double& timeSec_Predict, double& timeSec_Eval, double& timeSec_MGS,
   int reverse, int verbose )
{
   int n_point = 1;
   int n_step = 0;

   // Parameters
   cpu_inst_hom.n_eval_CPU = 0;
   cpu_inst_hom.n_mgs_CPU = 0;

   ComplexType delta_t = ComplexType(path_parameter.max_delta_t,0);

   // std::cout << "delta_t = " << delta_t << std::endl;

   ComplexType* tmp_t = workspace_cpu.t;
   ComplexType* tmp_t_last = workspace_cpu.t_last;

   // std::cout << "tmp_t = " << *tmp_t << std::endl;
   // std::cout << "tmp_t_last = " << *tmp_t_last << std::endl;

   int n_success = 0;
   string fail_reason;

   bool Debug = false;

   if(workspace_cpu.path_idx == 0)
   {
      //Debug = true;
   }

   bool Record = false;
   //Record = true;

   if(Record)
   {
      cpu_inst_hom.path_data.add_start_pt(workspace_cpu.x_last);
   }

   while(tmp_t_last->real < RealType(1))
   {
      if(Debug)
      {
         std::cout << "n_step = " << n_step  << ", n_point = "
                   << n_point  << std::endl;
      }
      if(delta_t.real + tmp_t_last->real < 1)
      {
         *tmp_t = *tmp_t_last + delta_t;
      }
      else
      {
         *tmp_t = ComplexType(1,0);
      }
      if(Debug)
      {
         std::cout << "delta_t = " << delta_t;
         std::cout << "      t = " << *tmp_t;
      }
      // clock_t begin_Predict = clock();
      int n_predictor = min(workspace_cpu.n_predictor, n_point);

      //predictor_newton<ComplexType>
      //   (workspace_cpu.x_array,workspace_cpu.t_array,
      //    workspace_cpu.x_t_idx,n_predictor,cpu_inst_hom.dim);

      predictor_divdif<ComplexType>
         (workspace_cpu.x_array,workspace_cpu.t_array,
          workspace_cpu.x_t_idx,n_predictor,cpu_inst_hom.dim,
          workspace_cpu.div_diff4pred, workspace_cpu.t_array4pred,
          workspace_cpu.t_diff4pred);

      // clock_t end_Predict = clock();
      // timeSec_Predict
      //   += (end_Predict - begin_Predict)
      //   / static_cast<double>( CLOCKS_PER_SEC );

      if(Debug)
      {
         std::cout << "Predict" << std::endl;
         int pr = 2 * sizeof(RealType);
         std::cout.precision(pr);
         for(int i=0; i<cpu_inst_hom.dim; i++)
         {
            std::cout << i << " " << workspace_cpu.x[i];
         }
      }
      if(Record)
      {
         cpu_inst_hom.path_data.add_step_empty();
         cpu_inst_hom.path_data.update_step_t(delta_t, *tmp_t);
         cpu_inst_hom.path_data.update_step_predict_pt(workspace_cpu.x);
      }
      bool newton_success = CPU_Newton<ComplexType,RealType>
         (workspace_cpu, cpu_inst_hom, path_parameter,
          timeSec_Eval, timeSec_MGS, reverse);

      if(newton_success == 1)
      {
         if(Debug)
         {
            std::cout << "Newton Success"<< std::endl;
            int pr = 2 * sizeof(RealType);
            std::cout.precision(pr);
            std::cout << "t = " << *tmp_t;
            for(int i=0; i<cpu_inst_hom.dim; i++)
            {
               std::cout << i << " " << workspace_cpu.x[i];
            }
         }
         if(tmp_t->real == 1)
         {
            CPU_Newton_Refine<ComplexType,RealType>
               (workspace_cpu,cpu_inst_hom,path_parameter,
                timeSec_Eval,timeSec_MGS,reverse);
         }
         n_point++;
         workspace_cpu.update_x_t_idx();
         tmp_t = workspace_cpu.t;
         tmp_t_last = workspace_cpu.t_last;
         n_success++;

         if(n_success > 1)
         {
            delta_t.real = delta_t.real*path_parameter.step_increase;
            // std::cout << "Increase delta_t = " << delta_t << std::endl;
         }
         RealType max_delta_t_real;
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
      }
      else
      {
         delta_t.real = delta_t.real*path_parameter.step_decrease;
         // std::cout << "Decrease delta_t = " << delta_t << std::endl;
         if(delta_t.real < path_parameter.min_delta_t)
         {
            fail_reason = "delta_t too small";
            break;
         }
         n_success = 0;
      }
      n_step++;
      if(n_step >= path_parameter.max_step)
      {
         fail_reason = "reach max step";
         break;
      }
      if(Debug)
      {
         std::cout << "---------------------"<< std::endl;
      }
      // std::cout << std::endl;
   }
   bool success = 0;
   if(verbose > 0)
      std::cout << "------------ Path Tracking Report -------------"
                << std::endl;
   if(tmp_t_last->real == 1)
   {
      success = 1;
      if(verbose > 0)
      {
         std::cout << "Success" << std::endl;
         std::cout << "t : " << *tmp_t_last;
      }
   }
   else
   {
      if(verbose > 0)
         std::cout << "Fail " << fail_reason << " t ="
                   << tmp_t_last->real<< std::endl;
   }
   if(Record)
   {
      cpu_inst_hom.path_data.add_end_pt(workspace_cpu.x_last);
      cpu_inst_hom.path_data.update_success(success);
   }
   cpu_inst_hom.success_CPU = success;
   cpu_inst_hom.t_CPU = *tmp_t_last;
   cpu_inst_hom.n_step_CPU  = n_step;
   cpu_inst_hom.n_point_CPU = n_point;

   if(verbose > 0)
   {
      std::cout << "n_point = " << n_point << std::endl;
      std::cout << "n_step = " << n_step << std::endl;
   }
   return success;
}
