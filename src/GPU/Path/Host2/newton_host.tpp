/* newton_host.tpp contains the templated code for newton_host.h */

double sum_norm ( double real, double imag )
{
   return abs(real) + abs(imag);
}

double max_norm ( complexH<double>* sol, int dim )
{
   double max_delta = sum_norm(sol[0].real,sol[0].imag);
   for(int k=1; k<dim; k++)
   {
      double tmp_delta = sum_norm(sol[k].real,sol[k].imag);
      if(tmp_delta>max_delta) max_delta = tmp_delta;
   }
   return max_delta;
}

double max_norm ( complexH<dd_real>* sol, int dim )
{
   double max_delta = sum_norm(sol[0].real.x[0],sol[0].imag.x[0]);
   for(int k=1; k<dim; k++)
   {
      double tmp_delta = sum_norm(sol[k].real.x[0],sol[k].imag.x[0]);
      if(tmp_delta>max_delta) max_delta = tmp_delta;
   }
   return max_delta;
}

double max_norm ( complexH<qd_real>* sol, int dim )
{
   double max_delta = sum_norm(sol[0].real.x[0],sol[0].imag.x[0]);
   for(int k=1; k<dim; k++)
   {
      double tmp_delta = sum_norm(sol[k].real.x[0],sol[k].imag.x[0]);
      if(tmp_delta>max_delta) max_delta = tmp_delta;
   }
   return max_delta;
}

template <class ComplexType, class RealType>
bool CPU_Newton 
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse, bool verbose )
{
   bool Debug = false; // Debug = true;
   // if(workspace_cpu.path_idx == 0) { // Debug = true; }
   bool Record = false; // Record = true;
   ComplexType* x = workspace_cpu.x;
   ComplexType t = *(workspace_cpu.t);
   ComplexType** V = (workspace_cpu.V);
   ComplexType** R = (workspace_cpu.R);
   ComplexType* sol = (workspace_cpu.sol);
   ComplexType* rhs = (workspace_cpu.rhs);
   int dim = cpu_inst_hom.dim;
   int n_eq = cpu_inst_hom.n_eq;
   // Parameters
   double err_round_off;
   if(t.real > 0.9)
      err_round_off = path_parameter.err_min_round_off_refine;
   else
      err_round_off = path_parameter.err_min_round_off;

   double max_x = max_norm(x,dim);
   bool success = 0;
   cpu_inst_hom.eval(workspace_cpu, x, t, reverse);
   cpu_inst_hom.n_eval_CPU++;
   double max_f_val = max_norm(V[dim],n_eq);
   double r_max_f_val = max_f_val/max_x;
   if(verbose) // if(Debug)
   {
      std::cout.precision(2);
      std::cout << scientific;
      std::cout << "          max_x : " << max_x << std::endl;
      std::cout << "   residual(a&r): " << max_f_val
                << " " << r_max_f_val << std::endl;
   }
   if(max_f_val < err_round_off || r_max_f_val < err_round_off)
   {
      success = 1;
      return success;
   }
   double last_max_f_val = max_f_val;
   double max_delta_x;
   for(int i=0; i<path_parameter.max_it; i++)
   {
      // if(Debug) std::cout << "Iteration " << i << std::endl;
      if(verbose) std::cout << "Iteration " << i << std::endl;

      CPU_mgs2qrls<ComplexType,RealType>(V,R,sol,n_eq,dim+1,rhs);
      cpu_inst_hom.n_mgs_CPU++;
      max_delta_x = max_norm(sol,dim);
      double r_max_delta_x = max_delta_x/max_x;
      if(verbose) // if(Debug)
         std::cout << " correction(a&r): " << max_delta_x  << " "
                   << r_max_delta_x;
      if(Record)
         cpu_inst_hom.path_data.add_iteration(max_delta_x, r_max_delta_x);

      // do the update anyway before break because of success
      for(int k=0; k<dim; k++) x[k] = x[k] - sol[k];

      if(max_delta_x < err_round_off || r_max_delta_x < err_round_off)
      {
         success = 1;
         break;
      }
      max_x = max_norm(x,dim);
      cpu_inst_hom.eval(workspace_cpu, x, t, reverse);
      cpu_inst_hom.n_eval_CPU++;
      max_f_val = max_norm(V[dim],n_eq);
      r_max_f_val = max_f_val/max_x;
      if(verbose) // if(Debug)
      {
         std::cout << std::endl << "          max_x : " << max_x << std::endl;
         std::cout << " residual(a&r): " << max_f_val << " "
                   << r_max_f_val << std::endl;
      }
      if(Record)
         cpu_inst_hom.path_data.update_iteration_res(max_f_val, r_max_f_val);

      if(max_f_val < err_round_off || r_max_f_val < err_round_off)
      {
         success = 1;
         break;
      }
      last_max_f_val = max_f_val;
   }
   if(verbose) // if(Debug)
   {
      std::cout << std::endl;
      std::cout << "Newton success = " << success << std::endl;
   }
   if(Record)
      cpu_inst_hom.path_data.update_step_correct_pt(x);

   return success;
}

template <class ComplexType, class RealType>
bool CPU_Newton_Refine
 ( Workspace<ComplexType>& workspace_cpu,
   CPUInstHom<ComplexType,RealType>& cpu_inst_hom,
   Parameter path_parameter, double& timeSec_Eval, double& timeSec_MGS,
   int reverse, bool verbose )
{
   // Parameters
   // eqs square to compare with normal square

   ComplexType* x = workspace_cpu.x;
   ComplexType t = *(workspace_cpu.t);
   ComplexType** V = (workspace_cpu.V);
   ComplexType** R = (workspace_cpu.R);
   ComplexType* sol = (workspace_cpu.sol);
   ComplexType* rhs = (workspace_cpu.rhs);
   int dim = cpu_inst_hom.dim;
   int n_eq = cpu_inst_hom.n_eq; // to be removed

   double last_max_delta_x = path_parameter.err_max_first_delta_x;
   double last_max_f_val   = path_parameter.err_max_res;
   bool success = 1;

   bool Debug = false; // Debug = true;

   for(int i=0; i<path_parameter.max_it_refine; i++)
   {
      // if(Debug) std::cout << "  Iteration " << i << std::endl;
      if(verbose) std::cout << "  Iteration " << i << std::endl;
      cpu_inst_hom.eval(workspace_cpu, x, t, reverse);
      cpu_inst_hom.n_eval_CPU++;
      double max_f_val = max_norm(V[dim],n_eq);
      if(verbose) // if(Debug)
         std::cout << "       max_f_value = " << max_f_val << std::endl;

      if(max_f_val > last_max_f_val || max_f_val != max_f_val)
      {
         success = 0;
         break;
      }
      if(max_f_val < path_parameter.err_min_round_off_refine)
      {
         last_max_delta_x = 0;
         break;
      }
      if(verbose) std::cout << "calling CPU_mgs2qrls ..." << std::endl;
      CPU_mgs2qrls<ComplexType,RealType>(V,R,sol,n_eq,dim+1,rhs);
      cpu_inst_hom.n_mgs_CPU++;
      double max_delta_x = max_norm(sol,dim);
      if(verbose) // if(Debug)
         std::cout << "       max_delta_x = " << max_delta_x << std::endl;

      // do update anyway, before breaking out because of success
      for(int k=0; k<dim; k++) x[k] = x[k] - sol[k];

      if(max_delta_x > last_max_delta_x || max_delta_x != max_delta_x)
      {
         success = 0;
         break;
      }
      if(max_delta_x < path_parameter.err_min_round_off)
      {
         last_max_delta_x = max_delta_x;
         break;
      }
      last_max_delta_x = max_delta_x;
      last_max_f_val = max_f_val;
   }
   if(success)
      if(last_max_delta_x > path_parameter.err_max_delta_x) success = 0;

   cpu_inst_hom.max_delta_x = sqrt(last_max_delta_x);

   return success;
}
