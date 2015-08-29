// C++ function called by Ada routine to accelerate the tracking of many paths.
// The prototypes of the functions are described in gpumanypaths_qd.h.

#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "ada_test_qd.h"
#include "parameter.h"
#include "path_host.h"
#include "path_gpu.h"
#include "poly.h"
#include "gpumanypaths_qd.h"

using namespace std;

extern "C" int gpumanypaths_qd
 ( int mode, int verbose, double regamma, double imgamma )
{
   int fail;
   PolySys ps;
   PolySys qs;
   PolySolSet sols;

   fail = copy_quaddobl_target_system_to_container();

   if(verbose > 0)
   {
      int dim, len;

      cout << endl;
      cout << "Acceleration of tracking many paths ..." << endl;
      cout << "gamma = " << setprecision(16)
           << regamma << " + i* " << imgamma << endl;
      cout << "Mode of execution : " << mode << endl;
      fail = syscon_number_of_quaddobl_polynomials(&dim);
      cout << "number of polynomials : " << dim << endl;
      fail = solcon_number_of_quaddobl_solutions(&len);
      cout << "number of solutions : " << len << endl;
   }

   ada_read_sys(verbose,ps);

   fail = copy_quaddobl_start_system_to_container();

   ada_read_sys(verbose,qs);
   ada_read_sols(qs,sols);

   fail = manytrack(mode,verbose,regamma,imgamma,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   ada_write_sols(sols);

   return 0;
}

int manytrack
 ( int mode, int verbose, double regamma, double imgamma,
   PolySys& p, PolySys& q, PolySolSet& s )
{
   double tpred,teval,tmgs;
   bool success;
   CT* sol = s.get_sol(0);
   CT alpha,t;
   CT *x_gpu;
   CPUInstHom cpu_inst_hom;
   Workspace workspace_cpu;
   Parameter path_parameter(N_PREDICTOR, MAX_STEP, MAX_IT, MAX_DELTA_T, \
      MAX_DELTA_T_END, MIN_DELTA_T, ERR_MAX_RES, ERR_MAX_DELTA_X, \
      ERR_MAX_FIRST_DELTA_X, ERR_MIN_ROUND_OFF, MAX_IT_REFINE, \
      ERR_MIN_ROUND_OFF_REFINE, STEP_INCREASE, STEP_DECREASE);
  
   if(verbose > 0)
   {
      cout << "The first solution on input :" << endl;
      for(int k=0; k<p.dim; k++)
      {
         cout << k << " :";
         cout << setw(72) << scientific << setprecision(64) << sol[k];
      }
   }
   int fail,n_path;
   fail = solcon_number_of_quaddobl_solutions(&n_path);

   alpha = CT(regamma,imgamma);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   if(mode == 0 || mode == 1)
   {
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         if(verbose > 0) cout << "tracking path " << path_idx << endl;
         CT* sol0 = s.get_sol(path_idx);
         CT* sol_new = NULL;
         t = CT(0.0,0);
         workspace_cpu.init_x_t_idx();
         workspace_cpu.update_x_t_value(sol0,t);
         workspace_cpu.update_x_t_idx();
         workspace_cpu.path_idx = path_idx;
         tpred = 0.0; teval = 0.0; tmgs = 0.0;
         if(verbose > 0) cout << "calling path_tracker ..." << endl;
         success = path_tracker(workspace_cpu,cpu_inst_hom,path_parameter,
                                tpred,teval,tmgs,0,verbose);
         if(verbose > 0) cout << "done with call to path_tracker." << endl;
         if(mode == 1) s.change_sol(path_idx,workspace_cpu.x_last);
         if(verbose > 0)
         {
            cout.precision(64);
            cout << "The solution " << path_idx << endl;
            for(int k=0; k<p.dim; k++)
               cout << k << " :" << setw(72) << workspace_cpu.x_last[k];
         }
      }
   }
   if(mode == 0 || mode == 2)
   {
      CT* sol0 = new CT[n_path*cpu_inst_hom.dim];
      CT* sol_tmp_mult = sol0;
      CT* tt = new CT[n_path];
      CT** x_gpu = NULL;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         CT* sol_tmp = s.get_sol(path_idx);
         for(int sol_idx=0; sol_idx<cpu_inst_hom.dim; sol_idx++)
         {
            *sol_tmp_mult++ = sol_tmp[sol_idx];
         }
         tt[path_idx] = CT(0.0,0.0);
      }
      success = GPU_Path_mult
         (cpu_inst_hom,path_parameter,sol0,tt,x_gpu,n_path,verbose);
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         s.change_sol(path_idx,x_gpu[path_idx]);
      }
      if(verbose > 0)
      {
         cout << "The first solution after GPU path tracker :" << endl;
         cout.precision(64);
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(72) << x_gpu[0][k];
      }
      delete[] tt;
   }
   return 0;
}
