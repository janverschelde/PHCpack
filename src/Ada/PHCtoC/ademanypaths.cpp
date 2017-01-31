// C++ function called by Ada routine to apply algorihtmic differentiation
// for the evaluation to the tracking of many paths

#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "lib2path.h"
// #include "parameter.h"
#include "path_host.h"
#include "ademanypaths.h"

using namespace std;

// Parameters
#define N_PREDICTOR           4

#define MAX_STEP              400
#define MAX_DELTA_T           1E-1
#define MAX_DELTA_T_END       1E-2
#define MIN_DELTA_T           1E-7

#define MAX_IT                3
#define ERR_MIN_ROUND_OFF     1E-9

#define MAX_IT_REFINE                   5
#define ERR_MIN_ROUND_OFF_REFINE    1E-11

#define ERR_MAX_RES           1E-2
#define ERR_MAX_DELTA_X       1E-1
#define ERR_MAX_FIRST_DELTA_X 1E-2

#define STEP_INCREASE   1.25
#define STEP_DECREASE   0.7

extern "C" int standard_ademanypaths
 ( int verbose, double regamma, double imgamma )
{
   int fail;
   PolySys<complexH<double>,double> ps;
   PolySys<complexH<double>,double> qs;
   PolySolSet<complexH<double>,double> sols;

   fail = copy_target_system_to_container();

   if(verbose > 0)
   {
      int dim, len;

      cout << endl;
      cout << "Tracking many paths ..." << endl;
      cout << "gamma = " << setprecision(16)
           << regamma << " + i* " << imgamma << endl;
      fail = syscon_number_of_standard_polynomials(&dim);
      cout << "number of polynomials : " << dim << endl;
      fail = solcon_number_of_standard_solutions(&len);
      cout << "number of solutions : " << len << endl;
   }

   lib2path_read_standard_sys(verbose,ps);

   fail = copy_start_system_to_container();

   lib2path_read_standard_sys(verbose,qs);
   lib2path_read_standard_sols(qs,sols);

   fail = standard_manytrack(verbose,regamma,imgamma,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   lib2path_write_standard_sols(sols);

   return 0;
}

int standard_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s )
{
   double tpred,teval,tmgs;
   bool success;
   complexH<double>* sol = s.get_sol(0);
   complexH<double> alpha,t;
   CPUInstHom<complexH<double>,double> cpu_inst_hom;
   Workspace< complexH<double> > workspace_cpu;
   Parameter pars(N_PREDICTOR, MAX_STEP, MAX_IT, MAX_DELTA_T, \
      MAX_DELTA_T_END, MIN_DELTA_T, ERR_MAX_RES, ERR_MAX_DELTA_X, \
      ERR_MAX_FIRST_DELTA_X, ERR_MIN_ROUND_OFF, MAX_IT_REFINE, \
      ERR_MIN_ROUND_OFF_REFINE, STEP_INCREASE, STEP_DECREASE);
  
   if(verbose > 0)
   {
      cout << "The first solution on input :" << endl;
      for(int k=0; k<p.dim; k++)
      {
         cout << k << " :";
         cout << setw(24) << scientific << setprecision(16) << sol[k];
      }
   }
   int fail,n_path;
   fail = solcon_number_of_standard_solutions(&n_path);

   alpha = complexH<double>(regamma,imgamma);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      if(verbose > 0) cout << "tracking path " << path_idx << endl;
      complexH<double>* sol0 = s.get_sol(path_idx);
      complexH<double>* sol_new = NULL;
      t = complexH<double>(0.0,0);
      workspace_cpu.init_x_t_idx();
      workspace_cpu.update_x_t_value(sol0,t);
      workspace_cpu.update_x_t_idx();
      workspace_cpu.path_idx = path_idx;
      tpred = 0.0; teval = 0.0; tmgs = 0.0;
      if(verbose > 0) cout << "calling path_tracker ..." << endl;
      success = path_tracker(workspace_cpu,cpu_inst_hom,pars,
                             tpred,teval,tmgs,0,verbose);
      if(verbose > 0) cout << "done with call to path_tracker." << endl;
      s.change_sol(path_idx,workspace_cpu.x_last);
      if(verbose > 0)
      {
         cout.precision(16);
         cout << "The solution " << path_idx << endl;
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(24) << workspace_cpu.x_last[k];
      }
   }
   return 0;
}
