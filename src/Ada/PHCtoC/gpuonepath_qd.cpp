// C++ function called by Ada routine to accelerate the tracking of one path

#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "ada_test_qd.h"
#include "parameter.h"
#include "path_gpu.h"
#include "path_host.h"
#include "poly.h"

using namespace std;

int track ( int mode, int verbose, PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   The execution mode is 0 (CPU+GPU), 1 (CPU), or 2 (GPU).
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int gpuonepath_qd ( int mode, int verbose )
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   If mode = 0, then both CPU and GPU will execute,
 *   if mode = 1, then only CPU runs Newton's method,
 *   if mode = 2, then only GPU runs Newton's method.
 *   If verbose > 0, then additional output will be written. */
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
      cout << "Acceleration of tracking one path ..." << endl;
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

   fail = track(mode,verbose,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   ada_write_sols(sols);

   return 0;
}

int track ( int mode, int verbose, PolySys& p, PolySys& q, PolySolSet& s )
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

   alpha = CT(1.0,0);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   if(mode == 0 || mode == 1)
   {
      workspace_cpu.init_x_t_idx();
      workspace_cpu.update_x_t_value(sol,t);
      workspace_cpu.update_x_t_idx();
      tpred = 0.0; teval = 0.0; tmgs = 0.0;
      success = path_tracker(workspace_cpu,cpu_inst_hom,path_parameter,
                             tpred,teval,tmgs,0,verbose);
      if(verbose > 0)
      {
         cout.precision(64);
         cout << "The first solution after CPU path tracker :" << endl;
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(64) << workspace_cpu.x_last[k];
      }
   }

   if(mode == 0 || mode == 2)
   {
      t = CT(0.0,0);
      success = GPU_Path(cpu_inst_hom,path_parameter,sol,t,x_gpu,verbose);

      if(verbose > 0)
      {
         cout << "The first solution after GPU path tracker :" << endl;
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(72) << x_gpu[k];
      }
   }

   if(mode == 1)
      s.change_sol(0,workspace_cpu.x_last);
   else
      s.change_sol(0,x_gpu);

   return 0;
}
