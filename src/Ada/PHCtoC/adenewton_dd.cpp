// C++ function called by Ada routine to execute Newton's method
// using algorithmic differentiation for its evaluation
// in double double precision

#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "ada_test.h"
#include "parameter.h"
#include "poly.h"
#include "newton_host.h"
#include "adenewton_dd.h"

using namespace std;

extern "C" int adenewton_dd ( int verbose )
{
   int fail;
   PolySys ps;
   PolySolSet sols;

   if(verbose > 0)
   {
      int dim, len;

      cout << endl;
      cout << "Newton's method ..." << endl;
      fail = syscon_number_of_dobldobl_polynomials(&dim);
      cout << "number of polynomials : " << dim << endl;
      fail = solcon_number_of_dobldobl_solutions(&len);
      cout << "number of solutions : " << len << endl;
   }

   ada_read_sys(verbose,ps);
   ada_read_sols(ps,sols);

   fail = newton(verbose,ps,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   ada_write_sols(sols);

   return 0;
}

int newton ( int verbose, PolySys& p, PolySolSet& s )
{
   double teval,tmgs;
   bool success;
   CT* sol = s.get_sol(0);
   CT alpha;
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
         cout << setw(40) << scientific << setprecision(32) << sol[k];
      }
   }

   alpha = CT(1,0);
   cpu_inst_hom.init(p,p,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   workspace_cpu.init_x_t_idx();
   workspace_cpu.update_x_t_value(sol,alpha);
   success = CPU_Newton
                (workspace_cpu,cpu_inst_hom,path_parameter,teval,tmgs);
   if(verbose > 0)
   {
      cout.precision(32);
      cout << "The first solution after CPU_Newton :" << endl;
      for(int k=0; k<p.dim; k++)
         cout << k << " :" << setw(40) << workspace_cpu.x[k];
   }

   s.change_sol(0,workspace_cpu.x);

   return 0;
}
