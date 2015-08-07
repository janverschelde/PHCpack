// C++ function called by Ada routine to accelerate Newton's method

#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "ada_test.h"
#include "parameter.h"
#include "path_gpu.h"
#include "poly.h"
#include "newton_host.h"

using namespace std;

int newton ( PolySys& p, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s. */

extern "C" int gpunewton_d ( void )
/*
 * DESCRIPTION :
 *   A C++ function to accelerate Newton's method,
 *   encapsulated as a C function for to be called from Ada. */
{
   int fail,dim,len;
   PolySys ps;
   PolySolSet sols;

   cout << endl;
   cout << "Acceleration of Newton's method ..." << endl;

   fail = syscon_number_of_polynomials(&dim);
   cout << "number of polynomials : " << dim << endl;
   fail = solcon_number_of_solutions(&len);
   cout << "number of solutions : " << len << endl;

   ada_read_sys(ps);
   ada_read_sols(ps,sols);

   fail = newton(ps,sols);

   cout << "writing the solutions to the container ..." << endl;
   ada_write_sols(sols);

   return 0;
}

int newton ( PolySys& p, PolySolSet& s )
{
   double teval,tmgs;
   bool success;
   CT* sol = s.get_sol(0);
   CT alpha;
   CT *x_gpu;
   CPUInstHom cpu_inst_hom;
   Workspace workspace_cpu;
   Parameter path_parameter(N_PREDICTOR, MAX_STEP, MAX_IT, MAX_DELTA_T, \
      MAX_DELTA_T_END, MIN_DELTA_T, ERR_MAX_RES, ERR_MAX_DELTA_X, \
      ERR_MAX_FIRST_DELTA_X, ERR_MIN_ROUND_OFF, MAX_IT_REFINE, \
      ERR_MIN_ROUND_OFF_REFINE, STEP_INCREASE, STEP_DECREASE);

   cout << "The first solution on input :" << endl;
   for(int k=0; k<p.dim; k++)
   {
      cout << k << " :";
      cout << setw(24) << scientific << setprecision(16) << sol[k];
   }

   alpha = CT(1,0);
   cpu_inst_hom.init(p,p,p.dim,p.n_eq,1,alpha);
   cpu_inst_hom.init_workspace(workspace_cpu);
   workspace_cpu.init_x_t_idx();
   workspace_cpu.update_x_t_value(sol,alpha);

   success = CPU_Newton(workspace_cpu,cpu_inst_hom,path_parameter,teval,tmgs);
   cout.precision(16);
   cout << "The first solution after CPU_Newton :" << endl;
   for(int k=0; k<p.dim; k++)
      cout << k << " :" << setw(24) << workspace_cpu.x[k];

   alpha = CT(1,0);
   success = GPU_Newton(cpu_inst_hom,path_parameter,sol,alpha,x_gpu);
   cout << "The first solution after GPU_Newton :" << endl;
   for(int k=0; k<p.dim; k++)
      cout << k << " :" << setw(24) << x_gpu[k];

   s.change_sol(0,x_gpu);

   return 0;
}
