// The file lib2path.cpp contains the definitions of the functions
// with prototypes in lib2path_d.h, to link the Path library into PHCpack.

#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "complexH.h"
#include "newton_host.h"
#include "path_host.h"
#include "lib2path_d.h"

void var_name ( char* x_string, int x_string_len, string*& x_names, int& dim )
{
   dim = 1;
   for(int pos=0; pos<x_string_len; pos++)
      if(x_string[pos]==' ') dim++;
   x_names = new string[dim];

   int begin=0;
   string mystring = x_string;
   int var_idx = 0;
   for(int pos=0; pos<x_string_len; pos++)
   {
      if(x_string[pos]==' ')
      {
         x_names[var_idx] = mystring.substr(begin, pos-begin);
         begin = pos+1;
         var_idx++;
      }
   }
   x_names[var_idx] = mystring.substr(begin, x_string_len-begin);
}

void lib2path_read_standard_sys
 ( int verbose, PolySys<complexH<double>,double>& sys )
{
   int fail,nbsym;

   fail = syscon_number_of_symbols(&nbsym);
   if(verbose > 0)
   {
      std::cout << "the system is .." << std::endl;
      fail = syscon_write_standard_system();
      std::cout << "the number of symbols : " << nbsym << std::endl;
   }
   int s_dim = 80*nbsym;
   char *s = (char*) calloc(s_dim,sizeof(char));
   fail = syscon_string_of_symbols(&s_dim, s);
   string* x_names;
   int dim = 0;
   var_name(s, s_dim, x_names, dim);
   int i = 1;
   if(verbose > 0) std::cout << "dim = " << dim << std::endl;

   double c[2];
   int d[dim];

   int n_eq = 0;
   fail = syscon_number_of_standard_polynomials(&n_eq);

   sys.n_eq = n_eq;
   sys.dim  = dim;
   sys.eq_space = new PolyEq<complexH<double>,double>[n_eq];
   sys.pos_var = x_names;

   PolyEq<complexH<double>,double>* tmp_eq = sys.eq_space;

   for(int i=1; i<n_eq+1; i++)
   {
      int nt;
      fail = syscon_number_of_standard_terms(i,&nt);
      if(verbose > 0)
         std::cout << " #terms in polynomial " << i << " : " << nt << std::endl;
      tmp_eq->n_mon = nt;
      tmp_eq->dim = dim;
      for(int j=1; j<=nt; j++)
      {
         fail = syscon_retrieve_standard_term(i,j,dim,d,c);
         if(verbose > 0)
         {
            std::cout << c[0] << " " << c[1] << std::endl;
            for (int k=0; k<dim; k++) std::cout << " " << d[k];
            std::cout << std::endl;
         }
         bool constant_term = true;
         for(int k=0; k<dim; k++)
            if(d[k]!=0) constant_term = false;

         if(constant_term==true)
         {
            tmp_eq->n_mon--;
            tmp_eq->constant += complexH<double>(c[0],c[1]);
         }
         else
         {
            PolyMon<complexH<double>,double>* a
               = new PolyMon<complexH<double>,double>(dim,d,c);
            tmp_eq->mon.push_back(a);
         }
      }
      if(verbose > 0) tmp_eq->print(x_names);
      sys.eq.push_back(tmp_eq);
      tmp_eq++;
   }
   if(verbose > 0) sys.print();
}

void lib2path_read_standard_sols
 ( PolySys<complexH<double>,double>& start_sys,
   PolySolSet<complexH<double>,double>& sols )
{
   int fail, len;

   fail = solcon_number_of_standard_solutions(&len);
   // printf("Number of start solutions : %d\n",len);
   int dim = start_sys.dim;
   sols.dim = dim;
   double sol[2*dim+5];

   for(int sol_idx=1; sol_idx<len+1; sol_idx++)
   {
      int mm,k;

      solcon_retrieve_next_standard_solution(dim,&k,&mm,sol);
      //solcon_retrieve_solution(dim,sol_idx,&mm,sol);
      /*std::cout << sol[0] << " " << sol[1] << std::endl;
      for(int var_idx=0; var_idx<4; var_idx++)
         std::cout << sol[2+2*var_idx] << " " << sol[2+2*var_idx] << std::endl;
      std::cout << sol[2+2*dim] 
                << " " << sol[3+2*dim]
                << " " << sol[4+2*dim] << std::endl;*/
      PolySol<complexH<double>,double>* tmp_sol
         = new PolySol<complexH<double>,double>
                 (dim, sol[0], sol[1], sol+2,
                  sol[2+2*dim], sol[3+2*dim], sol[4+2*dim]);
      //tmp_sol->print_info(start_sys.pos_var);
      sols.add_sol(tmp_sol);
   }
   // std::cout << "sol finished" << std::endl;
}

void lib2path_write_standard_sols
 ( PolySolSet<complexH<double>,double>& sols )
{
   int fail = solcon_clear_standard_solutions();
   if(fail != 0)
      std::cout << "failed to clear the solutions" << std::endl;
   int dim = sols.dim;
   int nbsols = sols.n_sol;
   // std::cout << "number of solutions : " << nbsols << std::endl;
   complexH<double>* sol = new complexH<double>[dim];
   complexH<double> tval;
   double csol[2*dim+5];

   for(int sol_idx=0; sol_idx<nbsols; sol_idx++)
   {
      sols.get_solt(sol_idx,sol,&tval);

      csol[0] = tval.real;
      csol[1] = tval.imag;
      int idx = 2;
      for(int k=0; k<dim; k++)
      {
         csol[idx++] = sol[k].real;
         csol[idx++] = sol[k].imag;
         // std::cout << sol[k].real << "  " << sol[k].imag << std::endl;
         // std::cout << csol[idx-2] << "  " << csol[idx-1] << std::endl;
      }
      csol[2*dim+2] = 0.0;
      csol[2*dim+3] = 0.0;
      csol[2*dim+4] = 0.0;
      fail = solcon_append_standard_solution(dim,1,csol);
      if(fail != 0)
         std::cout << "failed to append the solution" << std::endl;
   }
}

void lib2path_read_standard_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<double>,double>& start_sys,
   PolySys<complexH<double>,double>& target_sys,
   PolySolSet<complexH<double>,double>& sols)
{
   int fail;

   std::cout << target_file << " " << strlen(target_file) << std::endl;
   std::cout << start_file << " " << strlen(start_file) << std::endl;
   fail = read_standard_target_system_from_file
            (strlen(target_file), target_file);
   lib2path_read_standard_sys(0,target_sys);

   fail = read_standard_start_system_from_file(strlen(start_file),start_file);
   lib2path_read_standard_sys(0,start_sys);
   lib2path_read_standard_sols(start_sys, sols);
}

using namespace std;

int standard_ade_newton_with_pars ( int verbose, Parameter pars )
{
   int fail;
   PolySys<complexH<double>,double> ps;
   PolySolSet<complexH<double>,double> sols;

   if(verbose > 0)
   {
      int dim, len;

      cout << endl;
      cout << "Newton's method ..." << endl;
      fail = syscon_number_of_standard_polynomials(&dim);
      cout << "number of polynomials : " << dim << endl;
      fail = solcon_number_of_standard_solutions(&len);
      cout << "number of solutions : " << len << endl;
   }

   lib2path_read_standard_sys(verbose,ps);
   lib2path_read_standard_sols(ps,sols);

   fail = standard_newton_with_pars(verbose,pars,ps,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   lib2path_write_standard_sols(sols);

   return 0;
}

int standard_ade_newton ( int verbose )
{
   Parameter pars(16);

   return standard_ade_newton_with_pars(verbose,pars);
}

extern "C" int standard_adenewton ( int verbose )
{
   return standard_ade_newton(verbose);
}

extern "C" int standard_adenewton_with_pars
 ( int verbose, int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine )
{
   Parameter pars(16);

   pars.max_step = max_step;
   pars.n_predictor = n_predictor;
   pars.step_increase = step_increase;
   pars.step_decrease = step_decrease;
   pars.max_delta_t = max_delta_t;
   pars.max_delta_t_end = max_delta_t_end;
   pars.min_delta_t = min_delta_t;
   pars.err_max_res = err_max_res;
   pars.err_max_delta_x = err_max_delta_x;
   pars.err_max_first_delta_x = err_max_first_delta_x;
   pars.max_it = max_it;
   pars.err_min_round_off = err_min_round_off;
   pars.max_it_refine = max_it_refine;
   pars.err_min_round_off_refine = err_min_round_off_refine;

   return standard_ade_newton_with_pars(verbose,pars);
}

int standard_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars )
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
      cout << "Tracking one path ..." << endl;
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

   fail = standard_onetrack_with_pars(verbose,regamma,imgamma,pars,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   lib2path_write_standard_sols(sols);

   return 0;
}

int standard_ade_onepath
 ( int verbose, double regamma, double imgamma )
{
   Parameter pars(16);

   return standard_ade_onepath_with_pars(verbose,regamma,imgamma,pars);
}

extern "C" int standard_adeonepath
 ( int verbose, double regamma, double imgamma )
{
   return standard_ade_onepath(verbose,regamma,imgamma);
}

extern "C" int standard_adeonepath_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine )
{
   Parameter pars(16);

   pars.max_step = max_step;
   pars.n_predictor = n_predictor;
   pars.step_increase = step_increase;
   pars.step_decrease = step_decrease;
   pars.max_delta_t = max_delta_t;
   pars.max_delta_t_end = max_delta_t_end;
   pars.min_delta_t = min_delta_t;
   pars.err_max_res = err_max_res;
   pars.err_max_delta_x = err_max_delta_x;
   pars.err_max_first_delta_x = err_max_first_delta_x;
   pars.max_it = max_it;
   pars.err_min_round_off = err_min_round_off;
   pars.max_it_refine = max_it_refine;
   pars.err_min_round_off_refine = err_min_round_off_refine;

   return standard_ade_onepath_with_pars(verbose,regamma,imgamma,pars);
}

int standard_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars )
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

   fail = standard_manytrack_with_pars(verbose,regamma,imgamma,pars,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   lib2path_write_standard_sols(sols);

   return 0;
}

int standard_ade_manypaths
 ( int verbose, double regamma, double imgamma )
{
   Parameter pars(16);

   return standard_ade_manypaths_with_pars(verbose,regamma,imgamma,pars);
}

extern "C" int standard_ademanypaths
 ( int verbose, double regamma, double imgamma )
{
   return standard_ade_manypaths(verbose,regamma,imgamma);
}

extern "C" int standard_ademanypaths_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine )
{
   Parameter pars(16);

   pars.max_step = max_step;
   pars.n_predictor = n_predictor;
   pars.step_increase = step_increase;
   pars.step_decrease = step_decrease;
   pars.max_delta_t = max_delta_t;
   pars.max_delta_t_end = max_delta_t_end;
   pars.min_delta_t = min_delta_t;
   pars.err_max_res = err_max_res;
   pars.err_max_delta_x = err_max_delta_x;
   pars.err_max_first_delta_x = err_max_first_delta_x;
   pars.max_it = max_it;
   pars.err_min_round_off = err_min_round_off;
   pars.max_it_refine = max_it_refine;
   pars.err_min_round_off_refine = err_min_round_off_refine;

   return standard_ade_manypaths_with_pars(verbose,regamma,imgamma,pars);
}

int standard_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySolSet<complexH<double>,double>& s )
{
   double teval,tmgs;
   bool success;
   complexH<double>* sol = s.get_sol(0);
   complexH<double> alpha;
   CPUInstHom<complexH<double>,double> cpu_inst_hom;
   Workspace< complexH<double> > workspace_cpu;
  
   if(verbose > 0)
   {
      cout << "The first solution on input :" << endl;
      for(int k=0; k<p.dim; k++)
      {
         cout << k << " :";
         cout << setw(24) << scientific << setprecision(16) << sol[k];
      }
   }

   alpha = complexH<double>(1,0);
   cpu_inst_hom.init(p,p,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   workspace_cpu.init_x_t_idx();
   workspace_cpu.update_x_t_value(sol,alpha);
   success = CPU_Newton(workspace_cpu,cpu_inst_hom,pars,teval,tmgs);
   if(verbose > 0)
   {
      cout.precision(16);
      cout << "The first solution after CPU_Newton :" << endl;
      for(int k=0; k<p.dim; k++)
         cout << k << " :" << setw(24) << workspace_cpu.x[k];
   }

   s.change_sol(0,workspace_cpu.x);

   return 0;
}

int standard_newton
 ( int verbose, PolySys<complexH<double>,double>& p,
   PolySolSet<complexH<double>,double>& s )
{
   Parameter pars(16);

   return standard_newton_with_pars(verbose,pars,p,s);
}

int standard_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
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
  
   if(verbose > 0)
   {
      cout << "The first solution on input :" << endl;
      for(int k=0; k<p.dim; k++)
      {
         cout << k << " :";
         cout << setw(24) << scientific << setprecision(16) << sol[k];
      }
   }

   alpha = complexH<double>(regamma,imgamma);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   t = complexH<double>(0.0,0);
   workspace_cpu.init_x_t_idx();
   workspace_cpu.update_x_t_value(sol,t);
   workspace_cpu.update_x_t_idx();
   tpred = 0.0; teval = 0.0; tmgs = 0.0;
   if(verbose > 0) cout << "calling path_tracker ..." << endl;
   success = path_tracker(workspace_cpu,cpu_inst_hom,pars,
                          tpred,teval,tmgs,0,verbose);
   if(verbose > 0) cout << "done with call to path_tracker." << endl;
   if(verbose > 0)
   {
      cout.precision(16);
      cout << "The first solution after CPU path tracker :" << endl;
      for(int k=0; k<p.dim; k++)
         cout << k << " :" << setw(24) << workspace_cpu.x_last[k];
   }
   s.change_solt(0,workspace_cpu.x_last,workspace_cpu.t_last);

   return 0;
}

int standard_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s )
{
   Parameter pars(16);

   return standard_onetrack_with_pars(verbose,regamma,imgamma,pars,p,q,s);
}

int standard_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
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
      // s.change_sol(path_idx,workspace_cpu.x_last);
      s.change_solt(path_idx,workspace_cpu.x_last,workspace_cpu.t_last);
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

int standard_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s )
{
   Parameter pars(16);

   return standard_manytrack_with_pars(verbose,regamma,imgamma,pars,p,q,s);
}

int set_path_parameter_value ( Parameter pars, int idx, double val )
{
   if((idx < 1) || (idx > 14))
      return -1;
   else
   {
      pars.set_value(idx,val);
      return 0;      
   }
}

int get_path_parameter_value ( Parameter pars, int idx, double* val )
{
   if((idx < 1) || (idx > 14))
      return -1;
   else
   {
      pars.get_value(idx,val);
      return 0;
   }
}

int get_default_path_parameters
 ( int precision, int* max_step, int* n_predictor,
   double* step_increase, double* step_decrease,
   double* max_delta_t, double* max_delta_t_end, double* min_delta_t,
   double* err_max_res, double* err_max_delta_x, double* err_max_first_delta_x,
   int* max_it, double* err_min_round_off,
   int* max_it_refine, double* err_min_round_off_refine )
{
   Parameter pars(precision);

   *max_step = pars.max_step;
   *n_predictor = pars.n_predictor;
   *step_increase = pars.step_increase;
   *step_decrease = pars.step_decrease;
   *max_delta_t = pars.max_delta_t;
   *max_delta_t_end = pars.max_delta_t_end;
   *min_delta_t = pars.min_delta_t;
   *err_max_res = pars.err_max_res;
   *err_max_delta_x = pars.err_max_delta_x;
   *err_max_first_delta_x = pars.err_max_first_delta_x;
   *max_it = pars.max_it;
   *err_min_round_off = pars.err_min_round_off;
   *max_it_refine = pars.max_it_refine;
   *err_min_round_off_refine = pars.err_min_round_off_refine;

   return 0;
}

int standard_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine )
{
   Parameter pars(16);

   pars.max_step = max_step;
   pars.n_predictor = n_predictor;
   pars.step_increase = step_increase;
   pars.step_decrease = step_decrease;
   pars.max_delta_t = max_delta_t;
   pars.max_delta_t_end = max_delta_t_end;
   pars.min_delta_t = min_delta_t;
   pars.err_max_res = err_max_res;
   pars.err_max_delta_x = err_max_delta_x;
   pars.err_max_first_delta_x = err_max_first_delta_x;
   pars.max_it = max_it;
   pars.err_min_round_off = err_min_round_off;
   pars.max_it_refine = max_it_refine;
   pars.err_min_round_off_refine = err_min_round_off_refine;

   return standard_ade_manypaths_with_pars(verbose,regamma,imgamma,pars);
}
