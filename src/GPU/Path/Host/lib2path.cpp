#include <iostream>
#include <iomanip>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "complexH.h"
#include "path_host.h"
#include "lib2path.h"

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

void lib2path_read_dobldobl_sys
 ( int verbose, PolySys<complexH<dd_real>,dd_real>& sys )
{
   int fail,nbsym;

   fail = syscon_number_of_symbols(&nbsym);

   if(verbose > 0)
   {
      std::cout << "the system is .." << std::endl;
      fail = syscon_write_dobldobl_system();
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

   double c[4]; /* two consecutive double doubles are real and imag parts */
   int d[dim];

   int n_eq = 0;
   fail = syscon_number_of_dobldobl_polynomials(&n_eq);

   sys.n_eq = n_eq;
   sys.dim  = dim;
   sys.eq_space = new PolyEq<complexH<dd_real>,dd_real>[n_eq];
   sys.pos_var = x_names;

   PolyEq<complexH<dd_real>,dd_real>* tmp_eq = sys.eq_space;

   for(int i=1; i<n_eq+1; i++)
   {
      int nt;
      fail = syscon_number_of_dobldobl_terms(i,&nt);
      if(verbose > 0)
         std::cout << " #terms in polynomial " << i << " : " << nt << std::endl;
      tmp_eq->n_mon = nt;
      tmp_eq->dim = dim;
      for(int j=1; j<=nt; j++)
      {
         fail = syscon_retrieve_dobldobl_term(i,j,dim,d,c);
         if(verbose > 0)
         {
            std::cout << c[0] << " " << c[2] << std::endl;
            for (int k=0; k<dim; k++) std::cout << " " << d[k];
            std::cout << std::endl;
         }
         bool constant_term = true;
         for(int k=0; k<dim; k++)
            if(d[k]!=0) constant_term = false;

         if(constant_term==true)
         {
            tmp_eq->n_mon--;
            tmp_eq->constant += complexH<dd_real>(c[0],c[2]);
         }
         else
         {
            dd_real cffs[2];
            dd_real realpcff;
            dd_real imagpcff;
            realpcff.x[0] = c[0];
            realpcff.x[1] = c[1];
            imagpcff.x[0] = c[2];
            imagpcff.x[1] = c[3];
            cffs[0] = realpcff;
            cffs[1] = imagpcff;
            PolyMon<complexH<dd_real>,dd_real>* a
               = new PolyMon<complexH<dd_real>,dd_real>(dim,d,cffs);
            tmp_eq->mon.push_back(a);
         }
      }
      if(verbose > 0) tmp_eq->print(x_names);
      sys.eq.push_back(tmp_eq);
      tmp_eq++;
   }
   if(verbose > 0)
   {
      sys.print();
      std::cout << "End" << std::endl;
   }
}

void lib2path_read_quaddobl_sys
 ( int verbose, PolySys<complexH<qd_real>,qd_real>& sys )
{
   int fail,nbsym;
   
   fail = syscon_number_of_symbols(&nbsym);
   if(verbose > 0)
   {
      std::cout << "the system is .." << std::endl;
      fail = syscon_write_quaddobl_system();
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

   double c[8]; /* two consecutive quad doubles are real and imag parts */
   int d[dim];

   int n_eq = 0;
   fail = syscon_number_of_quaddobl_polynomials(&n_eq);

   sys.n_eq = n_eq;
   sys.dim  = dim;
   sys.eq_space = new PolyEq<complexH<qd_real>,qd_real>[n_eq];
   sys.pos_var = x_names;

   PolyEq<complexH<qd_real>,qd_real>* tmp_eq = sys.eq_space;

   for(int i=1; i<n_eq+1; i++)
   {
      int nt;
      fail = syscon_number_of_quaddobl_terms(i,&nt);
      if(verbose > 0)
         std::cout << " #terms in polynomial " << i << " : " << nt << std::endl;
      tmp_eq->n_mon = nt;
      tmp_eq->dim = dim;
      for(int j=1; j<=nt; j++)
      {
         fail = syscon_retrieve_quaddobl_term(i,j,dim,d,c);
         if(verbose > 0)
         {
            std::cout << c[0] << " " << c[2] << std::endl;
            for (int k=0; k<dim; k++) std::cout << " " << d[k];
            std::cout << std::endl;
         }
         bool constant_term = true;
         for(int k=0; k<dim; k++)
            if(d[k]!=0) constant_term = false;

         if(constant_term==true)
         {
            tmp_eq->n_mon--;
            tmp_eq->constant += complexH<qd_real>(c[0],c[4]);
            //std::cout << "constant " << c[0] << " " << c[1] << std::endl;
         }
         else
         {
            qd_real cffs[4];
            qd_real realpcff;
            qd_real imagpcff;
            realpcff.x[0] = c[0];
            realpcff.x[1] = c[1];
            realpcff.x[2] = c[2];
            realpcff.x[3] = c[3];
            imagpcff.x[0] = c[4];
            imagpcff.x[1] = c[5];
            imagpcff.x[2] = c[6];
            imagpcff.x[3] = c[7];
            cffs[0] = realpcff;
            cffs[1] = imagpcff;
            PolyMon<complexH<qd_real>,qd_real>* a
               = new PolyMon<complexH<qd_real>,qd_real>(dim,d,cffs);
            tmp_eq->mon.push_back(a);
         }
      }
      if(verbose > 0) tmp_eq->print(x_names);
      sys.eq.push_back(tmp_eq);
      tmp_eq++;
   }
   if(verbose > 0)
   {
      sys.print();
      std::cout << "End" << std::endl;
   }
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

void lib2path_read_dobldobl_sols
 ( PolySys<complexH<dd_real>,dd_real>& start_sys,
   PolySolSet<complexH<dd_real>,dd_real>& sols )
{
   int fail, len;

   fail = solcon_number_of_dobldobl_solutions(&len);
   // printf("Number of start solutions : %d\n",len);
   int dim = start_sys.dim;
   sols.dim = dim;
   double sol[4*dim+10];

   for(int sol_idx=1; sol_idx<len+1; sol_idx++)
   {
      int mm,k,pos;
      // dd_real dd_sol[2*dim];
      dd_real* dd_sol;
      dd_real t_real,t_imag,err,rco,res;

      dd_sol = (dd_real*)calloc(2*dim, sizeof(dd_real));

      solcon_retrieve_next_dobldobl_solution(dim,&k,&mm,sol);
      //solcon_retrieve_solution(dim,sol_idx,&mm,sol);
      /*std::cout << sol[0] << " " << sol[1] << std::endl;
      for(int var_idx=0; var_idx<4; var_idx++)
         std::cout << sol[2+2*var_idx] << " " << sol[2+2*var_idx] << std::endl;
      std::cout << sol[2+2*dim] 
                << " " << sol[3+2*dim]
                << " " << sol[4+2*dim] << std::endl;*/
      t_real.x[0] = sol[0];
      t_real.x[1] = sol[1];
      t_imag.x[0] = sol[2];
      t_imag.x[1] = sol[3];
      pos = 4;
      for(int dd_sol_idx=0; dd_sol_idx<2*dim; dd_sol_idx++)
      {
         dd_real sol_real,sol_imag;
         sol_real.x[0] = sol[pos++];
         sol_real.x[1] = sol[pos++];
         sol_imag.x[0] = sol[pos++];
         sol_imag.x[1] = sol[pos++];
         dd_sol[dd_sol_idx++] = sol_real;
         dd_sol[dd_sol_idx] = sol_imag;
      }
      err.x[0] = sol[4*dim+4];
      err.x[1] = sol[4*dim+5];
      rco.x[0] = sol[4*dim+6];
      rco.x[1] = sol[4*dim+7];
      res.x[0] = sol[4*dim+8];
      res.x[1] = sol[4*dim+9];
      PolySol<complexH<dd_real>,dd_real>* tmp_sol
         = new PolySol<complexH<dd_real>,dd_real>
                  (dim,t_real,t_imag,dd_sol,err,rco,res);
      //tmp_sol->print_info(start_sys.pos_var);
      sols.add_sol(tmp_sol);
   }

   // std::cout << "sol finished" << std::endl;
}

void lib2path_read_quaddobl_sols
 ( PolySys<complexH<qd_real>,qd_real>& start_sys,
   PolySolSet<complexH<qd_real>,qd_real>& sols )
{
   int fail, len;

   fail = solcon_number_of_quaddobl_solutions(&len);
   // printf("Number of start solutions : %d\n",len);
   int dim=start_sys.dim;
   sols.dim = dim;
   double sol[8*dim+20];

   for(int sol_idx=1; sol_idx<len+1; sol_idx++)
   {
      int mm,k,pos;
      // qd_real qd_sol[4*dim];
      qd_real* qd_sol;
      qd_real t_real,t_imag,err,rco,res;

      qd_sol = (qd_real*)calloc(4*dim,sizeof(qd_real));

      solcon_retrieve_next_quaddobl_solution(dim,&k,&mm,sol);
      //solcon_retrieve_solution(dim,sol_idx,&mm,sol);
      /*std::cout << sol[0] << " " << sol[1] << std::endl;
      for(int var_idx=0; var_idx<4; var_idx++)
         std::cout << sol[2+2*var_idx] << " " << sol[2+2*var_idx] << std::endl;
      std::cout << sol[2+2*dim] 
                << " " << sol[3+2*dim]
                << " " << sol[4+2*dim] << std::endl;*/
      t_real.x[0] = sol[0];
      t_real.x[1] = sol[1];
      t_real.x[2] = sol[2];
      t_real.x[3] = sol[3];
      t_imag.x[0] = sol[4];
      t_imag.x[1] = sol[5];
      t_imag.x[2] = sol[6];
      t_imag.x[3] = sol[7];
      pos = 8;
      for(int qd_sol_idx=0; qd_sol_idx<4*dim; qd_sol_idx++)
      {
         qd_real sol_real,sol_imag;
         sol_real.x[0] = sol[pos++];
         sol_real.x[1] = sol[pos++];
         sol_real.x[2] = sol[pos++];
         sol_real.x[3] = sol[pos++];
         sol_imag.x[0] = sol[pos++];
         sol_imag.x[1] = sol[pos++];
         sol_imag.x[2] = sol[pos++];
         sol_imag.x[3] = sol[pos++];
         qd_sol[qd_sol_idx++] = sol_real;
         qd_sol[qd_sol_idx] = sol_imag;
      }
      err.x[0] = sol[8*dim+8];
      err.x[1] = sol[8*dim+9];
      err.x[2] = sol[8*dim+10];
      err.x[3] = sol[8*dim+11];
      rco.x[0] = sol[8*dim+12];
      rco.x[1] = sol[8*dim+13];
      rco.x[2] = sol[8*dim+14];
      rco.x[3] = sol[8*dim+15];
      res.x[0] = sol[8*dim+16];
      res.x[1] = sol[8*dim+17];
      res.x[2] = sol[8*dim+18];
      res.x[3] = sol[8*dim+19];
      PolySol<complexH<qd_real>,qd_real>* tmp_sol
         = new PolySol<complexH<qd_real>,qd_real>
                  (dim,t_real,t_imag,qd_sol,err,rco,res);
      // tmp_sol->print_info(start_sys.pos_var);
      sols.add_sol(tmp_sol);
   }
   // sols.print_info(start_sys.pos_var);
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

   for(int sol_idx=0; sol_idx<nbsols; sol_idx++)
   {
      complexH<double>* sol = sols.get_sol(sol_idx);
      double csol[2*dim+5];
      csol[0] = 0.0;
      csol[1] = 0.0;
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

void lib2path_write_dobldobl_sols
 ( PolySolSet<complexH<dd_real>,dd_real>& sols )
{
   int fail = solcon_clear_dobldobl_solutions();
   if(fail != 0)
      std::cout << "failed to clear the solutions" << std::endl;
   int dim = sols.dim;
   int nbsols = sols.n_sol;

   for(int sol_idx=0; sol_idx<nbsols; sol_idx++)
   {
      complexH<dd_real>* sol = sols.get_sol(sol_idx);
      double csol[4*dim+10];
      csol[0] = 0.0; csol[1] = 0.0;
      csol[2] = 0.0; csol[3] = 0.0;
      int idx = 4;
      for(int k=0; k<dim; k++)
      {
         csol[idx++] = sol[k].real.x[0];
         csol[idx++] = sol[k].real.x[1];
         csol[idx++] = sol[k].imag.x[0];
         csol[idx++] = sol[k].imag.x[1];
      }
      csol[4*dim+4] = 0.0; csol[4*dim+5] = 0.0;
      csol[4*dim+6] = 0.0; csol[4*dim+7] = 0.0;
      csol[4*dim+8] = 0.0; csol[4*dim+9] = 0.0;
      fail = solcon_append_dobldobl_solution(dim,1,csol);
      if(fail != 0)
         std::cout << "failed to append the solution" << std::endl;
   }
}

void lib2path_write_quaddobl_sols
 ( PolySolSet<complexH<qd_real>,qd_real>& sols )
{
   int fail = solcon_clear_quaddobl_solutions();
   if(fail != 0)
      std::cout << "failed to clear the solutions" << std::endl;
   int dim = sols.dim;
   int nbsols = sols.n_sol;

   for(int sol_idx=0; sol_idx<nbsols; sol_idx++)
   {
      complexH<qd_real>* sol = sols.get_sol(sol_idx);
      double csol[8*dim+20];
      csol[0] = 0.0; csol[1] = 0.0; csol[2] = 0.0; csol[3] = 0.0;
      csol[4] = 0.0; csol[5] = 0.0; csol[6] = 0.0; csol[7] = 0.0;
      int idx = 8;
      for(int k=0; k<dim; k++)
      {
         csol[idx++] = sol[k].real.x[0]; csol[idx++] = sol[k].real.x[1];
         csol[idx++] = sol[k].real.x[2]; csol[idx++] = sol[k].real.x[3];
         csol[idx++] = sol[k].imag.x[0]; csol[idx++] = sol[k].imag.x[1];
         csol[idx++] = sol[k].imag.x[2]; csol[idx++] = sol[k].imag.x[3];
      }
      csol[8*dim+8] = 0.0; csol[8*dim+9] = 0.0;
      csol[8*dim+10] = 0.0; csol[8*dim+11] = 0.0;
      csol[8*dim+12] = 0.0; csol[8*dim+13] = 0.0;
      csol[8*dim+14] = 0.0; csol[8*dim+15] = 0.0;
      csol[8*dim+16] = 0.0; csol[8*dim+17] = 0.0;
      csol[8*dim+18] = 0.0; csol[8*dim+19] = 0.0;
      fail = solcon_append_quaddobl_solution(dim,1,csol);
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

void lib2path_read_dobldobl_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<dd_real>,dd_real>& start_sys,
   PolySys<complexH<dd_real>,dd_real>& target_sys,
   PolySolSet<complexH<dd_real>,dd_real>& sols )
{
   int fail;

   std::cout << target_file << " " << strlen(target_file) << std::endl;
   std::cout << start_file << " " << strlen(start_file) << std::endl;
   fail = read_dobldobl_target_system_from_file
            (strlen(target_file), target_file);
   lib2path_read_dobldobl_sys(0,target_sys);

   fail = read_dobldobl_start_system_from_file(strlen(start_file),start_file);
   lib2path_read_dobldobl_sys(0,start_sys);
   lib2path_read_dobldobl_sols(start_sys, sols);
}

void lib2path_read_quaddobl_homotopy
 ( char* start_file, char* target_file, \
   PolySys<complexH<qd_real>,qd_real>& start_sys,
   PolySys<complexH<qd_real>,qd_real>& target_sys,
   PolySolSet<complexH<qd_real>,qd_real>& sols)
{
   int fail;
   std::cout << target_file << " " << strlen(target_file) << std::endl;
   std::cout << start_file << " " << strlen(start_file) << std::endl;
   fail = read_quaddobl_target_system_from_file
            (strlen(target_file), target_file);
   lib2path_read_quaddobl_sys(0,target_sys);

   fail = read_quaddobl_start_system_from_file(strlen(start_file),start_file);
   lib2path_read_quaddobl_sys(0,start_sys);
   lib2path_read_quaddobl_sols(start_sys, sols);
}


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

extern "C" int dobldobl_ademanypaths
 ( int verbose, double regamma, double imgamma )
{
   int fail;
   PolySys<complexH<dd_real>,dd_real> ps;
   PolySys<complexH<dd_real>,dd_real> qs;
   PolySolSet<complexH<dd_real>,dd_real> sols;

   fail = copy_dobldobl_target_system_to_container();

   if(verbose > 0)
   {
      int dim, len;

      cout << endl;
      cout << "Tracking many paths ..." << endl;
      cout << "gamma = " << setprecision(16)
           << regamma << " + i* " << imgamma << endl;
      fail = syscon_number_of_dobldobl_polynomials(&dim);
      cout << "number of polynomials : " << dim << endl;
      fail = solcon_number_of_dobldobl_solutions(&len);
      cout << "number of solutions : " << len << endl;
   }

   lib2path_read_dobldobl_sys(verbose,ps);

   fail = copy_dobldobl_start_system_to_container();

   lib2path_read_dobldobl_sys(verbose,qs);
   lib2path_read_dobldobl_sols(qs,sols);

   fail = dobldobl_manytrack(verbose,regamma,imgamma,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   lib2path_write_dobldobl_sols(sols);

   return 0;
}

extern "C" int quaddobl_ademanypaths
 ( int verbose, double regamma, double imgamma )
{
   int fail;
   PolySys<complexH<qd_real>,qd_real> ps;
   PolySys<complexH<qd_real>,qd_real> qs;
   PolySolSet<complexH<qd_real>,qd_real> sols;

   fail = copy_quaddobl_target_system_to_container();

   if(verbose > 0)
   {
      int dim, len;

      cout << endl;
      cout << "Tracking many paths ..." << endl;
      cout << "gamma = " << setprecision(16)
           << regamma << " + i* " << imgamma << endl;
      fail = syscon_number_of_quaddobl_polynomials(&dim);
      cout << "number of polynomials : " << dim << endl;
      fail = solcon_number_of_quaddobl_solutions(&len);
      cout << "number of solutions : " << len << endl;
   }

   lib2path_read_quaddobl_sys(verbose,ps);

   fail = copy_quaddobl_start_system_to_container();

   lib2path_read_quaddobl_sys(verbose,qs);
   lib2path_read_quaddobl_sols(qs,sols);

   fail = quaddobl_manytrack(verbose,regamma,imgamma,ps,qs,sols);

   if(verbose > 0)
      cout << "writing the solutions to the container ..." << endl;

   lib2path_write_quaddobl_sols(sols);

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

int dobldobl_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s )
{
   double tpred,teval,tmgs;
   bool success;
   complexH<dd_real>* sol = s.get_sol(0);
   complexH<dd_real> alpha,t;
   CPUInstHom<complexH<dd_real>,dd_real> cpu_inst_hom;
   Workspace< complexH<dd_real> > workspace_cpu;
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
         cout << setw(40) << scientific << setprecision(32) << sol[k];
      }
   }
   int fail,n_path;
   fail = solcon_number_of_dobldobl_solutions(&n_path);

   alpha = complexH<dd_real>(regamma,imgamma);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      if(verbose > 0) cout << "tracking path " << path_idx << endl;
      complexH<dd_real>* sol0 = s.get_sol(path_idx);
      complexH<dd_real>* sol_new = NULL;
      t = complexH<dd_real>(0.0,0);
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
         cout.precision(32);
         cout << "The solution " << path_idx << endl;
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(40) << workspace_cpu.x_last[k];
      }
   }
   return 0;
}

int quaddobl_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s )
{
   double tpred,teval,tmgs;
   bool success;
   complexH<qd_real>* sol = s.get_sol(0);
   complexH<qd_real> alpha,t;
   CPUInstHom<complexH<qd_real>,qd_real> cpu_inst_hom;
   Workspace< complexH<qd_real> > workspace_cpu;
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
         cout << setw(72) << scientific << setprecision(64) << sol[k];
      }
   }
   int fail,n_path;
   fail = solcon_number_of_quaddobl_solutions(&n_path);

   alpha = complexH<qd_real>(regamma,imgamma);
   cpu_inst_hom.init(p,q,p.dim,p.n_eq,1,alpha,verbose);
   cpu_inst_hom.init_workspace(workspace_cpu);

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      if(verbose > 0) cout << "tracking path " << path_idx << endl;
      complexH<qd_real>* sol0 = s.get_sol(path_idx);
      complexH<qd_real>* sol_new = NULL;
      t = complexH<qd_real>(0.0,0);
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
         cout.precision(64);
         cout << "The solution " << path_idx << endl;
         for(int k=0; k<p.dim; k++)
            cout << k << " :" << setw(72) << workspace_cpu.x_last[k];
      }
   }
   return 0;
}
