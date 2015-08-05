#include "ada_test.h"

void var_name(char* x_string, int x_string_len, string*& x_names, int& dim)
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

void ada_read_sys(PolySys& sys)
{
   int fail,nbsym;

   std::cout << "testing reading and writing a system" << std::endl;
   //fail = syscon_read_system();
   std::cout << "the system is .." << std::endl;
   fail = syscon_write_system();
   fail = syscon_number_of_symbols(&nbsym);
   std::cout << "the number of symbols : " << nbsym << std::endl;
   // Get variable names
   int s_dim = 80*nbsym;
   char *s = (char*) calloc(s_dim,sizeof(char));
   fail = syscon_string_of_symbols(&s_dim, s);
   string* x_names;
   int dim = 0;
   var_name(s, s_dim, x_names, dim);
   int i = 1;
   std::cout << "dim = " << dim << std::endl;

   double c[2];
   int d[dim];

   int n_eq = 0;
   fail = syscon_number_of_polynomials(&n_eq);

   sys.n_eq = n_eq;
   sys.dim  = dim;
   sys.eq_space = new PolyEq[n_eq];
   sys.pos_var = x_names;

   PolyEq* tmp_eq = sys.eq_space;

   for(int i=1; i<n_eq+1; i++)
   {
      int nt;
      fail = syscon_number_of_terms(i,&nt);
      //std::cout << " #terms in polynomial " << i << " : " << nt << std::endl;
      tmp_eq->n_mon = nt;
      tmp_eq->dim = dim;
      for(int j=1; j<=nt; j++)
      {
         fail = syscon_retrieve_term(i,j,dim,d,c);
         //std::cout << c[0] << " " << c[1] << std::endl;
         //for (int k=0; k<n; k++) std::cout << " " << d[k];
         //std::cout << std::endl;
         bool constant_term = true;
         for(int k=0; k<dim; k++)
            if(d[k]!=0) constant_term = false;

         if(constant_term==true)
         {
            tmp_eq->n_mon--;
            tmp_eq->constant += CT(c[0],c[1]);
            //std::cout << "constant " << c[0] << " " << c[1] << std::endl;
         }
         else
         {
            PolyMon* a = new PolyMon(dim,d,c);
            tmp_eq->mon.push_back(a);
         }
      }
      tmp_eq->print(x_names);
      sys.eq.push_back(tmp_eq);
      tmp_eq++;
   }
   sys.print();
   std::cout << "End" << std::endl;
}

void ada_read_sols(PolySys& start_sys, PolySolSet& sols)
{
   int fail, len;
   //fail = copy_start_solutions_to_container();
   fail = solcon_number_of_solutions(&len);
   printf("Number of start solutions : %d\n",len);
   int dim=start_sys.dim;
   sols.dim = dim;
   double sol[2*dim+5];

   for(int sol_idx=1; sol_idx<len+1; sol_idx++)
   {
      int mm;
      int k;
      solcon_retrieve_next_standard_solution(dim,&k,&mm,sol);
      //solcon_retrieve_solution(dim,sol_idx,&mm,sol);
      /*std::cout << sol[0] << " " << sol[1] << std::endl;
      for(int var_idx=0; var_idx<4; var_idx++)
         std::cout << sol[2+2*var_idx] << " " << sol[2+2*var_idx] << std::endl;
      std::cout << sol[2+2*dim] 
                << " " << sol[3+2*dim]
                << " " << sol[4+2*dim] << std::endl;*/
      PolySol* tmp_sol = new PolySol(dim, sol[0], sol[1], sol+2, \
         sol[2+2*dim], sol[3+2*dim], sol[4+2*dim]);
      //tmp_sol->print_info(start_sys.pos_var);
      sols.add_sol(tmp_sol);
   }
   //sols.print_info(start_sys.pos_var);
   std::cout << "sol finished" << std::endl;
}

void ada_read_homotopy(char* start_file, char* target_file, \
		PolySys& start_sys, PolySys& target_sys, PolySolSet& sols)
{
   // adainit();
   int fail;
   std::cout << target_file << " " << strlen(target_file) << std::endl;
   std::cout << start_file << " " << strlen(start_file) << std::endl;
   fail = read_standard_target_system_from_file
            (strlen(target_file), target_file);
   ada_read_sys(target_sys);

   fail = read_standard_start_system_from_file(strlen(start_file),start_file);
   ada_read_sys(start_sys);
   ada_read_sols(start_sys, sols);
   // adafinal();
}
