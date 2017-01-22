#include "polysol.h"

bool compare_sol ( PolySol* sol1, PolySol* sol2 )
{
   if(*sol1 < *sol2)
      return true;
   else
      return false;
}

void PolySol::init
 ( CT* sol, int dim, T1 max_residual, T1 max_delta_x, int path_idx,
   string path_info )
{
   this->dim = dim;
   this->sol = new CT[dim];
   for(int i=0; i<dim; i++)
   {
      this->sol[i] = sol[i];
   }
   this->path_idx = path_idx;
   idx = 0;
   m = 0;
   t = CT(0.0,0.0);
   err = max_delta_x;
   rco = 0.0;
   res = max_residual;
   info = path_info;
}

void PolySol::init
 ( int dim, T1 t_real, T1 t_imag, T1* sol, 
   T1 max_delta_x, T1 rco, T1 max_residual,
   int m, int path_idx, string path_info )
{
   this->dim = dim;
   this->sol = new CT[dim];
   for(int i=0; i<dim; i++)
   {
      this->sol[i].real = sol[2*i];
      this->sol[i].imag = sol[2*i+1];
   }
   this->path_idx = path_idx;
   idx = 0;
   this->m = m;
   t = CT(t_real,t_imag);
   err = max_delta_x;
   this->rco = rco;
   res = max_residual;
   info = path_info;
}

void PolySol::init ( ifstream& myfile, int dim )
{
   this->dim = dim;
   string tmp_line;

   /*getline(myfile,tmp_line);
     std::cout << tmp_line << std::endl;
     std::cout << "dim = " << dim << std::endl;*/

   // Get idx
   getline(myfile,tmp_line, ' ');
   // std::cout << tmp_line << std::endl;
   myfile >> idx;
   // std::cout << "idx = " << idx << std::endl;
   getline(myfile,tmp_line);

   int tmp_line_l = tmp_line.size();
   if(tmp_line_l > 2)
   {
       unsigned found = tmp_line.find_last_of(" ");
       tmp_line = tmp_line.substr(found+1);
       tmp_line_l = tmp_line.size();
       if(tmp_line[tmp_line_l-1]<'a' or tmp_line[tmp_line_l-1]>'z')
       {
          tmp_line = tmp_line.substr(0,tmp_line_l-1);
       }
       info = tmp_line;
   }
   // Get t
   getline(myfile,tmp_line, ':');
   // std::cout << tmp_line << std::endl;
   t = get_complex_number(myfile);
   // std::cout << t << std::endl;
   getline(myfile,tmp_line);

   // Get m
   getline(myfile,tmp_line, ':');
   // std::cout << tmp_line << std::endl;
   int m;
   myfile >> m;
   // std::cout << "m = " << m << std::endl;
   getline(myfile,tmp_line);

   getline(myfile,tmp_line);
   // std::cout << tmp_line << std::endl;

   sol = new CT[dim];

   for(int i=0; i<dim; i++)
   {
      getline(myfile,tmp_line, ':');
      sol[i] = get_complex_number(myfile);
      // std::cout<< sol[i];
      getline(myfile,tmp_line);
   }
   // Get error
   getline(myfile,tmp_line, ':');
   myfile >> err;
   // std::cout << "error = " << err << endl;

   getline(myfile,tmp_line, ':');
   myfile >> rco;
   // std::cout << "rco = " << rco << endl;

   getline(myfile,tmp_line, ':');
   myfile >> res;
   // std::cout << "res = " << res << endl;
   getline(myfile,tmp_line);
}

void PolySol::init ( ifstream& myfile, int dim, VarDict& pos_dict )
{
   this->dim = dim;
   string tmp_line;

   /*getline(myfile,tmp_line);
     std::cout << tmp_line << std::endl;
     std::cout << "dim = " << dim << std::endl;*/

   // Get idx
   getline(myfile,tmp_line, ' ');
   // std::cout << tmp_line << std::endl;
   myfile >> idx;
   // std::cout << "idx = " << idx << std::endl;
   getline(myfile,tmp_line);

   int tmp_line_l = tmp_line.size();
   if(tmp_line_l > 2)
   {
      unsigned found = tmp_line.find_last_of(" ");
      tmp_line = tmp_line.substr(found+1);
      tmp_line_l = tmp_line.size();
      if(tmp_line[tmp_line_l-1]<'a' or tmp_line[tmp_line_l-1]>'z')
      {
         tmp_line = tmp_line.substr(0,tmp_line_l-1);
      }
      info = tmp_line;
   }
   // Get t
   getline(myfile,tmp_line, ':');
   // std::cout << tmp_line << std::endl;
   t = get_complex_number(myfile);
   // std::cout << t << std::endl;
   getline(myfile,tmp_line);

   // Get m
   getline(myfile,tmp_line, ':');
   // std::cout << tmp_line << std::endl;
   int m;
   myfile >> m;
   // std::cout << "m = " << m << std::endl;
   getline(myfile,tmp_line);

   getline(myfile,tmp_line);
   // std::cout << tmp_line << std::endl;

   sol = new CT[dim];

   // pos_dict.print();
   for(int i=0; i<dim; i++)
   {
      getline(myfile,tmp_line, ':');
      // std::cout << tmp_line << tmp_line.substr(1,tmp_line.length()-2);
      int var_idx = pos_dict.get(tmp_line.substr(1,tmp_line.length()-2));
      // std::cout << var_idx << std::endl;
      sol[var_idx] = get_complex_number(myfile);
      // std::cout<< sol[i];
      getline(myfile,tmp_line);
   }
   // Get error
   getline(myfile,tmp_line, ':');
   myfile >> err;
   // std::cout << "error = " << err << endl;

   getline(myfile,tmp_line, ':');
   myfile >> rco;
   // std::cout << "rco = " << rco << endl;

   getline(myfile,tmp_line, ':');
   myfile >> res;
   // std::cout << "res = " << res << endl;
   getline(myfile,tmp_line);
}

bool PolySol::operator == ( const PolySol& that )
{
   T1 dis = 0;
   for(int i=0; i<dim; i++)
   {
      T1 tmp_dis = abs(sol[i].real - that.sol[i].real);
      if(tmp_dis > dis)
      {
         dis = tmp_dis;
      }
      tmp_dis = abs(sol[i].imag - that.sol[i].imag);
      if(tmp_dis > dis)
      {
         dis = tmp_dis;
      }
   }
   if(dis > 1E-6)
   {
      return 0;
   }
   return 1;
}

void PolySol::print()
{
   std::cout << "dim = " << dim << std::endl;
   for(int i=0; i<dim; i++)
      std::cout << i << " " << sol[i];
   std::cout << std::endl;
}

void PolySol::print_info()
{
   std::cout << "Solution " << idx << ":" << std::endl;
   std::cout << "t : " << t;
   for(int i=0; i<dim; i++)
   {
      std::cout << i << " : " << sol[i];
   }
   std::cout << "== err : " << err \
             << " = rco : " << rco \
             << " = res : " << res \
             << " ==" << std::endl;
}

void PolySol::print_info ( string* pos_var )
{
   std::cout << "Solution " << idx << ":" << std::endl;
   std::cout << "t : " << t;
   for(int i=0; i<dim; i++)
   {
      std::cout << pos_var[i] << " : " << sol[i];
   }
   std::cout << "== err : " << err
             << " = rco : " << rco
             << " = res : " << res
             << " ==" << std::endl;
}

CT* PolySol::get_sol()
{
   CT* sol_tmp = new CT[dim];

   for(int i=0; i<dim; i++) sol_tmp[i] = sol[i];

   return sol_tmp;
}

void PolySol::print_short()
{
   std::cout << std::scientific
             << " " << err  << " " << res << " " << " x[0] = " << sol[0];
}

int to_int ( double a )
{
   int b = a;
   return b;
}

/*int to_int ( dd_real a )
{
   int b = a.x[0];
   return b;
}

int to_int ( qd_real a )
{
   int b = a.x[0];
   return b;
}*/

bool PolySol::operator< ( PolySol& that )
{
   if(dim < that.dim)
   {
      return true;
   }
   if(dim > that.dim)
   {
      return false;
   }
   CT* this_sol = new CT[dim];
   CT* that_sol = new CT[dim];

   for(int i=0; i<dim; i++)
   {
      this_sol[i] = sol[i];
      that_sol[i] = that.sol[i];
   }
   int digits_per_check = 2;

   T1 digits_multiplier = 1;
   T1 err_roundoff = 1E-4;

   for(int i=0; i<digits_per_check; i++)
   {
      digits_multiplier *= 10;
   }
   int digits_checked = 0;

   while(digits_checked < 8)
   {
      for(int i=0; i<dim; i++)
      {
         /*std::cout << i << std::endl;
           std::cout << this_sol[i].real << std::endl;
           std::cout << that_sol[i].real << std::endl;*/
         int this_digit;
         int that_digit;
         int tmp_number;

         // real check
         this_sol[i].real *= digits_multiplier;
         that_sol[i].real *= digits_multiplier;
         /*std::cout << digits_multiplier << std::endl;
           std::cout << this_sol[i].real << std::endl;
           std::cout << that_sol[i].real << std::endl;*/

         this_digit = to_int(this_sol[i].real);
         tmp_number = to_int(this_sol[i].real+err_roundoff);
         if(this_digit != tmp_number) this_digit = tmp_number;

         tmp_number = to_int(this_sol[i].real-err_roundoff);
         if(this_digit != tmp_number) this_digit = tmp_number;

         that_digit = to_int(that_sol[i].real);
         tmp_number = to_int(that_sol[i].real+err_roundoff);
         if(that_digit != tmp_number) that_digit = tmp_number;

         tmp_number = to_int(that_sol[i].real-err_roundoff);
         if(that_digit != tmp_number) that_digit = tmp_number;

         // std::cout << this_digit << " " << that_digit<< std::endl;

         if(this_digit < that_digit)
         {
            return true;
         }
         else if(this_digit > that_digit)
         {
            return false;
         }
         this_sol[i].real -= this_digit;
         that_sol[i].real -= that_digit;

         // imag check
         this_sol[i].imag *= digits_multiplier;
         that_sol[i].imag *= digits_multiplier;
         this_digit = to_int(this_sol[i].imag);
         tmp_number = to_int(this_sol[i].imag+err_roundoff);
         if(this_digit != tmp_number) this_digit = tmp_number;

         tmp_number = to_int(this_sol[i].imag-err_roundoff);
         if(this_digit != tmp_number) this_digit = tmp_number;

         that_digit = to_int(that_sol[i].imag);
         tmp_number = to_int(that_sol[i].imag+err_roundoff);
         if(that_digit != tmp_number) that_digit = tmp_number;

         tmp_number = to_int(that_sol[i].imag-err_roundoff);
         if(that_digit != tmp_number) that_digit = tmp_number;

         // std::cout << this_digit << " " << that_digit<< std::endl;
         if(this_digit < that_digit)
         {
            return true;
         }
         else if(this_digit > that_digit)
         {
            return false;
         }
         this_sol[i].imag -= this_digit;
         that_sol[i].imag -= that_digit;
      }
      digits_checked += digits_per_check;
   }
   delete[] this_sol;
   delete[] that_sol;

   return false;
}

void PolySolSet::init ( ifstream& myfile )
{
   n_sol = 0;
   dim = 0;

   string prefix_two[2] = {"START SOLUTIONS : ", "THE SOLUTIONS :"};
   read_until_line(myfile, prefix_two, 2);

   myfile >> n_sol;
   myfile >> dim;

   string prefix = "======";
   read_until_line(myfile, prefix);

   sols.reserve(n_sol);

   //std::cout << "n_sol = " << n_sol << std::endl;
   for(int i=0; i<n_sol; i++)
   {
      // std::cout << i << std::endl;
      PolySol* tmp_sol= new PolySol(myfile, dim);
      sols.push_back(tmp_sol);
   }
}

void PolySolSet::init ( ifstream& myfile, VarDict& pos_dict )
{
   n_sol = 0;
   dim = 0;

   string prefix_two[2] = {"START SOLUTIONS : ", "THE SOLUTIONS :"};
   read_until_line(myfile, prefix_two, 2);

   myfile >> n_sol;
   myfile >> dim;

   string prefix = "======";
   read_until_line(myfile, prefix);

   sols.reserve(n_sol);

   // std::cout << "n_sol = " << n_sol << std::endl;
   for(int i=0; i<n_sol; i++)
   {
      // std::cout << i << std::endl;
      PolySol* tmp_sol= new PolySol(myfile, dim, pos_dict);
      sols.push_back(tmp_sol);
   }
}

bool PolySolSet::find_same_sol ( PolySol* tmp_sol )
{
   for(int i=0; i<n_sol; i++)
   {
      if(*tmp_sol == *sols[i])
      {
         return true;
      }
   }
   return false;
}

int PolySolSet::count_same_sol ( PolySol* tmp_sol )
{
   int n_same_sol = 0;
   for(int i=0; i<n_sol; i++)
   {
      if(*tmp_sol == *sols[i])
      {
         n_same_sol++;
      }
   }
   return n_same_sol;
}

bool PolySolSet::add_diff_sol ( CT* new_sol )
{
   PolySol* tmp_sol= new PolySol(new_sol, dim);
   if(find_same_sol(tmp_sol)==true)
   {
      return false;
   }
   std::cout << "Add New Solution" << std::endl;
   sols.push_back(tmp_sol);
   n_sol++;
   return true;
}

void PolySolSet::add_sol ( PolySol* tmp_sol )
{
   sols.push_back(tmp_sol);
   n_sol++;
}

void PolySolSet::change_sol ( int idx, CT* coords )
{
   PolySol* idxsol = sols[idx];
   for(int k=0; k<dim; k++) idxsol->sol[k] = coords[k];
}

void PolySolSet::add_sol
 ( CT* new_sol, T1 max_residual, T1 max_delta_x, int path_idx,
   string path_info )
{
   PolySol* tmp_sol = new PolySol(new_sol, dim, max_residual, max_delta_x, 
                                  path_idx, path_info);
   add_sol(tmp_sol);
}


void PolySolSet::print()
{
   std::cout << "dim   = " << dim << std::endl
             << "n_sol = " << n_sol << std::endl;

   for(int i=0; i<n_sol; i++) sols[i]->print();
}

void PolySolSet::print_info ( string* pos_var )
{
   std::cout << "dim   = " << dim << std::endl
             << "n_sol = " << n_sol << std::endl;

   for(int i=0; i<n_sol; i++) sols[i]->print_info(pos_var);
}

void PolySolSet::print_short()
{
   std::cout << "n_sol = " << n_sol << std::endl;
   for(int i=0; i<n_sol; i++)
      std::cout << i << " " << sols[i]->info << " x[0] = " << sols[i]->sol[0];
}

CT* PolySolSet::get_sol ( int idx )
{
   return sols[idx]->get_sol();
}

void PolySolSet::sort_set()
{
   sort(sols.begin(), sols.end(), compare_sol);
}

void PolySolSet::compare ( PolySolSet& that )
{
   vector<PolySol*> that_sols(that.sols);
   vector<PolySol*> this_sols(this->sols);

   sort(this_sols.begin(), this_sols.end(), compare_sol);
   sort(that_sols.begin(), that_sols.end(), compare_sol);

   vector<PolySol*>::iterator this_pointer = this_sols.begin();
   vector<PolySol*>::iterator that_pointer = that_sols.begin();

   vector<PolySol*> this_sols_only;
   vector<PolySol*> that_sols_only;

   int same_sol = 0;
   int this_sol = 0;
   int that_sol = 0;

   while(true)
   {
      if((**this_pointer)==(**that_pointer))
      {
         this_pointer++;
         that_pointer++;
         same_sol++;
      }
      else
      {
         if((**this_pointer)<(**that_pointer))
         {
            this_sols_only.push_back(*this_pointer);
            this_pointer++;
            this_sol++;
         }
         else
         {
            /* std::cout << this_pointer - this_sols.begin() << std::endl;
               std::cout << that_pointer - that_sols.begin() << std::endl;
               (*this_pointer)->print();
               (*that_pointer)->print();
             */
            that_sols_only.push_back(*that_pointer);
            that_pointer++;
            that_sol++;
         }
      }
      if(this_pointer==this_sols.end() || that_pointer==that_sols.end())
      {
         break;
      }
   }

   if(this_pointer != this_sols.end())
   {
      this_sol += this_sols.end() - this_pointer;
      while(this_pointer < this_sols.end())
      {
         this_sols_only.push_back(*this_pointer);
         this_pointer++;
      }
   }

   if(that_pointer != that_sols.end())
   {
      that_sol += that_sols.end() - that_pointer;
      while(that_pointer < that_sols.end())
      {
         that_sols_only.push_back(*that_pointer);
         that_pointer++;
      }
   }
   std::cout << "same_sol = " << same_sol << std::endl;

   if(this_sol > 0) std::cout << "this_sol = " << this_sol << std::endl;
   if(that_sol > 0) std::cout << "that_sol = " << that_sol << std::endl;

   that_pointer = that_sols_only.begin();
   int i=0;
   while(that_pointer < that_sols_only.end())
   {
      bool find_same_this = find_same_sol(*that_pointer);
      int that_n_same_sol = that.count_same_sol(*that_pointer);
      std::cout << std::setw(4) << i++ << " "
                << std::setw(8) << (*that_pointer)->path_idx << " "
                << std::setw(4) << find_same_this << " "
                << std::setw(4) << that_n_same_sol;
      (*that_pointer)->print_short();
      that_pointer++;
   }
}
