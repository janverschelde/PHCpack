// The templated definitions are included in polysol.h.

template <class ComplexType, class RealType>
bool compare_sol
 ( PolySol<ComplexType,RealType>* sol1, PolySol<ComplexType,RealType>* sol2 )
{
   if(*sol1 < *sol2)
      return true;
   else
      return false;
}

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::init
 ( ComplexType* sol, int dim, RealType max_residual,
   RealType max_delta_x, int path_idx, string path_info )
{
   this->dim = dim;
   this->sol = new ComplexType[dim];
   for(int i=0; i<dim; i++)
   {
      this->sol[i] = sol[i];
   }
   this->path_idx = path_idx;
   idx = 0;
   m = 0;
   t = ComplexType(0.0,0.0);
   err = max_delta_x;
   rco = 0.0;
   res = max_residual;
   info = path_info;
}

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::init
 ( int dim, RealType t_real, RealType t_imag, RealType* sol, 
   RealType max_delta_x, RealType rco, RealType max_residual,
   int m, int path_idx, string path_info )
{
   this->dim = dim;
   this->sol = new ComplexType[dim];
   for(int i=0; i<dim; i++)
   {
      this->sol[i].real = sol[2*i];
      this->sol[i].imag = sol[2*i+1];
   }
   this->path_idx = path_idx;
   idx = 0;
   this->m = m;
   t = ComplexType(t_real,t_imag);
   err = max_delta_x;
   this->rco = rco;
   res = max_residual;
   info = path_info;
}

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::init ( ifstream& myfile, int dim )
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
   t = get_complex_number<ComplexType,RealType>(myfile);
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

   sol = new ComplexType[dim];

   for(int i=0; i<dim; i++)
   {
      getline(myfile,tmp_line, ':');
      sol[i] = get_complex_number<ComplexType,RealType>(myfile);
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

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::init
 ( ifstream& myfile, int dim, VarDict& pos_dict )
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
   t = get_complex_number<ComplexType,RealType>(myfile);
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

   sol = new ComplexType[dim];

   // pos_dict.print();
   for(int i=0; i<dim; i++)
   {
      getline(myfile,tmp_line, ':');
      // std::cout << tmp_line << tmp_line.substr(1,tmp_line.length()-2);
      int var_idx = pos_dict.get(tmp_line.substr(1,tmp_line.length()-2));
      // std::cout << var_idx << std::endl;
      sol[var_idx] = get_complex_number<ComplexType,RealType>(myfile);
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

template <class ComplexType, class RealType>
bool PolySol<ComplexType,RealType>::operator == ( const PolySol& that )
{
   RealType dis = 0;
   for(int i=0; i<dim; i++)
   {
      RealType tmp_dis = abs(sol[i].real - that.sol[i].real);
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

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::print()
{
   std::cout << "dim = " << dim << std::endl;
   for(int i=0; i<dim; i++)
      std::cout << i << " " << sol[i];
   std::cout << std::endl;
}

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::print_info()
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

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::print_info ( string* pos_var )
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

template <class ComplexType, class RealType>
ComplexType* PolySol<ComplexType,RealType>::get_sol()
{
   ComplexType* sol_tmp = new ComplexType[dim];

   for(int i=0; i<dim; i++) sol_tmp[i] = sol[i];

   return sol_tmp;
}

template <class ComplexType, class RealType>
void PolySol<ComplexType,RealType>::print_short()
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

template <class ComplexType, class RealType>
bool PolySol<ComplexType,RealType>::operator< ( PolySol& that )
{
   if(dim < that.dim)
   {
      return true;
   }
   if(dim > that.dim)
   {
      return false;
   }
   ComplexType* this_sol = new ComplexType[dim];
   ComplexType* that_sol = new ComplexType[dim];

   for(int i=0; i<dim; i++)
   {
      this_sol[i] = sol[i];
      that_sol[i] = that.sol[i];
   }
   int digits_per_check = 2;

   RealType digits_multiplier = 1;
   RealType err_roundoff = 1E-4;

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
