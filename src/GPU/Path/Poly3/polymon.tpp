// The file polymon.tpp provides the definitions of the methods of the
// class PolyMon, with prototypes in the file polymon.h.

inline int var_to_pos ( string var )
{
   string var1 = var.substr(1,var.length()-1);
   // return atoi(var1.c_str());
   // return atoi(var1.c_str())-1;
   return (atoi(var1.c_str())+6)%8;
}

template <class ComplexType>
ComplexType pow ( const ComplexType x, int exp )
{
   ComplexType tmp = x;

   for(int i=1; i<exp; i++) tmp *= x;

   return tmp;
}

int not_empty_line ( const string& line )
{
   int not_empty = -1;
   for(unsigned int j=0; j < line.length(); j++)
   {
      char c = line[j];

      if(c!=' ' && c!=';' && c!='\n') not_empty = j;
   }
   return not_empty;
}

std::vector<std::string> &split
 ( const std::string &s, char delim, std::vector<std::string> &elems )
{
   std::stringstream ss(s);
   std::string item;
   std::string empy_str = "";

   while (std::getline(ss, item, delim))
      if(item != empy_str) elems.push_back(item);

   return elems;
}

std::vector<std::string> split ( const std::string &s, char delim )
{
   std::vector<std::string> elems;
   split(s, delim, elems);
   return elems;
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::read
 ( const string& mon_string, VarDict& pos_dict, int verbose )
{
   if(verbose > 0) cout << "entering the first PolyMon.read()" << endl;

   int end = mon_string.length();
   int loc = 0;
   ComplexType coef = get_coef_complex<ComplexType>(mon_string, loc);
   read(mon_string, pos_dict, loc, end, coef, verbose);

   if(verbose > 0) cout << "... leaving the first PolyMon.read()" << endl;
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::read
 ( const string& mon_string, VarDict& pos_dict, ComplexType coef, int verbose )
{
   if(verbose > 0) cout << "entering the second PolyMon.read()" << endl;
   int end = mon_string.length();
   read(mon_string, pos_dict, 0, end, coef, verbose);
   if(verbose > 0) cout << "... leaving the second PolyMon.read()" << endl;
}

void bubblesortinfo ( int index, int* positions, int* exponents )
/*
 * Writes the information in positions and exponents,
 * for when the verbose flag in bubblesort is raised.
 */
{
   cout << "-> the positions :"; 
   for(int posidx=0; posidx <= index; posidx++)
      cout << " " << positions[posidx];
   cout << endl;
   cout << "-> the exponents :"; 
   for(int posidx=0; posidx <= index; posidx++)
      cout << " " << exponents[posidx];
   cout << endl;
}

void bubblesort ( int index, int* positions, int* exponents, int verbose=0 )
/*
 * Applies bubble sort to insert the element at position index.
 *
 * ON ENTRY :
 *   index     last index in positions and exponents, index >= 0;
 *   positions is a sorted array, sorted up to position index-1;
 *   exponents contains the exponents of the variables at positions;
 *   verbose   if > 0, then extra output is written to screen.
 *
 * ON RETURN :
 *   positions is a sorted array, sorted up to position index;
 *   exponents contains the exponents of the variables at positions.
 */
{
   if(verbose > 0)
   {
      cout << "Entering bubblesort with index = " << index << endl;
      bubblesortinfo(index,positions,exponents);
   }
   if(index > 0)
   {
      int posidx = index;
      bool unsorted = true;
      do
      {
         if(positions[posidx] > positions[posidx-1])
            unsorted = false;
         else // if(positions[posidx] < positions[posidx-1])
         {
            int tmp = positions[posidx];
            positions[posidx] = positions[posidx-1];
            positions[posidx-1] = tmp;
            tmp = exponents[posidx];
            exponents[posidx] = exponents[posidx-1];
            exponents[posidx-1] = tmp;
            posidx--;
         }
      }
      while(posidx > 0 and unsorted);
   }
   if(verbose > 0)
   {
      cout << "Leaving bubblesort with index = " << index << endl;
      bubblesortinfo(index,positions,exponents);
   }
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::read
 ( const string& eq_string, VarDict& pos_dict, int start, int end,
   ComplexType coef, int verbose )
{
   if(verbose > 0) cout << "entering the third PolyMon.read()" << endl;

   n_var = 1;

   this->coef = coef;

   for(int i = start; i<end; ++i)   // counts the number of variables
      if(eq_string[i] == '*')
      {
          if(eq_string[i+1] == '*') // ** is exponentiation
             i++;
          else
             n_var++;
      }

   if(verbose > 0) cout << "n_var : " << n_var << endl;

   pos = new int[n_var];
   exp = new int[n_var];

   string var;
   int cur_type = 0, pos_ind = 0;

   for(int i=start, next_type=1; i<end; ++i)
   {
      char c = eq_string[i];
      int new_var = 0;
      // handle symbols
      if(c == '*')
      {
         if(eq_string[i+1] == '*')
         {
            next_type = 2; // exponent of type **
            i++;
         }
         else
         {
            next_type = 1;
         }
         new_var++;
      }
      else if(c=='^')
      {
         next_type = 2;   // exponent of type ^
         new_var++;
      }
      else if(c == ' ' || c ==';')
      {
         new_var++;
      }

      if(verbose > 0)
         cout << "  new_var : " << new_var
              << "  cur_type : " << cur_type
              << "  next_type : " << next_type
              << "  pos_ind : " << pos_ind << endl;

      // handle vars
      if(new_var == 1)
      {
         if(cur_type == 1)
         {
            if(verbose > 0)
            {
               string var1 = var.substr(1,var.length()-1);
               cout << var << " " << var1 << " " << atoi(var1.c_str()) << endl;
               cout << "assiging to pos_ind : " << pos_ind << endl;
               cout << "pos_dict.get(var) : " << pos_dict.get(var) << endl;
               cout << "next_type : " << next_type << endl;
            }
            // pos[pos_ind] = var_to_pos(var);
            // pos[pos_ind] = atoi(var1.c_str())-1;
            // In a sparse representation of the monomial, the position
            // of the variable is not used in its position.
            pos[pos_ind] = pos_dict.get(var);
            exp[pos_ind] = 1;
            if(next_type == 1) bubblesort(pos_ind,pos,exp,verbose);
            pos_ind++;
            cur_type = 0;
         }
         else if(cur_type == 2)
         {
            int tmp_exp = atoi(var.c_str());
            if(verbose > 0)
            {
               cout << "assigning exponent " << tmp_exp
                    << " at position " << pos_ind-1 << endl;
               cout << "next_type : " << next_type << endl;
            }
            exp[pos_ind-1] = tmp_exp;
            bubblesort(pos_ind-1,pos,exp,verbose);
            cur_type = 0;
         }
      }
      else // new_var == 0
      {
         if(cur_type == 0)
         {
            var = c;
            cur_type = next_type;            
         }
         else
         {
            var += c;
         }
      }        
   } // end for loop

   if(cur_type == 1)
   {
      if(verbose > 0)
      {
         string var1 = var.substr(1,var.length()-1);
         cout << var << " " << var1 << " " << atoi(var1.c_str()) << endl;
         cout << "assiging to pos_ind : " << pos_ind << endl;
         cout << "pos_dict.get(var) : " << pos_dict.get(var) << endl;
      }
      pos[pos_ind] = pos_dict.get(var);
      exp[pos_ind] = 1;
      bubblesort(pos_ind,pos,exp,verbose);
   }
   else if(cur_type == 2)
   {
      exp[pos_ind-1] = atoi(var.c_str()); 
      if(verbose > 0)
      {
         cout << "assigning exponent " << atoi(var.c_str())
              << " at position " << pos_ind-1 << endl;
      }
      bubblesort(pos_ind-1,pos,exp,verbose);
   }
   update_base();

   if(verbose > 0) 
   {
      int dim = pos_dict.n_job;
      cout << "the monomial read : ";
      this->print_tableau(dim);
      cout << endl;
      cout << "pos :";
      for(int idx=0; idx<n_var; idx++) cout << " " << pos[idx];
      cout << endl;
      cout << "exp :";
      for(int idx=0; idx<n_var; idx++) cout << " " << exp[idx];
      cout << endl;
   }

   if(verbose > 0) cout << "... leaving the third PolyMon.read()" << endl;
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::update_base()
{
   n_base = 0;
   for(int var_idx=0; var_idx<n_var; var_idx++)
   {
      if(exp[var_idx]>1) n_base++;
   }
   pos_base = new int[n_base];
   exp_base = new int[n_base];
   exp_tbl_base = new int[n_base];
   int base_idx = 0;
   for(int var_idx=0; var_idx<n_var; var_idx++)
   {
      if(exp[var_idx]>1)
      {
         pos_base[base_idx] = pos[var_idx];
         exp_base[base_idx] = exp[var_idx];
         exp_tbl_base[base_idx] = exp[var_idx]-2;
         base_idx++;
      }
   }
}

template <class ComplexType, class RealType>
ComplexType PolyMon<ComplexType,RealType>::speel
 ( const ComplexType* x_val, ComplexType* deri )
{
   deri[1] = x_val[pos[0]];

   for(int i=1; i<n_var-1; i++) deri[i+1] = deri[i]*x_val[pos[i]];

   ComplexType tmp = coef;
   for(int i=n_var-1; i>0; i--)
   {
      deri[i] *= tmp;
      tmp *= x_val[pos[i]];
   }
   deri[0] = tmp;

   return tmp*x_val[pos[0]];
}

template <class ComplexType, class RealType>
ComplexType PolyMon<ComplexType,RealType>::speel_with_base
 ( const ComplexType* x_val, ComplexType* deri, ComplexType base )
{
   deri[1] = x_val[pos[0]];

   for(int i=1; i<n_var-1; i++) deri[i+1] = deri[i]*x_val[pos[i]];

   ComplexType tmp = base;
   for(int i=n_var-1; i>0; i--)
   {
      if(exp[i]>1)
         deri[i] *= tmp*exp[i];
      else
         deri[i] *= tmp;
    	
      tmp *= x_val[pos[i]];
   }
   deri[0] = tmp*exp[0];

   return tmp*x_val[pos[0]];
}

template <class ComplexType, class RealType>
ComplexType PolyMon<ComplexType,RealType>::eval ( const ComplexType* x_val )
{
   ComplexType val = coef;

   for(int i=0; i<n_var; i++) val *= pow(x_val[pos[i]],exp[i]);

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyMon<ComplexType,RealType>::eval_base
 ( const ComplexType* x_val, ComplexType** deg_table )
{
   ComplexType val = coef;

   for(int i=0; i<n_base; i++)
      val *= deg_table[pos_base[i]][exp_tbl_base[i]];

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyMon<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri )
{
   ComplexType val = speel(x_val, deri);
   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyMon<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri, ComplexType** deg_table )
{
   ComplexType base = eval_base(x_val, deg_table);
   ComplexType val = speel_with_base(x_val, deri, base);

   return val;
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::print ( const string* pos_var )
{
   // print coef
   print_coef_complex<ComplexType,RealType>(coef);
   // cout << endl;
   cout << pos_var[pos[0]];
   if(exp[0]!= 1) cout << '^' << exp[0];
   for(int i=1; i<n_var; i++)
   {
      cout << " * " << pos_var[pos[i]];
      if(exp[i]!= 1) cout << '^' << exp[i];
   }
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::print_tableau ( int dim )
{
   cout << coef.real << "  " << coef.imag;

   int varidx = 0;
   for(int posidx=0; posidx<n_var; posidx++)
   {
      while(varidx < pos[posidx])
      {
         cout << " 0";
         varidx = varidx+1;
      }
      cout << " " << exp[posidx];
      varidx = varidx+1;
   }
   while(varidx < dim)
   {
      cout << " 0";
      varidx = varidx + 1;
   }
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::update_max_deg ( int* max_deg )
{
   for(int i=0; i<n_var; i++)
   {
      // cout << "var " << i << endl;
      if(exp[i]>max_deg[pos[i]]) max_deg[pos[i]] = exp[i];
   }
}
