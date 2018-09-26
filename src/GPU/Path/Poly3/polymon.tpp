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
 ( const string& mon_string, VarDict& pos_dict )
{
   int end = mon_string.length();
   int loc = 0;
   ComplexType coef = get_coef_complex<ComplexType>(mon_string, loc);
   read(mon_string, pos_dict, loc, end, coef);
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::read
 ( const string& mon_string, VarDict& pos_dict, ComplexType coef )
{
   int end = mon_string.length();
   read(mon_string, pos_dict, 0, end, coef);
}

template <class ComplexType, class RealType>
void PolyMon<ComplexType,RealType>::read
 ( const string& eq_string, VarDict& pos_dict, int start, int end,
   ComplexType coef )
{
   n_var = 1;

   this->coef = coef;

   for(int i = start; i<end; ++i)
      if(eq_string[i] == '*')
      {
          if(eq_string[i+1] == '*')
             i++;
          else
             n_var++;
      }

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
            next_type = 2;
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
         next_type = 2;
         new_var++;
      }
      else if(c == ' ' || c ==';')
      {
         new_var++;
      }

      // handle vars
      if(new_var == 1)
      {
         if(cur_type == 1)
         {
            // pos[pos_ind] = var_to_pos(var);
            // string var1 = var.substr(1,var.length()-1);
            // cout << var << " " << var1 << " " << atoi(var1.c_str()) << endl;
            // pos[pos_ind] = atoi(var1.c_str())-1;
            pos[pos_ind] = pos_dict.get(var);
            exp[pos_ind] = 1;
            pos_ind++;
            cur_type = 0;
         }
         else if(cur_type == 2)
         {
            int tmp_exp =  atoi(var.c_str());
            exp[pos_ind-1] = tmp_exp;
            cur_type = 0;
         }
      }
      else
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
   }

   if(cur_type == 1)
   {
      string var1 = var.substr(1,var.length()-1);
      // cout << var << " " << var1 << " " << atoi(var1.c_str()) << endl;
      // pos[pos_ind] = var_to_pos(var);
      // pos[pos_ind] = atoi(var1.c_str())-1;
      pos[pos_ind] = pos_dict.get(var);
      exp[pos_ind] = 1;
   }
   else if(cur_type == 2)
   {
      exp[pos_ind-1] = atoi(var.c_str()); 
   }
   update_base();
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

