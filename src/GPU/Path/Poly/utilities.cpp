# include "utilities.h"

int log2ceil(int n)
{
   n = 2*n-1;
   int log2n = 0;
   while(n>1)
   {
      log2n++;
      n/=2;
   }
   return log2n;
}

/**
 * Stop string at certain single symbol.
 * Used by get_coef_complex
 */

void string_stop
 ( const string& mon_string, int& loc, const char* symbols,
   const int n_symbol, int l )
{
   for(int i = loc; i<l; i++)
   {
      char tmp = mon_string[i];
      int stop = 1;
      for(int j=0; j<n_symbol; j++)
      {
         if(tmp == symbols[j])
         {
            loc++;
            stop = 0;
            break;
         }
      }
      if(stop == 1)
      {
         break;
      }
   }
}

/**
 * Read number string from string.
 * Used by get_coef_complex
 */

string get_number_string( const string& mon_string, int& loc, int l )
{
   string number_string;

   // Get symbol
   if((mon_string[loc] == '+') || (mon_string[loc] == '-'))
   {
       number_string += mon_string[loc];
       loc++;
   }
   // Get number
   int read_to_end = 1;
   for(int i=loc, coef_exp = 0; i<l; i++)
   {
      char tmp = mon_string[i];
      if(coef_exp && (tmp == '+' || tmp == '-'))
      {
         number_string += tmp;
         coef_exp = 0;
      }
      else if((tmp>'9' || tmp<'.') && tmp!='e' && tmp!='E' && tmp != ' ')
      {
         loc = i;
         read_to_end = 0;
         break;
      }
      else if(tmp!= ' ')
      {
         number_string += tmp;
         if(tmp == 'E' || tmp == 'e')
         {
            coef_exp = 1;
         }
      }
   }
   if(read_to_end == 1) loc = l;

   return number_string;
}

CT get_complex_number ( ifstream& myfile )
{
   T1 tmp_real, tmp_imag;
   myfile >> tmp_real;
   myfile >> tmp_imag;
   return CT(tmp_real, tmp_imag);
}

void read_until_line ( ifstream& myfile, string prefix )
{
   string tmp_line;
   while(getline(myfile,tmp_line))
   {
      if(tmp_line.substr(0,prefix.size()) == prefix)
      {
         break;
      }
   }
   //std::cout << tmp_line << std::endl;
}

void read_until_line ( ifstream& myfile, string* prefix, int n_prefix )
{
   string tmp_line;
   while(getline(myfile,tmp_line))
   {
      for(int i=0; i<n_prefix; i++)
      {
         if(tmp_line.substr(0,prefix[i].size()) == prefix[i])
         {
            return;
         }
      }
   }
   std::cout << tmp_line << std::endl;
}

/**
 * Read complex number from string.
 */

CT get_coef_complex ( const string& mon_string, int& loc )
{
   int l = mon_string.length();
   CT coef(0.0,0.0);

   const char skip_symbols[] = {'+','(', ' '};
   int n_skip_symbols = 3;
   string_stop(mon_string, loc, skip_symbols, n_skip_symbols, l);

   string coef_string;

   while(true)
   {
      const char skip_symbols2[] = {' '};
      int n_skip_symbols2 = 1;
      string_stop(mon_string, loc, skip_symbols2, n_skip_symbols2, l);

      // number part
      coef_string = get_number_string(mon_string, loc, l);

      if(coef_string.length() == 0)
      {
         coef = CT(1.0);
         break;
      }
      else if(coef_string.length()== 1)
      {
         if(coef_string[0] == '+')
         {
            coef = CT(1.0);
            break;
         }
         else if(coef_string[0] == '-')
         {
            coef = CT(-1.0);
            break;
         }
      }
      // real or imag
      if( ((loc < l) && (mon_string[loc]=='i'))\
         || ((loc+1 < l) && (mon_string[loc]=='*' and mon_string[loc+1]=='i')))
      {
         coef.imag += read_number(coef_string.c_str());
         loc += 2;
      }
      else
      {
         coef.real += read_number(coef_string.c_str());
      }
      // break when get to the end or meet other symbol
      if(loc >= l)
      {
         break;
      }
      else
      {
         char tmp = mon_string[loc];
         if((tmp>'9' || tmp<'0') && tmp != '.' 
             && tmp!='e' && tmp!='E' && tmp!='+' && tmp!='-')
         {
            break;
         }
      }
   }

   const char skip_symbols3[] = {')', ' ', '*', ';'};
   int n_skip_symbols3 = 4;
   string_stop(mon_string, loc, skip_symbols3, n_skip_symbols3, l);

   // std::cout << "real " << coef_string_real << std::endl;
   // std::cout << "imag " << coef_string_imag << std::endl;

   return coef;
}

CT get_coef_real ( const string& mon_string, int& loc )
{
   int l = mon_string.length();

   int start = loc;
   CT coef;
   loc = l;
   string coef_string;

   if((mon_string[start] == '+') || (mon_string[start] == '-'))
   {
      coef_string += mon_string[start];
      start++;
   }

   for(int i=start, coef_exp = 0; i<l; i++)
   {
      char tmp = mon_string[i];
      if(coef_exp && (tmp == '+' || tmp == '-'))
      {
         coef_string += tmp;
         coef_exp = 0;
      }
      else if((tmp>'9' || tmp<'.') && tmp!='E' && tmp != ' ')
      {
         loc = i;
         break;
      }
      else if(tmp!= ' ')
      {
         coef_string += tmp;
         if(tmp == 'E')
         {
            coef_exp = 1;
         }
      }
   }
   //cout << coef_string << endl;

   if(coef_string.length() == 0)
   {
      coef = CT(1.0);
   }
   else if(coef_string.length()== 1)
   {
      if(coef_string[0] == '+')
      {
         coef = CT(1.0);
      }
      else if(coef_string[0] == '-')
      {
         coef = CT(-1.0);
      }
      else
      {
         coef = CT(coef_string.c_str());
         //cout << coef_string << " " << coef <<endl;
      }
   }
   else
   {
      coef = CT(coef_string.c_str());
      //cout << coef_string << " " << coef <<endl;
   }
   for(int i = loc; i<l; i++)
   {
      char tmp = mon_string[i];
      if(tmp == ' ')
      {
         loc++;
      }
      else if (tmp == '*')
      {
         loc++;
      }
      else if (tmp == ';')
      {
         loc++; break;
      }
      else
      {
         break;
      }
   }
   return coef;
}

template <class T1>
void print_matrix ( T1** m, int row, int col )
{
   // int width = 6;
   for(int i=0; i<row; i++)
   {
      for(int j=0; j<col; j++)
      {
         std::cout << m[i][j];
      }
      std::cout << endl;
   }
};

template <class T1>
void print_vector ( T1* v, int n )
{
   // int width = 6;
   for(int i=0; i<n; i++)
   {
      std::cout << v[i];
   }
   std::cout << endl;
};


template <class T1>
void print_result ( T1* v, T1** m, int dim, int n_eq )
{
   std::cout << "f =" << endl;
   print_vector<T1>(v, n_eq);
   std::cout << "deri ="<< endl;
   print_matrix<T1>(m, n_eq, dim);
};

void print_size ( size_t tmp_size_B )
{
   size_t tmp_size_KB = tmp_size_B/1000;
   size_t tmp_size_MB = tmp_size_KB/1000;
   size_t tmp_size_GB = tmp_size_MB/1000;
   if(tmp_size_GB != 0)
   {
      tmp_size_B -= 1000*tmp_size_KB;
      tmp_size_KB -= 1000*tmp_size_MB;
      tmp_size_MB -= 1000*tmp_size_GB;
      std::cout << setw(3)<< tmp_size_GB << " GB ";
      std::cout << setw(3)<<tmp_size_MB << " MB ";
      std::cout << setw(3)<<tmp_size_KB << " KB ";
      std::cout << setw(3)<<tmp_size_B << " B" << std::endl;
   }
   else if(tmp_size_MB > 0)
   {
      tmp_size_B -= 1000*tmp_size_KB;
      tmp_size_KB -= 1000*tmp_size_MB;
      std::cout <<"       ";
      std::cout <<setw(3)<< tmp_size_MB << " MB ";
      std::cout <<setw(3)<<tmp_size_KB << " KB ";
      std::cout <<setw(3)<< tmp_size_B << " B" << std::endl;
   }
   else if(tmp_size_KB > 0)
   {
      tmp_size_B -= 1000*tmp_size_KB;
      std::cout <<"       ";
      std::cout <<"       ";
      std::cout << setw(3)<<tmp_size_KB << " KB ";
      std::cout << setw(3)<<tmp_size_B << " B" << std::endl;
   }
   else
   {
      std::cout <<"       ";
      std::cout <<"       ";
      std::cout <<"       ";
      std::cout << setw(3)<<tmp_size_B << " B" << std::endl;
   }
};

double* rand_val ( int dim )
{
   double* val = new double[dim];
   for(int i=0; i< dim; i++)
   {
      val[i] = i+1;
      // val[i] = (double)rand() / (double)RAND_MAX;
   }
   return val;
}

CT* rand_val_complex ( int dim )
{
   CT* val = new CT[dim];
   T1 zero = 0.0;
   for(int i=0; i< dim; i++)
   {
      // T1 tmp = T1(rand())/T1(RAND_MAX);
      // val[i] = CT(double(i+1), 0.0);
      T1 tmp = T1(1)/T1(i+2);
      T1 tmp1 = T1(i+2)/T1(i+3);
      val[i] = CT(tmp, tmp1);
   }
   return val;
}

CT* rand_val_complex_frac ( int dim )
{
   CT* val = new CT[dim];
   T1 zero = 0.0;
   for(int i=0; i< dim; i++)
   {
      T1 tmp = T1(1)/T1(i+2);
      T1 tmp1 = T1(i+2)/T1(i+3);
      val[i] = CT(tmp, tmp1);
   }
   return val;
}

CT* rand_val_complex_one ( int dim )
{
   CT* val = new CT[dim];
   T1 zero = 0.0;
   for(int i=0; i< dim; i++)
   {
      val[i] = CT(1.0, 0.0);
   }
   return val;
}


CT* rand_val_complex_n ( int dim )
{
   CT* val = new CT[dim];
   T1 zero = 0.0;
   for(int i=0; i< dim; i++)
   {
      val[i] = CT(double(i+1), 0.0);
   }
   return val;
}

CT* rand_val_complex_unit(int dim)
{
   CT* val = new CT[dim];
   T1 zero = 0.0;
   for(int i=0; i<dim; i++)
   {
      int r = rand();
      T1 tmp = T1(r);
      val[i].init(sin(tmp),cos(tmp));
   }
   return val;
}

CT* rand_val_complex_unit_n( int dim )
{
   CT* val = new CT[dim];
   T1 zero = 0.0;
   for(int i=0; i<dim; i++)
   {
      T1 tmp = T1(i+1);
      val[i].init(sin(tmp),cos(tmp));
   }
   return val;
}

string* var_list ( int dim, string v )
{
   string* vars = new string[dim];
   for(int i=0; i<dim; i++)
   {
      std::stringstream tmp;
      tmp << v << i;
      vars[i] = tmp.str();
   }
   return vars;
}

void print_number ( T1 v )
{
   if(v > 0)
   {
      cout << " + " << v;
   }
   else if(v < 0)
   {
      cout << " - " << abs(v);
   }
}

void print_number_complex ( CT v )
{
   T1 v_real = v.real;
   T1 v_imag = v.imag;
   if(v_real != 0 && v_imag!= 0)
   {
      cout << " + (";
      print_number(v_real);
      print_number(v_imag);
      cout << "*i )";
   }
   else if(v_real != 0)
   {
      print_number(v_real);
   }
   else if(v.imag != 0)
   {
      print_number(v_imag);
   }
}

void print_coef_complex ( CT coef )
{
   if(coef.real != 1.0 or coef.imag != 0.0)
   {
      //cout << coef;
      print_number_complex(coef);
      cout <<  " * ";
   }
   else
   {
      cout << " + ";
   }
}

void cpu_speel0
 ( const CT* x_val, unsigned short* pos, CT* deri, const CT& coef )
{
   //int n_var = pos[0];
   deri[0] = coef*x_val[pos[1]];
   deri[1] = coef;
}

void cpu_speel_with_base0
 ( const CT* x_val, unsigned short* pos, unsigned short* exp, CT* deri,
   const CT& coef )
{
   //int n_var = pos[0];
   deri[0] = coef*x_val[pos[1]]*exp[1];
   deri[1] = coef;
}

// speelpanning on CPU
void cpu_speel
 ( const CT* x_val, unsigned short* pos, CT* deri, const CT& coef )
{
   int n_var = pos[0];
   CT tmp = x_val[pos[1]];

   //std::cout << pos[1] << " ";
   CT* deri_tmp = deri + 1;
   deri_tmp[1] = tmp;
   for(int i=2; i<n_var; i++)
   {
      tmp *=x_val[pos[i]];
      deri_tmp[i] = tmp;
      //std::cout << pos[i] << " ";
   }
   /*if(n_var > 1){
        std::cout << pos[n_var] << std::endl;
   }*/

   tmp = coef;
   //CT tmp2 = x_val[pos[n_var]];

   for(int i=n_var; i>1; i--)
   {
      deri[i] *= tmp;
      tmp *= x_val[pos[i]];
   }
   deri[1] = tmp;
   deri[0] = x_val[pos[1]]*tmp;

   //std::cout << "deri[0] = " << deri[0];
}

// speelpanning on CPU
void cpu_speel_with_base
 ( const CT* x_val, unsigned short* pos, unsigned short* exp, CT* deri,
   const CT& coef)
{
   int n_var = pos[0];
   CT tmp = x_val[pos[1]];

   //std::cout << pos[1] << " ";
   CT* deri_tmp = deri + 1;
   deri_tmp[1] = tmp;
   for(int i=2; i<n_var; i++)
   {
      tmp *=x_val[pos[i]];
      deri_tmp[i] = tmp;
      //std::cout << pos[i] << " ";
   }
   /*if(n_var > 1){
        std::cout << pos[n_var] << std::endl;
   }*/

   tmp = coef;
   //CT tmp2 = x_val[pos[n_var]];
   int base_start = exp[0];
   for(int i=n_var; i>base_start; i--)
   {
      deri[i] *= tmp*exp[i];
      tmp *= x_val[pos[i]];
   }
   for(int i=base_start; i>1; i--)
   {
      deri[i] *= tmp;
      tmp *= x_val[pos[i]];
   }
   deri[1] = tmp*exp[1];
   deri[0] = x_val[pos[1]]*tmp;

   //std::cout << "deri[0] = " << deri[0];
}

double time_interval1
 ( struct timeval start, struct timeval end )
{
   long seconds, useconds;   // mtime,
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   return seconds*1000+ (double)useconds/1000;
}
