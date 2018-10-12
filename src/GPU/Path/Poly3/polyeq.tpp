template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval ( const ComplexType* x_val )
{
   ComplexType val = constant;

   for(int i=0; i<n_mon; i++) val += mon[i]->eval(x_val);

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri, ComplexType *monderi )
{
   for(int i=0; i<dim; i++) deri[i].init(0.0,0.0);

   ComplexType val = constant;
    
   for(int i=0; i<n_mon; i++)
   {
      PolyMon<ComplexType,RealType>* m = mon[i];
      val += m->eval(x_val, monderi);

      for(int j=0; j<m->n_var; j++) deri[m->pos[j]] += monderi[j];
   }

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri )
{
   ComplexType* mon_deri = new ComplexType[dim];

   ComplexType val = eval(x_val, deri, mon_deri);

   delete[] mon_deri;

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri, ComplexType* monderi,
   ComplexType** deg_table )
{
   for(int i=0; i<dim; i++) deri[i].init(0.0,0.0);

   ComplexType val = constant;

   for(int mon_idx=0; mon_idx<n_mon; mon_idx++)
   {
      PolyMon<ComplexType,RealType>* m = mon[mon_idx];
      val += m->eval(x_val, monderi, deg_table);
      for(int j=0; j<m->n_var; j++) deri[m->pos[j]] += monderi[j];
   }

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri, ComplexType** deg_table )
{
   ComplexType* mon_deri = new ComplexType[dim];

   ComplexType val = eval(x_val,deri,mon_deri,deg_table);

   delete [] mon_deri;

   return val;
}

template <class ComplexType, class RealType>
void PolyEq<ComplexType,RealType>::read
 ( const string& eq_string, VarDict& pos_dict, int verbose )
{
   if(verbose > 0) cout << "entering the first PolyEq.read()" << endl;
   int l = eq_string.length();
   read(eq_string, pos_dict, 0, l, verbose);
   if(verbose > 0) cout << "... leaving the first PolyEq.read()" << endl;
}

template <class ComplexType, class RealType>
void PolyEq<ComplexType,RealType>::read
 ( const string& eq_string, VarDict& pos_dict, int start, int end,
   int verbose )
{
   if(verbose > 0) cout << "entering the second PolyEq.read()" << endl;

   n_mon = 0;

   // Get the starting position of first monomial
   int mon_first = start;
   for(int i=start; i< end; i++)
   {
      char c = eq_string[i];
      if(c != ' ')
      {
         mon_first = i;
         break;
      }
   }
   // Count the number of monomials
   // Generate monomials link list, include coef, start and end pos in string

   int mon_start = mon_first;
   int i=mon_first;
   int parentheses_open = 0;
   while(i<end)
   {
      char c = eq_string[i];
      if(c == '(')
      {
         parentheses_open = 1;
      }
      else if (c == ')')
      {
         parentheses_open = 0;
      }
      // If a new monomial appears,
      // record the starting and ending position of last one.
      if(((c== '+' || c == '-') && eq_string[i-1]!='e' && eq_string[i-1]!='E' \
         && eq_string[i-1]!='(') && (i != mon_first) && parentheses_open == 0)
      {
         int tmp_start = mon_start;
         ComplexType tmp_coef 
           = get_coef_complex<ComplexType>(eq_string, mon_start);
         parentheses_open = 0;

         if(mon_start < i)
         {
            n_mon++;
            // std::cout << "n_mon = " << n_mon << std::endl;
            PolyMon<ComplexType,RealType>* mm
               = new PolyMon<ComplexType,RealType>;
            mm->read(eq_string, pos_dict, mon_start, i, tmp_coef, verbose);
            mon.push_back(mm);
            mon_start = i;
         }
         else
         {
            // std::cout << tmp_coef;
            constant += tmp_coef;
         }
      }
      i++;
   }
   if(mon_start < end)
   {
      ComplexType tmp_coef
        = get_coef_complex<ComplexType>(eq_string, mon_start);
      if(mon_start < end)
      {
         n_mon++;
         PolyMon<ComplexType,RealType>* mm = new PolyMon<ComplexType,RealType>;
         mm->read(eq_string, pos_dict, mon_start, end, tmp_coef, verbose);
         mon.push_back(mm);
      }
      else
      {
         constant += tmp_coef;
      }
   }
   dim = pos_dict.n_job;

   if(verbose > 0) cout << "... leaving the second PolyEq.read()" << endl;
}

template <class ComplexType, class RealType>
void PolyEq<ComplexType,RealType>::print(const string* pos_var)
{
   // std::cout << "n_mon = " << n_mon << std::endl;
   for(int i=0; i<n_mon; i++) mon[i]->print(pos_var);

   print_number_complex<ComplexType,RealType>(constant);
   cout << endl;
}

template <class ComplexType, class RealType>
void PolyEq<ComplexType,RealType>::update_max_deg ( int* max_deg )
{
   for(int i=0; i<n_mon; i++)
   {
      // cout << "mon " << i << endl;
      mon[i]->update_max_deg(max_deg);
   }
}
