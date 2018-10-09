// The file cpuinsthommon.tpp provides the definitions for the methods
// of the class CPUInstHomMon, with prototypes in cpuinsthommon.h.

template <class ComplexType>
void CPUInstHomMon<ComplexType>::init
 ( MonSet<ComplexType>* hom_monset, int n_monset, int total_n_mon,
   int n_constant, int* max_deg_base, int verbose )
{
   this->max_deg_base = max_deg_base;
   int max_n_var = hom_monset[n_monset-1].get_n();
   level = 1;
   for(int i=1; i<max_n_var; i<<=1) level++;
   n_mon_level = new int[level];
   for(int i=0; i<level; i++) n_mon_level[i] = 0;
   mon_pos_size = 0;
   n_mon = total_n_mon - n_constant;
   if(verbose > 0)
      std::cout << "n_constant = " << n_constant << std::endl;
   int constant_exist;
   if(n_constant == 0)
      constant_exist = 0;
   else
      constant_exist = 1;
   // Write monomial start position
   int tmp_level = 0;
   int tmp_level_size = 1;
   mon_pos_start = new int[n_mon];

   int mon_idx = 0;
   for(int i=constant_exist; i<n_monset; i++)
   {
      int tmp_n = hom_monset[i].get_n();
      int tmp_n_mon = hom_monset[i].get_n_mon();
      if(tmp_n > 1)
      {
         flops_multiple += tmp_n_mon*3*(tmp_n-1);
      }
      else
      {
         flops_multiple += tmp_n_mon;
      }
      for(int j=0; j<tmp_n_mon; j++)
      {
         mon_pos_start[mon_idx++] = mon_pos_size;
         mon_pos_size += tmp_n+1;
      }
      while(tmp_level_size < tmp_n)
      {
         tmp_level_size *= 2;
         tmp_level++;
      }
      n_mon_level[tmp_level]+=tmp_n_mon;
   }
   // Write position instruction
   mon_pos = new unsigned short[mon_pos_size];
   mon_exp = new unsigned short[mon_pos_size];

   unsigned short* tmp_pos = mon_pos;
   unsigned short* tmp_exp = mon_exp;
   for(int i=constant_exist; i<n_monset; i++)
   {
      // Number of variable each term
      int tmp_n_mon = hom_monset[i].get_n_mon();
      for(int j=0; j<tmp_n_mon; j++)
      {
         int n_var = hom_monset[i].get_n();
         *tmp_pos++ = n_var;
         int base_start = hom_monset[i].get_base_start();
         *tmp_exp++ = base_start;
         // Write position
         hom_monset[i].write_pos(tmp_pos);
         hom_monset[i].write_exp(tmp_exp);
      }
   }
   n_mon_base_start = n_mon;
   for(int mon_idx=0; mon_idx<n_mon; mon_idx++)
   {
      int tmp_mon_pos_start = mon_pos_start[mon_idx];
      int tmp_n_mon = mon_pos[tmp_mon_pos_start];
      if(mon_exp[tmp_mon_pos_start+tmp_n_mon]>1)
      {
         n_mon_base_start = mon_idx;
         break;
      }
   }
   if(verbose > 0)
   {
      std::cout << "*** flops_multiple = " << flops_multiple << std::endl;
      std::cout << "*** flops          = "
                << flops_multiple*4+flops_multiple*2 << std::endl;
   }
}

template <class ComplexType>
void CPUInstHomMon<ComplexType>::eval
 ( int dim, const ComplexType* x_val, ComplexType* mon,
   ComplexType* coef, ComplexType** deg_table )
{
   if(deg_table != NULL)
   {
      eval_deg_table(dim, x_val, deg_table);
      eval_base(deg_table, coef);
      for(int j=0; j<n_mon_level[0]; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel_with_base0(x_val,mon_pos+tmp_idx,mon_exp+tmp_idx,
            mon+tmp_idx, coef[j]);
      }
      for(int j=n_mon_level[0]; j<n_mon; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel_with_base(x_val,mon_pos+tmp_idx,mon_exp+tmp_idx,
            mon+tmp_idx, coef[j]);
      }
   }
   else
   {
      for(int j=0; j<n_mon_level[0]; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel0(x_val,mon_pos+tmp_idx,mon+tmp_idx,coef[j]);
      }
      for(int j=n_mon_level[0]; j<n_mon; j++)
      {
         int tmp_idx = mon_pos_start[j];
         cpu_speel(x_val,mon_pos+tmp_idx,mon+tmp_idx,coef[j]);
      }
   }
}

template <class ComplexType>
void CPUInstHomMon<ComplexType>::eval_deg_table
 ( int dim, const ComplexType* x_val, ComplexType** deg_table )
{
   for(int var_idx=0; var_idx<dim; var_idx++)
   {
      if(max_deg_base[var_idx]>0)
      {
         ComplexType* tmp_deg_table = deg_table[var_idx];
         ComplexType tmp_var = x_val[var_idx];
         tmp_deg_table[0] = tmp_var;
         for(int deg_idx=1; deg_idx<max_deg_base[var_idx]; deg_idx++)
         {
            tmp_deg_table[deg_idx] = tmp_deg_table[deg_idx-1]*tmp_var;
         }
      }
   }
}

template <class ComplexType>
void CPUInstHomMon<ComplexType>::eval_base
 ( ComplexType** deg_table, ComplexType* coef )
{
   for(int mon_idx=n_mon_base_start; mon_idx<n_mon; mon_idx++)
   {
      int tmp_idx = mon_pos_start[mon_idx];
      unsigned short* tmp_mon_pos = mon_pos + tmp_idx;
      unsigned short* tmp_mon_exp = mon_exp + tmp_idx;
      int tmp_var_start = *tmp_mon_exp++;
      int tmp_n_var = *tmp_mon_pos++;
      if(tmp_var_start < tmp_n_var)
      {
         ComplexType tmp_val = coef[mon_idx];
         for(int var_idx=tmp_var_start; var_idx<tmp_n_var; var_idx++)
            tmp_val *= deg_table[tmp_mon_pos[var_idx]][tmp_mon_exp[var_idx]-2];
         coef[mon_idx] = tmp_val;
      }
   }
}

template <class ComplexType>
void CPUInstHomMon<ComplexType>::print()
{
   std::cout << "level = " << level << std::endl;
   for(int i=0; i<level; i++)
   {
      std::cout << i << " " << n_mon_level[i] << std::endl;
   }
   std::cout << "n_mon = " << n_mon << std::endl;
   // Print monomial with position
   for(int i=0; i<n_mon; i++)
   {
      int tmp_n = mon_pos[mon_pos_start[i]];
      int tmp_exp_start = mon_exp[mon_pos_start[i]];
      std::cout << i << " n = " << tmp_n << ": ";
      for(int j=0; j<tmp_n; j++)
         std::cout << mon_pos[mon_pos_start[i]+j+1] << ", ";
      std::cout << " exp_start = " << tmp_exp_start << ": ";
      for(int j=0; j<tmp_n; j++)
         std::cout << mon_exp[mon_pos_start[i]+j+1] << ", ";
      std::cout << " mon_pos_start = " << mon_pos_start[i] << std::endl;
   }
}
