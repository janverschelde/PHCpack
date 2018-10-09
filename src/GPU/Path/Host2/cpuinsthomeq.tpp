// The file cpuinsthomeq.tpp provides definitions for the methods in
// the class CPUInstHomEq, with prototypes in the file cpuinsthomeq.h.

template <class ComplexType>
void CPUInstHomEq<ComplexType>::init
 ( MonSet<ComplexType>* hom_monset, int n_monset, int n_eq, int n_constant )
{
   this->n_eq = n_eq;
   n_mon_eq = NULL;
   eq_pos_start = NULL;
   n_mon_total = 0;
   mon_pos_start_eq = NULL;
   n_pos_total = 0;
   mon_pos_eq = NULL;
   coef = NULL;

   n_mon_eq = new int[n_eq];
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++) n_mon_eq[eq_idx] = 0;

   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      if(hom_monset[set_idx].get_n()!=0)
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
            int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
            n_mon_eq[eq_idx]++;
         }
      }
   }
   eq_pos_start = new int[n_eq];
   int tmp_pos=0;
   eq_pos_start[0] = 0;

   for(int eq_idx=0; eq_idx<n_eq-1; eq_idx++)
      eq_pos_start[eq_idx+1] = eq_pos_start[eq_idx] + (n_mon_eq[eq_idx]+1);

   n_mon_total = eq_pos_start[n_eq-1] + n_mon_eq[n_eq-1]+1;

   int* eq_mon_pos_size = new int[n_eq];
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++) eq_mon_pos_size[eq_idx] = 0;

   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      if(hom_monset[set_idx].get_n()!=0)
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
            int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
            eq_mon_pos_size[eq_idx] += hom_monset[set_idx].get_n() + 1;
         }
      }
   }
   n_pos_total = 0;

   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
      n_pos_total += eq_mon_pos_size[eq_idx];

   mon_pos_start_eq = new int[n_mon_total];
   coef = new ComplexType[n_mon_total*2];
   mon_pos_start_eq[0] = 0;

   int* eq_mon_idx_tmp = new int[n_eq];
   int* eq_mon_pos_tmp = new int[n_eq];
   eq_mon_pos_tmp[0] = 0;
   eq_mon_idx_tmp[0] = 0;
   for(int eq_idx=0; eq_idx<n_eq-1; eq_idx++)
   {
      eq_mon_pos_tmp[eq_idx+1] = eq_mon_pos_tmp[eq_idx]
                               + eq_mon_pos_size[eq_idx];
   }
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      eq_mon_idx_tmp[eq_idx] = eq_pos_start[eq_idx];
   }
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      mon_pos_start_eq[eq_pos_start[eq_idx]] = n_mon_eq[eq_idx];
      eq_mon_idx_tmp[eq_idx]++;
   }
   // init coef
   for(int eq_idx=0; eq_idx<n_eq; eq_idx++)
   {
      coef[2*eq_pos_start[eq_idx]] = ComplexType(0.0,0.0);
      coef[2*eq_pos_start[eq_idx]+1] = ComplexType(0.0,0.0);
   }
   mon_pos_eq = new unsigned short[n_pos_total];
   for(int set_idx=0; set_idx<n_monset; set_idx++)
   {
      if(hom_monset[set_idx].get_n()!=0)
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
            int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
            mon_pos_start_eq[eq_mon_idx_tmp[eq_idx]] = eq_mon_pos_tmp[eq_idx];
            ComplexType* coef_tmp = coef + 2*eq_mon_idx_tmp[eq_idx];
            hom_monset[set_idx].write_coef(coef_tmp, mon_idx);
            int tmp_pos = eq_mon_pos_tmp[eq_idx];
            mon_pos_eq[tmp_pos] = hom_monset[set_idx].get_n();
            unsigned short* tmp_pointer = mon_pos_eq+tmp_pos+1;
            hom_monset[set_idx].write_pos(tmp_pointer);
            eq_mon_pos_tmp[eq_idx] += hom_monset[set_idx].get_n() + 1;
            eq_mon_idx_tmp[eq_idx]++;
         }
      }
      else
      {
         int n_mon_set = hom_monset[set_idx].get_n_mon();
         for(int mon_idx=0; mon_idx<n_mon_set; mon_idx++)
         {
             int eq_idx = hom_monset[set_idx].get_eq_idx(mon_idx);
             ComplexType* coef_tmp = coef + 2*eq_pos_start[eq_idx];
             hom_monset[set_idx].write_coef(coef_tmp, mon_idx);
         }
      }
   }
}

template <class ComplexType>
void CPUInstHomEq<ComplexType>::print()
{
   std::cout << "CPUInstHomMonEq" << std::endl;
   std::cout << "n_eq = " << n_eq << std::endl;
   std::cout << "n_mon       = " << n_mon_total-n_eq << std::endl;
   std::cout << "n_mon+n_eq  = " << n_mon_total << std::endl;
   std::cout << "n_pos       = " << n_pos_total-(n_mon_total-n_eq) << std::endl;
   std::cout << "n_pos+n_mon = " << n_pos_total << std::endl;
}
