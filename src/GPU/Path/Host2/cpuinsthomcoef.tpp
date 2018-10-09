// The file cpuinsthomcoef.tpp provides the definitions of the methods
// of the class CPUInstHomCoef, defined in the file cpuinsthomcoef.h.

template <class ComplexType, class RealType>
void CPUInstHomCoef<ComplexType,RealType>::init
 ( MonSet<ComplexType>* hom_monset, int total_n_mon, int n_monset,
   int n_constant, ComplexType alpha, int verbose )
{
   this->alpha = alpha;
   n_coef = total_n_mon;
   coef_orig = new ComplexType[n_coef*2];
   ComplexType* tmp_coef_orig = coef_orig;

   int constant_exist = 0;
   if(n_constant > 0) constant_exist = 1;

   // write start coefficient and target coefficient together
   for(int i=constant_exist; i<n_monset; i++)
      hom_monset[i].write_coef(tmp_coef_orig);

   if(n_constant > 0)
      hom_monset[0].write_coef(tmp_coef_orig);

   // write start coefficient and target coefficient seperately
   tmp_coef_orig = new ComplexType[n_coef*2];
   for(int coef_idx=0; coef_idx<n_coef; coef_idx++)
   {
      tmp_coef_orig[coef_idx] = coef_orig[2*coef_idx];
      tmp_coef_orig[coef_idx+n_coef] = coef_orig[2*coef_idx+1];
   }
   delete[] coef_orig;
   coef_orig = tmp_coef_orig;
}

template <class ComplexType, class RealType>
void CPUInstHomCoef<ComplexType,RealType>::print()
{
   for(int i=0; i<n_coef; i++)
   {
      std::cout << i << std::endl
                << coef_orig[i]
                << coef_orig[i+n_coef]<< std::endl;
   }
}

template <class ComplexType, class RealType>
void CPUInstHomCoef<ComplexType,RealType>::eval
 ( const ComplexType t, ComplexType* coef, int reverse )
{
   ComplexType one_minor_t(1.0- t.real, -t.imag);

   int k = 1;
   ComplexType t_power_k = t;
   ComplexType one_minor_t_power_k = one_minor_t;
   for(int i=1; i<k; i++)
   {
      t_power_k *= t;
      one_minor_t_power_k *= one_minor_t;
   }
   ComplexType t0, t1;
   if(reverse == 0)
   {
      t0 = one_minor_t_power_k*alpha;
      t1 = t_power_k;
   }
   else
   {
      t0 = t_power_k*alpha;
      t1 = one_minor_t_power_k;
   }
   for(int i=0; i<n_coef; i++)
      coef[i] = coef_orig[i+n_coef]*t0 + coef_orig[i]*t1;
}

template <class ComplexType, class RealType>
void CPUInstHomCoef<ComplexType,RealType>::update_alpha ( ComplexType alpha )
{
   if(alpha.real == 0 && alpha.imag == 0)
   {
      int r = rand();
      RealType tmp = RealType(r);
      this->alpha = ComplexType(sin(tmp),cos(tmp));
   }
   else
      this->alpha = alpha;
}
