// The file cpuinsthomcoef.h defines the class CPUInstHomCoef.

#ifndef __CPUINSTHOMCOEF_H__
#define __CPUINSTHOMCOEF_H__

template <class ComplexType, class RealType>
class CPUInstHomCoef
{
   public:

      int n_coef;
      ComplexType* coef_orig;
      ComplexType alpha;

      CPUInstHomCoef()
      {
         n_coef = 0;
         coef_orig = NULL;
         alpha = ComplexType(0.0,0);
      }

      CPUInstHomCoef ( MonSet<ComplexType>* hom_monset,
                       int total_n_mon, int n_monset,
                       int n_constant, ComplexType alpha )
      {
         init(hom_monset, total_n_mon, n_monset, n_constant, alpha);
      }

      ~CPUInstHomCoef()
      {
         // std::cout << "Delete CPUInstHomCoef" << std::endl;
         if(coef_orig != NULL) delete[] coef_orig;
      }

      void init ( MonSet<ComplexType>* hom_monset,
                  int total_n_mon, int n_monset,
                  int n_constant, ComplexType alpha, int verbose = 0 );

      void print();

      void eval ( const ComplexType t, ComplexType* coef, int reverse=0 );

      void update_alpha ( ComplexType alpha=ComplexType(0.0,0.0) );
      /*
       * If alpha on input is zero,
       * then a new random value for alpha will be generated,
       * else the value on input for alpha will be stored.
       */
};

#include "cpuinsthomcoef.tpp"

#endif /* __CPUINSTHOMCOEF_H__ */
