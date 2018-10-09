// The file eqidxcoef.h defines the class EqIdxCoef to store two
// complex constants and an equation index.

#ifndef __EQIDXCOEF_H__
#define __EQIDXCOEF_H__

#include<iostream>

template <class ComplexType>
class EqIdxCoef
{
   int eq_idx;
   ComplexType coef[2];

   void init( const int eq_idx, const ComplexType& coef0,
              const ComplexType& coef1 )
   {
      this->eq_idx = eq_idx;
      coef[0] = coef0;
      coef[1] = coef1;
   }

   public:

      EqIdxCoef()
      {
         eq_idx = 0;
      }

      EqIdxCoef ( const int eq_idx, const ComplexType& coef0,
                  const ComplexType& coef1 )
      {
         init(eq_idx, coef0, coef1);
      }

      EqIdxCoef& operator= ( const EqIdxCoef& original )
      {
         init(original.eq_idx, original.coef[0], original.coef[1]);
         return *this;
      }

      EqIdxCoef( const int eq_idx, const ComplexType& coef, const bool sys )
      {
         this->eq_idx = eq_idx;
         this->coef[sys] = coef;
         this->coef[(!sys)] = 0.0;
      }

      void print()
      {
         std::cout << "Eq " << eq_idx << std::endl;
         std::cout << coef[0];
         std::cout << coef[1];
      }

      void write_coef ( ComplexType*& tmp_coef )
      {
         *tmp_coef++ = coef[0];
         *tmp_coef++ = coef[1];
      }

      int get_eq_idx()
      {
         return eq_idx;
      }

      friend std::ostream& operator << ( std::ostream& o, const EqIdxCoef& c )
      {
         return o << "Eq " << c.eq_idx << std::endl 
                  << c.coef[0] << c.coef[1];
      }
};

#endif /* __EQIDXCOEF_H__ */
