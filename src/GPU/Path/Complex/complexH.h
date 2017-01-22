// The class complexH define complex arithmetic on the host,
// the real type is a templated parameter T.

#ifndef __COMPLEXH_H__
#define __COMPLEXH_H__

#include <iostream>
#include <qd/qd_real.h>
#include <stdlib.h>

using namespace std;

#include <iostream>
#include <fstream>

template <class T>
class complexH
{
   public:

      complexH() {};

      complexH operator= ( const complexH& );
      complexH operator= ( const int& );
      complexH operator= ( const double& );
      complexH operator= ( const dd_real& );
      complexH operator= ( const qd_real& );

      complexH ( const double&, const double& );
      complexH ( const dd_real&, const dd_real& );
      complexH ( const qd_real&, const qd_real& );

      complexH ( const int& );
      complexH ( const double& );

      complexH ( const char* );
      complexH ( const char*, const char* );
      complexH ( const complexH& );

      void init ( const double&, const double& );
      // define a complex number with a tuple of doubles

      void init ( const dd_real&, const dd_real& );
      // define a complex number with a tuple of double doubles

      void init ( const qd_real&, const qd_real& );
      // define a complex number with a tuple of quad doubles

      // operators for complex arithmetic are below
      complexH operator+ ( const complexH& );
      void     operator+= ( const complexH& );
      complexH operator- ( const complexH& );
      void     operator-= ( const complexH& );
      complexH operator* ( const complexH& ) const;
      complexH operator* ( const T& );
      complexH operator* ( const int& );
      void     operator*= ( const complexH& );
      complexH operator/ ( const complexH& );
      complexH operator/ ( const T& );

      T absv() {return sqrt(real*real+imag*imag);};
      // returns the modulus or radius of a complex number

      complexH adj() {return complexH(real,-imag);};
      // complex conjugate of complex number

      T real; // real part of complex number
      T imag; // imaginary part of complex number

      template <class T1>
      friend std::ostream& operator<<
         ( std::ostream& os, const complexH<T1>& number );

      template <class T1>
      friend std::ifstream& operator>>
         ( std::ifstream& is, const complexH<T1>& number );
};

#include "complexH.cpp"
#endif
