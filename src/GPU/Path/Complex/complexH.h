#ifndef __COMPLEXH_H__
#define __COMPLEXH_H__

/*#define __CUDAC__

#ifdef __CUDAC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif*/

#include <iostream>
#include <qd/qd_real.h>
#include <stdlib.h>

using namespace std;

/**
 * The class complexH defines a complex type for CPU computations,
 * the 'H' stands for 'the host'
 * and the template parameter value T specifies the precision level.
 */

#include <iostream>
#include <fstream>

template <class T>
class complexH
{
	public:
		complexH() {};

		complexH operator=(const complexH&);
		complexH operator=(const int&);
		complexH operator=(const double&);
		complexH operator=(const dd_real&);
		complexH operator=(const qd_real&);

		complexH(const double&, const double&);
		complexH(const dd_real&, const dd_real&);
		complexH(const qd_real&, const qd_real&);

		complexH(const int&);
		complexH(const double&);

		complexH(const char*);
		complexH(const char*, const char*);
		complexH(const complexH&);

		/**
		 * Update value of complex number by (double, double)
		 */
		void init(const double&, const double&);
		/**
		 * Update value of complex number by (dd_real,dd_real)
		 */
		void init(const dd_real&, const dd_real&);
		/**
		 * Update value of complex number by (qd_real,qd_real)
		 */
		void init(const qd_real&, const qd_real&);

		/* Complex number computation operators*/
		complexH operator+(const complexH&);
		void     operator+=(const complexH&);
		complexH operator-(const complexH&);
		void     operator-=(const complexH&);
		complexH operator*(const complexH&) const;
		complexH operator*(const T&);
		complexH operator*(const int&);
		void     operator*=(const complexH&);
		complexH operator/(const complexH&);
		complexH operator/(const T&);

		/**
		 * Compute norm2 of complex number
		 */
		T absv() {return sqrt(real*real+imag*imag);};

		/**
		 * Compute conjugation of complex number
		 */
		complexH adj() {return complexH(real,-imag);};

		T real;
		T imag;

	/**
	 * Input stream operator for format a + b*i
	 */
	template <class T1>
	friend std::ostream& operator<<
	  ( std::ostream& os, const complexH<T1>& number );

	template <class T1>
	friend std::ifstream& operator>>
	  ( std::ifstream& is, const complexH<T1>& number );
};

#include "complexH.cpp"


#endif

/*template<class T, int size>
class storage
{
   public:

      complexH<T> array[size];

      void print()
      {
         for(int i=0; i<size; i++) cout << array[i];
      }

      void copyToStor(complexH<T>* a)
      {
         for(int i=0; i<size; i++) array[i]=a[i];
      }
};*/
