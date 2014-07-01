/* The class complex defines a complex type for GPU computations.
   The template parameter value T specifies precision level */ 

#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#define __CUDAC__

#ifdef __CUDAC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

#include <iostream>
#include <gqd_type.h>

using namespace std;

template <class T>
class complex
{
   public:

      DEVICE complex<T> operator=(complex<T>);
      DEVICE complex(T,T,bool);
      DEVICE complex();
      DEVICE complex(double,double);
      DEVICE complex(double,double,double,double);
      DEVICE complex(double,double,double,double,
                     double,double,double,double);
      DEVICE void init(double,double);
      DEVICE void init(double,double,double,double);
      DEVICE void init(double,double,double,double,
                       double,double,double,double);
      void initH(double,double);
      void initH(double,double,double,double);
      void initH(double,double,double,double,
                 double,double,double,double);
 
      DEVICE complex operator+(complex);
      DEVICE complex operator-(complex);
      DEVICE complex operator*(complex);
      DEVICE complex operator*(T);
      DEVICE complex operator*(int);
      DEVICE complex operator/(complex);
      DEVICE complex operator/(T);
      DEVICE T absv() {return sqrt(real*real+imag*imag);};
      DEVICE complex adj() {return complex(real,0.0-imag,(bool)1);};
 
      T real;
      T imag;
};

template <class T>
DEVICE complex<T> complex<T>::operator=(complex<T> a)
{
   real=a.real;
   imag=a.imag;

   return *this;
}

template <class T>
inline DEVICE complex<T>::complex(T a, T b, bool c)
{
   real = a;
   imag = b;
}

template <class T>
DEVICE complex<T>::complex(double a, double b)
{
   real.x = a; real.y = 0.0;
   imag.x = b; imag.y = 0.0;
}

template <>
inline DEVICE complex<gdd_real>::complex(double a, double b)
{
   real.x = a; real.y = 0.0;
   imag.x = b; imag.y = 0.0;
}

template <>
inline DEVICE complex<gqd_real>::complex(double a, double b)
{
   real.x = a; real.y = 0.0; real.z = 0.0; real.w = 0.0;
   imag.x = b; imag.y = 0.0; imag.z = 0.0; imag.w = 0.0;
}

template <>
inline DEVICE complex<double>::complex(double a, double b)
{
   real = a; 
   imag = b; 
}

template <class T>
DEVICE void complex<T>::init(double a, double b)
{
   complex<T> temp(a,b);
   real = temp.real; imag = temp.imag;   
}

template <class T>
DEVICE void complex<T>::init(double ahi, double alo,
                             double bhi, double blo)
{
   complex<T> temp(ahi,alo,bhi,blo);
   real = temp.real; imag = temp.imag;   
}

template <class T>
DEVICE void complex<T>::init(double ahihi, double ahilo,
                             double alohi, double alolo,
                             double bhihi, double bhilo,
                             double blohi, double blolo)
{
   complex<T> temp(ahihi,ahilo,alohi,alolo,bhihi,bhilo,blohi,blolo);
   real = temp.real; imag = temp.imag;   
}

template <>
inline void complex<gqd_real>::initH(double a, double b)
{
   real.x = a; real.y = 0.0; real.z = 0.0; real.w = 0.0;
   imag.x = b; imag.y = 0.0; imag.z = 0.0; imag.w = 0.0; 
}

template <>
inline void complex<gqd_real>::initH(double ahihi, double alohi,
                                     double ahilo, double alolo,
                                     double bhihi, double blohi,
                                     double bhilo, double blolo)
{
   real.x = ahihi; real.y = alohi; real.z = ahilo; real.w = alolo;
   imag.x = bhihi; imag.y = blohi; imag.z = bhilo; imag.w = blolo; 
}

template <>
inline void complex<gdd_real>::initH(double a, double b)
{
   real.x = a; real.y = 0.0; imag.x = b; imag.y = 0.0;
}

template <>
inline void complex<gdd_real>::initH(double ahi, double alo,
                                     double bhi, double blo)
{
   real.x = ahi; real.y = alo; imag.x = bhi; imag.y = blo;
}

template <>
inline void complex<double>::initH(double a, double b)
{
   real = a;  imag = b;
}

template <class T>
void complex<T>::initH(double a, double b)
{
  //complex<T> temp(a,b);
  //real=temp.real; imag=temp.imag;
  //real.x=a; real.y=0.0; imag.x=b; imag.y=0.0;
}

template <class T>
DEVICE complex<T>::complex() {}

template <class T>
DEVICE complex<T> complex<T>::operator+(complex<T> a)
{
   return complex(real+a.real,imag+a.imag,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator-(complex<T> a)
{
   return complex(real-a.real,imag-a.imag,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator*(complex<T> a)
{
   return complex(real*a.real-imag*a.imag,imag*a.real+real*a.imag,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator*(T a)
{
   return complex(real*a,imag*a,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator*(int a)
{
   return complex(real*a,imag*a,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator/(T a)
{
   return complex(real/a,imag/a,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator/(complex<T> a)
{
   return complex((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
                  (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag),1);
}
#endif
