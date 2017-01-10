/* Defining complex type for GPU computations.
   Template parameter value specifies precision level */ 

#ifndef __COMPLEXD_H__
#define __COMPLEXD_H__

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
class complexD
{
   public:

      DEVICE complexD<T> operator=(complexD<T>);
      DEVICE complexD(T,T,bool);
      DEVICE complexD();
      DEVICE complexD(double,double);
      DEVICE void init(double,double);
      void initH(double,double);

      DEVICE complexD operator+(complexD);
      DEVICE complexD operator-(complexD);
      DEVICE complexD operator*(complexD);
      DEVICE complexD operator*(T);
      DEVICE complexD operator*(int);
      DEVICE complexD operator/(complexD);
      DEVICE T absv() {return sqrt(real*real+imag*imag);};
      DEVICE complexD adj() {return complexD(real,0.0-imag,(bool)1);};

      T real;
      T imag;
};

template <class T>
DEVICE complexD<T> complexD<T>::operator=(complexD<T> a)
{
   real=a.real;
   imag=a.imag;

   return *this;
}

template <class T>
inline DEVICE complexD<T>::complexD(T a, T b, bool c)
{
   real = a;
   imag = b;
}

template <class T>
DEVICE complexD<T>::complexD(double a, double b)
{
   real.x = a; real.y = 0.0;
   imag.x = b; imag.y = 0.0;
}

template <>
inline DEVICE complexD<gqd_real>::complexD(double a, double b)
{
   real.x = a; real.y = 0.0; real.z = 0.0; real.w = 0.0;
   imag.x = b; imag.y = 0.0; imag.z = 0.0; imag.w = 0.0;
}

template <>
inline DEVICE complexD<double>::complexD(double a, double b)
{
   real = a; 
   imag = b; 
}

template <class T>
DEVICE void complexD<T>::init(double a, double b)
{
   complexD<T> temp(a,b);
   real = temp.real;
   imag = temp.imag;   
}

template <>
inline void complexD<gqd_real>::initH(double a, double b)
{
   real.x = a; real.y = 0.0; real.z = 0.0; real.w = 0.0;
   imag.x = b; imag.y = 0.0; imag.z = 0.0; imag.w = 0.0;
}

template <>
inline void complexD<gdd_real>::initH(double a, double b)
{
   real.x=a; real.y=0.0; imag.x=b; imag.y=0.0;
}

template <>
inline void complexD<double>::initH(double a, double b)
{
   real=a;
   imag=b;
}

template <class T>
void complexD<T>::initH(double a, double b)
{
  //complex<T> temp(a,b);
  //real=temp.real; image=temp.image;
  //real.x=a; real.y=0.0; image.x=b; image.y=0.0;
}

template <class T>
DEVICE complexD<T>::complexD() {}

template <class T>
DEVICE complexD<T> complexD<T>::operator+(complexD<T> a)
{
   return complexD(real + a.real, imag + a.imag,1);
}

template <class T>
DEVICE complexD<T> complexD<T>::operator-(complexD<T> a)
{
   return complexD(real-a.real,imag-a.imag,1);
}

template <class T>
DEVICE complexD<T> complexD<T>::operator*(complexD<T> a)
{
   return complexD(real*a.real-imag*a.imag, imag*a.real+real*a.imag,1);
}

template <class T>
DEVICE complexD<T> complexD<T>::operator*(T a)
{
   return complexD(real*a,imag*a,1);
}

template <class T>
DEVICE complexD<T> complexD<T>::operator*(int a)
{
   return complexD(real*a,imag*a,1);
}

template <class T>
DEVICE complexD<T> complexD<T>::operator/(complexD<T> a)
{
   return complexD((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
		   (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag),1);
}

#endif
