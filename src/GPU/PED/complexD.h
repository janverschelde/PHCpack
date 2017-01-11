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

template <class real>
class complexD
{
   public:

      DEVICE complexD<real> operator=(complexD<real>);
      DEVICE complexD(real,real,bool);
      DEVICE complexD();
      DEVICE complexD(double,double);
      DEVICE void init(double,double);
      void initH(double,double);

      DEVICE complexD operator+(complexD);
      DEVICE complexD operator-(complexD);
      DEVICE complexD operator*(complexD);
      DEVICE complexD operator*(real);
      DEVICE complexD operator*(int);
      DEVICE complexD operator/(complexD);
      DEVICE real absv() {return sqrt(re*re+im*im);};
      DEVICE complexD adj() {return complexD(re,0.0-im,(bool)1);};

      real re;
      real im;
};

template <class real>
DEVICE complexD<real> complexD<real>::operator=(complexD<real> a)
{
   re=a.re;
   im=a.im;

   return *this;
}

template <class real>
inline DEVICE complexD<real>::complexD(real a, real b, bool c)
{
   re = a;
   im = b;
}

template <class real>
DEVICE complexD<real>::complexD(double a, double b)
{
   re.x = a; re.y = 0.0;
   im.x = b; im.y = 0.0;
}

template <>
inline DEVICE complexD<gqd_real>::complexD(double a, double b)
{
   re.x = a; re.y = 0.0; re.z = 0.0; re.w = 0.0;
   im.x = b; im.y = 0.0; im.z = 0.0; im.w = 0.0;
}

template <>
inline DEVICE complexD<double>::complexD(double a, double b)
{
   re = a; 
   im = b; 
}

template <class real>
DEVICE void complexD<real>::init(double a, double b)
{
   complexD<real> temp(a,b);
   re = temp.re;
   im = temp.im;   
}

template <>
inline void complexD<gqd_real>::initH(double a, double b)
{
   re.x = a; re.y = 0.0; re.z = 0.0; re.w = 0.0;
   im.x = b; im.y = 0.0; im.z = 0.0; im.w = 0.0;
}

template <>
inline void complexD<gdd_real>::initH(double a, double b)
{
   re.x=a; re.y=0.0; im.x=b; im.y=0.0;
}

template <>
inline void complexD<double>::initH(double a, double b)
{
   re=a;
   im=b;
}

template <class real>
void complexD<real>::initH(double a, double b)
{
  //complex<real> temp(a,b);
  //re=temp.re; ime=temp.ime;
  //re.x=a; re.y=0.0; ime.x=b; ime.y=0.0;
}

template <class real>
DEVICE complexD<real>::complexD() {}

template <class real>
DEVICE complexD<real> complexD<real>::operator+(complexD<real> a)
{
   return complexD(re + a.re, im + a.im,1);
}

template <class real>
DEVICE complexD<real> complexD<real>::operator-(complexD<real> a)
{
   return complexD(re-a.re,im-a.im,1);
}

template <class real>
DEVICE complexD<real> complexD<real>::operator*(complexD<real> a)
{
   return complexD(re*a.re-im*a.im, im*a.re+re*a.im,1);
}

template <class real>
DEVICE complexD<real> complexD<real>::operator*(real a)
{
   return complexD(re*a,im*a,1);
}

template <class real>
DEVICE complexD<real> complexD<real>::operator*(int a)
{
   return complexD(re*a,im*a,1);
}

template <class real>
DEVICE complexD<real> complexD<real>::operator/(complexD<real> a)
{
   return complexD((re*a.re+im*a.im)/(a.re*a.re+a.im*a.im),
		   (im*a.re-re*a.im)/(a.re*a.re+a.im*a.im),1);
}

#endif
