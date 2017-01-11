/* Defining complex type for CPU computations (the host).
   Template parameter value specifies precision level */

#ifndef __COMPLEXH_H__
#define __COMPLEXH_H__

#define __CUDAC__

#ifdef __CUDAC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

#include <iostream>
#include <qd/qd_real.h>

using namespace std;

template <class realH>
class complexH
{
   template <class realD>
   friend std::ostream& operator<<
      ( std::ostream& os, const complexH<realD>& number );

   template <class realD>
   friend std::ifstream& operator>>
      ( std::ifstream& is, const complexH<realD>& number );

   public:

      complexH<realH> operator=(complexH<realH>);
      complexH(realH,realH,bool);
      complexH(double,double);
      complexH(const complexH<realH> &);
      void init(double,double);

      complexH() {};
      complexH operator+(complexH);
      complexH operator-(complexH);
      complexH operator*(complexH);
      complexH operator*(realH);
      complexH operator*(int);
      complexH operator/(complexH);
      realH absv() {return sqrt(re*re+im*im);};
      complexH adj() {return complexH(re,-im,(bool)1);};

      realH re;
      realH im;
};

template<class realH, int size>
class storage
{
   public:

      complexH<realH> array[size];

      void print()
      {
         for(int i=0; i<size; i++) cout << array[i];
      }

      void copyToStor( complexH<realH>* a)
      {
         for(int i=0; i<size; i++) array[i]=a[i];
      }

};

template <class realH>
complexH<realH>::complexH ( const complexH<realH>& source )
{
   re = source.re;
   im = source.im;
}

template <class realH>
std::ostream& operator<< ( std::ostream& os, const complexH<realH>& number )
{
   int pr = 2*sizeof(realH);
   cout.precision(pr);
   os << number.re << " + i*" << number.im << endl;
   return os;
}

template <class realH>
std:: ifstream& operator >> ( std::ifstream& is, const complexH<realH>& number )
{
   is >> number.re >> "+i*" >> number.im;
   return is;
}

template <class realH>
complexH<realH> complexH<realH>::operator= ( complexH<realH> a )
{
   re = a.re;
   im = a.im;

   return *this;
}

template <class realH>
complexH<realH>::complexH(realH a, realH b, bool c)
{
   re = a;
   im = b;
}

template <class realH>
complexH<realH>::complexH(double a, double b)
{
   re.x[0] = a; re.x[1] = 0.0;
   im.x[0] = b; im.x[1] = 0.0;
}

template <>
inline complexH<qd_real>::complexH(double a, double b)
{
   re.x[0] = a; re.x[1] = 0.0; re.x[2] = 0.0; re.x[3] = 0.0;
   im.x[0] = b; im.x[1] = 0.0; im.x[2] = 0.0; im.x[3] = 0.0;
}

template <>
inline complexH<double>::complexH ( double a, double b )
{
   re = a; 
   im = b; 
}

template <class realH>
void complexH<realH>::init(double a, double b)
{
   complexH<realH> temp(a,b);

   re = temp.re;
   im = temp.im;
}

template <class realH>
__device__ complexH<realH> complexH<realH>::operator+(complexH<realH> a)
{
   return complexH(re + a.re, im + a.im,1);
}

template <class realH>
complexH<realH> complexH<realH>::operator-(complexH<realH> a)
{
   return complexH(re - a.re, im - a.im,1);
}

template <class realH>
complexH<realH> complexH<realH>::operator*(complexH<realH> a)
{
   return complexH(re*a.re-im*a.im, im*a.re+re*a.im,1);
}

template <class realH>
complexH<realH> complexH<realH>::operator*(realH a)
{
  return complexH(re*a, im*a,1);
}

template <class realH>
complexH<realH> complexH<realH>::operator*(int a)
{
   return complexH(re*a, im*a,1);
}

template <class realH>
complexH<realH> complexH<realH>::operator/(complexH<realH> a)
{
   return complexH((re*a.re+im*a.im)/(a.re*a.re+a.im*a.im),
                   (im*a.re-re*a.im)/(a.re*a.re+a.im*a.im),1);
}

#endif
