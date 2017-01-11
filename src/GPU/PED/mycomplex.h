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

template <class real>
class mycomplex
{
   public:

      mycomplex<real> operator=(mycomplex<real>);
      mycomplex(real,real,bool);
      mycomplex(double,double);
      mycomplex(const mycomplex<real> &);
      void init(double,double);

      mycomplex() {};
      mycomplex operator+(mycomplex);
      mycomplex operator-(mycomplex);
      mycomplex operator*(mycomplex);
      mycomplex operator*(real);
      mycomplex operator*(int);
      mycomplex operator/(mycomplex);
      real absv() {return sqrt(re*re+im*im);};
      mycomplex adj() {return mycomplex(re,-im,(bool)1);};

      real re;
      real im;
};

template<class real, int size>
class storage
{
   public:

      mycomplex<real> array[size];

      void print()
      {
         for(int i=0; i<size; i++) cout << array[i];
      }

      void copyrealoStor( mycomplex<real>* a)
      {
         for(int i=0; i<size; i++) array[i]=a[i];
      }

};

template <class real>
mycomplex<real>::mycomplex ( const mycomplex<real>& source )
{
   re = source.re;
   im = source.im;
}

template <class real>
std::ostream& operator<< ( std::ostream& os, const mycomplex<real>& number )
{
   int pr = 2*sizeof(real);
   cout.precision(pr);
   os << number.re << " + i*" << number.im << endl;
   return os;
}

template <class real>
std:: ifstream& operator >>
 ( std::ifstream& is, const mycomplex<real>& number )
{
   is >> number.re >> "+i*" >> number.im;
   return is;
}

template <class real>
mycomplex<real> mycomplex<real>::operator= ( mycomplex<real> a )
{
   re = a.re;
   im = a.im;

   return *this;
}

template <class real>
mycomplex<real>::mycomplex(real a, real b, bool c)
{
   re = a;
   im = b;
}

template <class real>
mycomplex<real>::mycomplex(double a, double b)
{
   re.x[0] = a; re.x[1] = 0.0;
   im.x[0] = b; im.x[1] = 0.0;
}

template <>
inline mycomplex<qd_real>::mycomplex ( double a, double b )
{
   re.x[0] = a; re.x[1] = 0.0; re.x[2] = 0.0; re.x[3] = 0.0;
   im.x[0] = b; im.x[1] = 0.0; im.x[2] = 0.0; im.x[3] = 0.0;
}

template <>
inline mycomplex<double>::mycomplex ( double a, double b )
{
   re = a; 
   im = b; 
}

template <class real>
void mycomplex<real>::init ( double a, double b )
{
   mycomplex<real> temp(a,b);

   re = temp.re;
   im = temp.im;   
}

template <class real>
mycomplex<real> mycomplex<real>::operator+ ( mycomplex<real> a )
{
   return mycomplex(re + a.re,im + a.im,1);
}

template <class real>
mycomplex<real> mycomplex<real>::operator- ( mycomplex<real> a )
{
   return mycomplex(re - a.re, im - a.im,1);
}

template <class real>
mycomplex<real> mycomplex<real>::operator* ( mycomplex<real> a )
{
   return mycomplex(re*a.re-im*a.im, im*a.re+re*a.im,1);
}

template <class real>
mycomplex<real> mycomplex<real>::operator* ( real a )
{
   return mycomplex(re*a, im*a,1);
}

template <class real>
mycomplex<real> mycomplex<real>::operator* ( int a )
{
   return mycomplex(re*a, im*a,1);
}

template <class real>
mycomplex<real> mycomplex<real>::operator/ ( mycomplex<real> a )
{

  return mycomplex((re*a.re+im*a.im)/(a.re*a.re+a.im*a.im),
		   (im*a.re-re*a.im)/(a.re*a.re+a.im*a.im),1);
}

#endif
