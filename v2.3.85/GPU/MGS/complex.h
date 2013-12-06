/* Defining complex type for GPU computations.
Template parameter value specifies precision level */ 


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

//#include <qd/qd_config.h>
//#include <qd/dd_real.h>

using namespace std;

template <class T>
class complex{
/*
template <class T1>
  friend std::ostream& operator<< (std::ostream& os, const complex<T1>& number);

template <class T1>
  friend std::ifstream& operator>> (std::ifstream& is, con;

st complex<T1>& number);
*/
 public:

DEVICE complex<T> operator=(complex<T>);
//complex<T> operator=(complex<T>, bool b=0);
DEVICE complex(T,T,bool);
DEVICE complex();
DEVICE complex(double,double);
DEVICE void init(double,double);
void initH(double,double);

DEVICE complex operator+(complex);
DEVICE complex operator-(complex);
DEVICE complex operator*(complex);
DEVICE complex operator*(T);
DEVICE complex operator*(int);
DEVICE complex operator/(complex);
DEVICE T absv() {return sqrt(real*real+image*image);};
DEVICE complex adj() {return complex(real,0.0-image,(bool)1);};
 // void init(double,double);
 //void init(qd_real,qd_real);

 T real;
 T image;
};

/*
template <class T>
std::ostream& operator<< (std::ostream& os, const complex<T>& number)

{
  os << number.real << " + i*" << number.image << endl;
  return os;
}

template <class T>
std:: ifstream& operator >> (std::ifstream& is, const complex<T>& number)

{
  is >> number.real >> "+i*" >> number.image;
  return is;
}
*/


template <class T>
DEVICE complex<T> complex<T>::operator=(complex<T> a)
{
  real=a.real;
  image=a.image;

  return *this;
}

template <class T>
inline DEVICE complex<T>::complex(T a, T b, bool c){

  real=a;
  image=b;

}


template <class T>
DEVICE complex<T>::complex(double a, double b)
{
//if (sizeof(T)==8) { real.x=a; cout <<"privet"<< endl;}
//if (sizeof(T)==8) { T c(a); T d(b); real=c; image=d;}

//else{
real.x = a; real.y=0.0;
image.x =b; image.y=0.0;
//}
}

template <>
inline DEVICE complex<gqd_real>::complex(double a, double b)
{
//if (sizeof(T)==8) { real.x=a; cout <<"privet"<< endl;}
//if (sizeof(T)==8) { T c(a); T d(b); real=c; image=d;}

//else{
real.x = a; real.y=0.0; real.z=0.0; real.w=0.0;
image.x =b; image.y=0.0; image.z=0.0; image.w=0.0;
//}
}

template <>
inline DEVICE complex<double>::complex(double a, double b)
{
//if (sizeof(T)==8) { real.x=a; cout <<"privet"<< endl;}
//if (sizeof(T)==8) { T c(a); T d(b); real=c; image=d;}

//else{
real = a; 
image =b; 
//}
}

template <class T>
DEVICE void complex<T>::init(double a, double b)
{
  complex<T> temp(a,b);
  real=temp.real; image=temp.image;   
}

template <>
inline void complex<gqd_real>::initH(double a, double b)
{
  //complex<T> temp(a,b);
  //real=temp.real; image=temp.image;
  real.x=a; real.y=0.0; real.z=0; real.w=0;
  image.x=b; image.y=0.0; image.z=0.0; image.w=0; 
}

template <>
inline void complex<gdd_real>::initH(double a, double b)
{
  //complex<T> temp(a,b);
  //real=temp.real; image=temp.image;
  real.x=a; real.y=0.0; image.x=b; image.y=0.0;
}

template <>
inline void complex<double>::initH(double a, double b)
{
  //complex<T> temp(a,b);
  //real=temp.real; image=temp.image;
  real=a;  image=b;
}


template <class T>
void complex<T>::initH(double a, double b)

{
  //complex<T> temp(a,b);
  //real=temp.real; image=temp.image;
  //real.x=a; real.y=0.0; image.x=b; image.y=0.0;
}

template <class T>
DEVICE complex<T>::complex() {}

template <class T>
DEVICE complex<T> complex<T>::operator+(complex<T> a)
{
  return complex( real + a.real, image + a.image,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator-(complex<T> a)
{
  return complex( real - a.real, image - a.image,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator*(complex<T> a)
{
  return complex( real*a.real-image*a.image, image*a.real+real*a.image,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator*(T a)
{
  return complex( real*a, image*a,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator*(int a)
{
  return complex( real*a, image*a,1);
}

template <class T>
DEVICE complex<T> complex<T>::operator/(complex<T> a)
{
  return complex( (real*a.real+image*a.image)/(a.real*a.real+a.image*a.image),
		  (image*a.real-real*a.image)/ (a.real*a.real+a.image*a.image),1);

}

#endif
