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

template <class T>
class mycomplex
{
   //template <class T1>
   //friend std::ostream& operator<<
   //   ( std::ostream& os, const mycomplex& number );
      // ( std::ostream& os, const mycomplex<T1>& number );

   // template <class T1>
   //friend std::ifstream& operator>>
   //   ( std::ifstream& is, const mycomplex& number );
      // ( std::ifstream& is, const mycomplex<T1>& number );

   public:

      mycomplex<T> operator=(mycomplex<T>);
      mycomplex(T,T,bool);
      mycomplex(double,double);
      mycomplex(const mycomplex<T> &);
      void init(double,double);

      mycomplex() {};
      mycomplex operator+(mycomplex);
      mycomplex operator-(mycomplex);
      mycomplex operator*(mycomplex);
      mycomplex operator*(T);
      mycomplex operator*(int);
      mycomplex operator/(mycomplex);
      T absv() {return sqrt(real*real+image*image);};
      mycomplex adj() {return mycomplex(real,-image,(bool)1);};
      // void init(double,double);
      //void init(qd_real,qd_real);

      T real;
      T image;
};

template<class T, int size>
class storage
{
   public:

      mycomplex<T> array[size];

      void print()
      {
         for(int i=0; i<size; i++) cout << array[i];
      }

      void copyToStor( mycomplex<T>* a)
      {
         for(int i=0; i<size; i++) array[i]=a[i];
      }

};

template <class T>
mycomplex<T>::mycomplex ( const mycomplex<T>& source )
{
   real = source.real;
   image = source.image;
}

template <class T>
std::ostream& operator<< ( std::ostream& os, const mycomplex<T>& number )
{
   int pr = 2*sizeof(T);
   cout.precision(pr);
   os << number.real << " + i*" << number.image << endl;
   return os;
}

template <class T>
std:: ifstream& operator >> ( std::ifstream& is, const mycomplex<T>& number )
{
   is >> number.real >> "+i*" >> number.image;
   return is;
}

template <class T>
mycomplex<T> mycomplex<T>::operator= ( mycomplex<T> a )
{
   real = a.real;
   image = a.image;

   return *this;
}

template <class T>
mycomplex<T>::mycomplex(T a, T b, bool c)
{
   real = a;
   image = b;
}

/*
template <class T>
mycomplex<T>::mycomplex(double a, double b)
{
real.x[0] = a;
image.x[0] =b;
}
*/

template <class T>
mycomplex<T>::mycomplex(double a, double b)
{
real.x[0] = a; real.x[1]=0.0;
image.x[0] =b; image.x[1]=0.0;
}

template <>
inline mycomplex<qd_real>::mycomplex(double a, double b)
{
real.x[0] = a; real.x[1]=0.0; real.x[2]=0.0; real.x[3]=0.0;
image.x[0] =b; image.x[1]=0.0; image.x[2]=0.0; image.x[3]=0.0;
}

template <>
inline mycomplex<double>::mycomplex ( double a, double b )
{
   real = a; 
   image = b; 
}

template <class T>
void mycomplex<T>::init(double a, double b)

{
   mycomplex<T> temp(a,b);

   real = temp.real;
   image = temp.image;   
}

//template <class T>
//mycomplex<T>::mycomplex() {}

template <class T>
mycomplex<T> mycomplex<T>::operator+(mycomplex<T> a)
{
   return mycomplex( real + a.real, image + a.image,1);
}

template <class T>
mycomplex<T> mycomplex<T>::operator-(mycomplex<T> a)
{
   return mycomplex( real - a.real, image - a.image,1);
}

template <class T>
mycomplex<T> mycomplex<T>::operator*(mycomplex<T> a)
{
   return mycomplex( real*a.real-image*a.image, image*a.real+real*a.image,1);
}

template <class T>
mycomplex<T> mycomplex<T>::operator*(T a)
{
  return mycomplex( real*a, image*a,1);
}

template <class T>
mycomplex<T> mycomplex<T>::operator*(int a)
{
   return mycomplex( real*a, image*a,1);
}

template <class T>
mycomplex<T> mycomplex<T>::operator/(mycomplex<T> a)
{

  return mycomplex( (real*a.real+image*a.image)/(a.real*a.real+a.image*a.image),
		  (image*a.real-real*a.image)/ (a.real*a.real+a.image*a.image),1);

}

#endif
