/* complexD.cu defines the functions with prototypes in complexD.h */

#ifndef __COMPLEXD_CU__
#define __COMPLEXD_CU__

#include "complexD.h"

template <class T>
__device__ complexD<T>::complexD<T> ( T a, T b )
{
   real = a;
   imag = b;
}

template <class T>
__device__ complexD<T> complexD<T>::operator+ ( complexD<T> a )
{
   return complexD<T>(real+a.real,imag+a.imag);
}

template <class T>
__device__ complexD<T> complexD<T>::operator* ( complexD<T> a )
{
   return complexD<T>(real*a.real-imag*a.imag,imag*a.real+real*a.imag);
}

template <class T>
__device__ complexD<T> complexD<T>::operator- ( complexD<T> a )
{
   return complexD<T>(real-a.real,imag-a.imag);
}

template <class T>
__device__ complexD<T> complexD<T>::operator/( complexD<T> a )
{
   return complexD<T>((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
                      (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag));
}

template <class T>
__device__ void complexD<T>::operator*= ( complexD<T> a )
{
   double real_tmp = real;
   real = real*a.real - imag*a.imag;
   imag = imag*a.real + real_tmp*a.imag;
}

template <class T>
__device__ void complexD<T>::operator+=( complexD<T> a )
{
   real = real + a.real;
   imag = imag + a.imag;
}

template <class T>
__device__ void complexD<T>::operator/= ( double a )
{
   real = real/a;
   imag = imag/a;
}

template <class T>
__device__ double complexD<T>::norm_double()
{
   return real*real + imag*imag;
}

template <class T>
__device__ double complexD<T>::norm1_double()
{
   return abs(real) + abs(imag);
}

template <class T>
__device__ void complexD<T>::init_imag()
{
   imag = 0.0;
}

template <class T>
__device__ complexD<T> complexD<T>::adj()
{
   return complexD<T>(real,0.0-imag);
}

template <class T>
__device__ complexD<T> complexD<T>::adj_multiple ( complexD<T> a )
{
   return complexD<T>(real*a.real+imag*a.imag,real*a.imag-imag*a.real);
}

#endif /* __COMPLEXD_CU__ */
