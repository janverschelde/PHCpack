/* complex.cu created by yxc on Feb 1, 2015 with edits by jv */

#ifndef COMPLEX_GD_CU_
#define COMPLEX_GD_CU_

#include "complex.h"

__device__ gd_complex::gd_complex ( double a, double b )
{
   real = a;
   imag = b;
}

__device__ gd_complex gd_complex::operator+ ( gd_complex a )
{
   return gd_complex(real+a.real,imag+a.imag);
}

__device__ gd_complex gd_complex::operator* ( gd_complex a )
{
   return gd_complex(real*a.real-imag*a.imag,imag*a.real+real*a.imag);
}

__device__ gd_complex gd_complex::operator- ( gd_complex a )
{
   return gd_complex(real-a.real,imag-a.imag);
}

__device__ gd_complex gd_complex::operator/( gd_complex a )
{
   return gd_complex((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
                     (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag));
}

__device__ void gd_complex::operator*= ( gd_complex a )
{
   double real_tmp = real;
   real = real*a.real - imag*a.imag;
   imag = imag*a.real + real_tmp*a.imag;
}

__device__ void gd_complex::operator+=( gd_complex a )
{
   real = real + a.real;
   imag = imag + a.imag;
}

__device__ void gd_complex::operator/= ( double a )
{
   real = real/a;
   imag = imag/a;
}

__device__ double gd_complex::norm_double()
{
   return real + imag*imag;
}

__device__ double gd_complex::norm1_double()
{
   return abs(real) + abs(imag);
}

__device__ void gd_complex::init_imag()
{
   imag = 0.0;
}

__device__ gd_complex gd_complex::adj()
{
   return gd_complex(real,0.0-imag);
}

__device__ gd_complex gd_complex::adj_multiple ( gd_complex a )
{
   return gd_complex(real*a.real+imag*a.imag,real*a.imag-imag*a.real);
}

#endif /* COMPLEX_GD_CU_ */
