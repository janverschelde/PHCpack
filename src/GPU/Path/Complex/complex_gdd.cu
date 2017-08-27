/* complex.cu created by yxc on Feb 1, 2015, with edits by jv */

#ifndef COMPLEX_GDD_CU_
#define COMPLEX_GDD_CU_

#include <gdd_basic.cu>
#include <gdd_sqrt.cu>

#include "complex_gdd.h"

__device__ gdd_complex::gdd_complex ( double a, double b )
{
   real.x = a; real.y = 0.0;
   imag.x = b; imag.y = 0.0;
}

__device__ gdd_complex::gdd_complex ( gdd_real a, gdd_real b )
{
   real = a;
   imag = b;
}

__device__ gdd_complex gdd_complex::operator+ ( gdd_complex a )
{
   return gdd_complex(real+a.real,imag+a.imag);
}

__device__ gdd_complex gdd_complex::operator* ( gdd_complex a )
{
   return gdd_complex(real*a.real-imag*a.imag,imag*a.real+real*a.imag);
}

__device__ gdd_complex gdd_complex::operator- ( gdd_complex a )
{
   return gdd_complex(real-a.real,imag-a.imag);
}

__device__ gdd_complex gdd_complex::operator/ ( gdd_complex a )
{
   return gdd_complex((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
                      (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag));
}

__device__ void gdd_complex::operator*= ( gdd_complex a )
{
   gdd_real real_tmp = real;
   real = real*a.real - imag*a.imag;
   imag = imag*a.real + real_tmp*a.imag;
}

__device__ void gdd_complex::operator+= ( gdd_complex a )
{
   real = real + a.real;
   imag = imag + a.imag;
}

__device__ void gdd_complex::operator/= ( gdd_real a )
{
   real = real/a;
   imag = imag/a;
}

__device__ double gdd_complex::norm_double()
{
   return real.x*real.x + imag.x*imag.x;
}

__device__ double gdd_complex::norm1_double()
{
   return abs(real.x) + abs(imag.x);
}

__device__ void gdd_complex::init_imag()
{
   imag.x = 0.0; imag.y = 0.0;
}

__device__ gdd_complex gdd_complex::adj()
{
   return gdd_complex(real,0.0-imag);
}

__device__ gdd_complex gdd_complex::adj_multiple ( gdd_complex a )
{
   return gdd_complex(real*a.real+imag*a.imag,real*a.imag-imag*a.real);
}

#endif /* COMPLEX_GDD_CU_ */
