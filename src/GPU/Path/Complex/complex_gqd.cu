/* complex.cu created by yxc on Feb 1, 2015, with edits by jv */

#ifndef COMPLEX_GQD_CU_
#define COMPLEX_GQD_CU_

#include <gqd_basic.cu>
#include <gqd_sqrt.cu>

#include "complex_gqd.h"

__device__ gqd_complex::gqd_complex ( double a, double b )
{
   real.x = a; real.y = 0.0; real.z = 0.0; real.w = 0.0;
   imag.x = b; imag.y = 0.0; imag.z = 0.0; imag.w = 0.0;
}

__device__ gqd_complex::gqd_complex ( gqd_real a, gqd_real b )
{
   real = a;
   imag = b;
}

__device__ gqd_complex gqd_complex::operator+ ( gqd_complex a )
{
   return gqd_complex(real+a.real,imag+a.imag);
}

__device__ gqd_complex gqd_complex::operator* ( gqd_complex a )
{
   return gqd_complex(real*a.real-imag*a.imag,imag*a.real+real*a.imag);
}

__device__ gqd_complex gqd_complex::operator- ( gqd_complex a )
{
   return gqd_complex(real-a.real,imag-a.imag);
}

__device__ gqd_complex gqd_complex::operator/ ( gqd_complex a )
{
   return gqd_complex((real*a.real+imag*a.imag)/(a.real*a.real+a.imag*a.imag),
                      (imag*a.real-real*a.imag)/(a.real*a.real+a.imag*a.imag));
}

__device__ void gqd_complex::operator*= ( gqd_complex a )
{
   gqd_real real_tmp = real;
   real = real*a.real - imag*a.imag;
   imag = imag*a.real + real_tmp*a.imag;
}

__device__ void gqd_complex::operator+= ( gqd_complex a )
{
   real = real + a.real;
   imag = imag + a.imag;
}

__device__ void gqd_complex::operator/= ( gqd_real a )
{
   real = real/a;
   imag = imag/a;
}

__device__ double gqd_complex::norm_double()
{
   return real.x*real.x + imag.x*imag.x;
}

__device__ double gqd_complex::norm1_double()
{
   return abs(real.x) + abs(imag.x);
}

__device__ void gqd_complex::init_imag()
{
   imag.x = 0.0; imag.y = 0.0; imag.z = 0.0; imag.w = 0.0;
}

__device__ gqd_complex gqd_complex::adj()
{
   return gqd_complex(real,0.0-imag);
}

__device__ gqd_complex gqd_complex::adj_multiple ( gqd_complex a )
{
   return gqd_complex(real*a.real+imag*a.imag,real*a.imag-imag*a.real);
}

#endif /* COMPLEX_GQD_CU_ */
