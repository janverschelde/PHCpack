/* complex_gd.h created by yxc on Feb 2, 2015, edits by jv */

#ifndef COMPLEX_GD_H_
#define COMPLEX_GD_H_

#include "gqd_type.h"

class gd_complex
{
   public:

      __device__ gd_complex(){};
      __device__ gd_complex(double,double);
      __device__ gd_complex operator+(gd_complex);
      __device__ gd_complex operator*(gd_complex);
      __device__ gd_complex operator-(gd_complex);
      __device__ gd_complex operator/(gd_complex);
      __device__ void operator*=(gd_complex);
      __device__ void operator+=(gd_complex);
      __device__ void operator/=(double);
      __device__ gd_complex adj();
      __device__ gd_complex adj_multiple(gd_complex a);
      __device__ double norm_double();
      __device__ double norm1_double();
      __device__ void init_imag();

     double real;
     double imag;
};

#endif /* COMPLEX_GD_H_ */
