/* complex_gqd.h created by yxc on Feb 2, 2015, with edits by jv */

#ifndef COMPLEX_GQD_H_
#define COMPLEX_GQD_H_

#include "gqd_type.h"

class gqd_complex
{
   public:

      __device__ gqd_complex(){};
      __device__ gqd_complex(double,double);
      __device__ gqd_complex(gqd_real,gqd_real);
      __device__ gqd_complex operator+(gqd_complex);
      __device__ gqd_complex operator*(gqd_complex);
      __device__ gqd_complex operator-(gqd_complex);
      __device__ gqd_complex operator/(gqd_complex);
      __device__ void operator*=(gqd_complex);
      __device__ void operator+=(gqd_complex);
      __device__ void operator/=(gqd_real);
      __device__ gqd_complex adj();
      __device__ gqd_complex adj_multiple(gqd_complex a);
      __device__ double norm_double();
      __device__ double norm1_double();
      __device__ void init_imag();

      gqd_real real;
      gqd_real imag;
};

#endif /* COMPLEX_GQD_H_ */
