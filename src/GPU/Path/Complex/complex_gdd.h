/* complex_gdd.h created by yxc on Feb 1, 2015, with edits by jv */

#ifndef COMPLEX_GDD_H_
#define COMPLEX_GDD_H_

#include "gqd_type.h"

class gdd_complex
{
   public:

      __device__ gdd_complex(){};
      __device__ gdd_complex(double,double);
      __device__ gdd_complex(gdd_real,gdd_real);
      __device__ gdd_complex operator+(gdd_complex);
      __device__ gdd_complex operator*(gdd_complex);
      __device__ gdd_complex operator-(gdd_complex);
      __device__ gdd_complex operator/(gdd_complex);
      __device__ void operator*=(gdd_complex);
      __device__ void operator+=(gdd_complex);
      __device__ void operator/=(gdd_real);
      __device__ gdd_complex adj();
      __device__ gdd_complex adj_multiple(gdd_complex a);
      __device__ double norm_double();
      __device__ double norm1_double();
      __device__ void init_imag();

      gdd_real real;
      gdd_real imag;
};

#endif /* COMPLEX_GDD_H_ */
