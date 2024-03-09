// The file dbl16_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl16_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "hexa_double_gpufun.cu"
#include "dbl16_sqrt_kernels.h"

__global__ void dbl16_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int max_steps )
{
}

void GPU_dbl16_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int maxstp )
{
}
