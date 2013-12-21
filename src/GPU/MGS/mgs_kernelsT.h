#ifndef __KERNELS__
#define __KERNELS__

#include <iostream>
#include <cmath>
#include <gqd_type.h>
#include "DefineType.h"
#include "complex.h"
#include "vector_types.h"
#include "vector_functions.h"

void GPU_GS 
 ( complex<T>* v_h, complex<T>* R_h, complex<T>* sol_h, 
   int dim, int dimR, int k, int r, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the linear
 *   system into the allocated memory, launched the modified Gram-Schmidt
 *   method, the back subsitution, and transfers the solution of the system
 *   from global memory of the card to the memory of the host.
 *
 * ON ENTRY :
 *   v_h
 *   R_h
 *   sol_h
 *   dim
 *   dim_R
 *   k
 *   r
 *   BS      block size, number of threads per block.
 *
 * ON RETURN :
 *   v_h
 *   R_h
 *   sol_h                                                            */

#endif
