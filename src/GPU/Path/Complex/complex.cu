/* complex.cu created by yxc on Feb 1, 2015, with edits by jv */

#ifndef COMPLEX_CU_
#define COMPLEX_CU_

#include "complex.h"

#if(path_precision == 0)
#include "complex_gd.cu"
#elif(path_precision == 1)
#include "complex_gdd.cu"
#elif(path_precision == 2)
#include "complex_gqd.cu"
#endif

#endif /* COMPLEX_CU_ */
