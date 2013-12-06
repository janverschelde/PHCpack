/* file "dc_inverse.h" contains prototypes for computing the inverse
   matrix of an n-by-n matrix of complex numbers based on LU-decomposition*/
#include "dcmplx.h"


void dcinverse ( int n, dcmplx a[n][n] );
/* computes the inverse matrix of an n-by-n matrix a,
   the matrix a on return contains the result */
