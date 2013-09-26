/* file "c2ada_poly_dcmplx.h" contains the prototype for transfering a polynomial matrix from
   c to ada format */
#include "poly_dcmplx.h"

double* c2ada_poly_matrix( int n, int m, POLY c[n][m], int l, int a[l], int *sum );
/* given a polynomial matrix c, save the degree information in integer array a
   and return an double array of coefficients */
