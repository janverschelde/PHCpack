/* file "c2ada_dc_dcmplx.h" contains the prototype for transfering a polynomial matrix from
   c to ada format */

#include "dcmplx.h"

void c2ada_dc_matrix( int n, int m, dcmplx c[n][m], int db, double b[db], int start);
