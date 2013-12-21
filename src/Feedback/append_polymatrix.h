/* append_dcmatrix.h contains the prototype for two functions that return
   the result of appending two matrices horizontally and vertically */

#include "poly_dcmplx.h"

void poly_h_append ( int n, int m1, int m2, POLY a[n][m1], POLY b[n][m2], POLY c[n][m1+m2]);
/* appends matrices a and b horizontally and saves the result in matrix c */

void poly_v_append ( int n1, int n2, int m, POLY a[n1][m], POLY b[n2][m], POLY c[n1+n2][m]);
/* appends matrices a and b vertically and saves the result in matrix c */ 
