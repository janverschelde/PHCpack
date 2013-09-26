/* append_dcmatrix.h contains the prototype for two functions that return
   the result of appending two matrices horizontally and vertically */

#include "dcmplx.h"

void h_append ( int n, int m1, int m2, dcmplx a[n][m1], dcmplx b[n][m2], dcmplx c[n][m1+m2]);
/* appends matrices a and b horizontally and saves the result in matrix c */

void v_append ( int n1, int n2, int m, dcmplx a[n1][m], dcmplx b[n2][m], dcmplx c[n1+n2][m]);
/* appends matrices a and b vertically and saves the result in matrix c */ 
