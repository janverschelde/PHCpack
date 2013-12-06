#include "dcmplx.h"
#include "poly_dcmplx.h" 

int pivot_row (int n, int m, POLY a[n][m], int k, int z);
/* for all a[i][k], where i>=k and k is the index of nonzero columns,
   z is the number of zero columns, returns -1 if column j is zero, 
   otherwise return the first nonzero index. */

void Permute_row ( int n, int m, POLY a[n][m], POLY p[n][n], int r1, int r2);
/* interchange the rows r1 and r2 of matrix a and the corresponding
   left multiplier matrix p */

void Interchange_Rows1 ( int n, int m, POLY a[n][m], int r1, int r2);
/* Interchange the rows r1 and r2 of matrix a */

void Eliminate_Col1 ( int n, int m, POLY a[n][m], POLY p[n][n], int r, int c );
/* make zeros under a[r][c] */

void Hermite ( int n, int m, POLY a[n][m], POLY p[n][n]);
/* get the Hermite form of a n-by-m matrix */



















