#include "dcmplx.h"
#include "poly_dcmplx.h" 

void Find_pivots (int n, int m, POLY a[n][m], int *pivrow, int *pivcol);
/* find the first nonzero entry in the matrix a, starting at current pivot */

void Permute ( int n, int m, POLY a[n][m], POLY p[n][n], POLY q[m][m],
			   int pivrow, int row, int pivcol, int col);
/*Switches the pivot rows and columns in a to the current row and column.
  The unimodular matrices p and q are pivoted along.*/

void Interchange_Rows ( int n, int m, POLY a[n][m], int r1, int r2);
/* Interchange the rows r1 and r2 of matrix a */

void Interchange_Cols ( int n, int m, POLY a[n][m], int c1, int c2);
/* Interchange the columns c1 and c2 of matrix a */

void Eliminate_Col ( int n, int m, POLY a[n][m], POLY p[n][n], int r, int c );
/* make zeros under a[r][c] */

void Eliminate_Row ( int n, int m, POLY a[n][m], POLY q[m][m], int r, int c );
/* make all elemnets after a[r][c] zeros in the row r */


void Smith_Diagonal ( int n, int m, POLY a[n][m], POLY p[n][n], POLY q[m][m]);
/* use extended gcd method to get the diagonal form of a n-by-m matrix, which
   is a subroutine of Smith function */

int Diagonal ( int n, int m , POLY a[n][m]);
/* Returns 1 if the matrix is diagonal, returns 0 otherwise. */

int poly_divide ( POLY a, POLY b );
/* Returns 1 if a divides into b without remainder, returns 0 otherwise. */ 

void Smith ( int n, int m, POLY a[n][m], POLY p[n][n], POLY q[m][m]);
/* get the Smith form of a n-by-m matrix */


















