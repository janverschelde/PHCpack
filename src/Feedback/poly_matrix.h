/* file poly_matrix.h contains the prototype for the functions 
 * of polynomial matrices operations*/
#include "poly_dcmplx.h"

void read_matrix ( int n, int m, POLY a[n][m] );
/* asks the user to provide the entries of an n-by-m POLY matrix */

void read_matrix1 ( int n, int m, POLY a[n][m] );
/* asks the user to provide the entries of an n-by-m POLY matrix */

void random_matrix ( int n, int m, POLY a[n][m] );
/* generate a random n-by-m POLY matrix */

void I_matrix ( int n, POLY POLY[n][n] );
/* assign matrix p as an identity matrix */

void zero_matrix ( int n, int m, POLY p[n][m]);
/* assing a n by m zero matrix */

void copy ( int n, int m, POLY a[n][m], POLY t[n][m]);
/* copy the original matrix to a temporary matrix t */

void print ( int n, int m, POLY a[n][m] );
/* writes n-by-m Polynomial matrix to the screen */

void print1 ( int n, int m, POLY a[n][m] );
/* writes n-by-m Polynomial matrix to the screen */

void Multiply ( int n, int m, int l, POLY a[n][m], POLY b[m][l], POLY c[n][l] );
/* multiply matrices a and b, put the result in c matrix */

void Transpose ( int n, int m, POLY a[n][m], POLY b[m][n] );
/* get the transpose of matrix a and put it into matrix b*/

void free_matrix ( int n, int m, POLY a[n][m]);
/* free the memory for each entry in a polynomial matrix */

POLY Inverse_Poly ( int n, POLY (*M)[n] );
/* get the inverse matrix of polynomial matrix M, M_Inverse = 1/ds * M, where the ds is the
least 
  common denominator of the inverse matrix, M is the corresponding polynomial matrix */

void dcmatrix_Multiply ( int n, int m, int l, dcmplx a[n][m], POLY b[m][l], POLY c[n][l] );
/* multiply complex matrix a and polynomial matrix b, save the result in polynomial matrix c */

void add_polymatrix( int n, int m, POLY a[n][m], POLY b[n][m], POLY c[n][m] );
/* add matrices a and b, save the result in c matrix */

void sub_polymatrix( int n, int m, POLY a[n][m], POLY b[n][m], POLY c[n][m] );
/* matrix a subtracts matrix b, save the result in c matrix */

void evaluate_matrix( int n, int m, POLY a[n][m], dcmplx b[n][m], dcmplx x);
/* evaluates the polynomial matrix a at the complex number x, save the 
   result in the b matrix */

void neg_polymatrix( int n, int m, POLY a[n][m] );
/*  returns a matrix with each elment has the negative coefficient of the original matrix */  
