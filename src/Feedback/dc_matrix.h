/* file "dc_matrix.h" defines the matrices operation 
 * on complex number of doubles*/
#include "dcmplx.h"

void print_dcmatrix ( int n, int m, dcmplx a[n][m] );
/* writes an n-by-m dcmplx matrix to the screen */

void print_dcmatrix1 ( int n, int m, dcmplx a[n][m], FILE *ofp);
/* writes an n-by-m dcmplx matrix to a file*/

void read_dcmatrix ( int n, int m, dcmplx a[n][n] );
/* reads an n-by-m matrix a of double floating-point complex numbers */

void read_dcmatrix0 ( int n, int m, dcmplx a[n][m], FILE *ifp );
/* reads an n-by-m matrix a of double floating-point complex numbers from an input file with real format*/

void read_dcmatrix1 ( int n, int m, dcmplx a[n][n] );
/* reads an n-by-m matrix a of double floating-point complex numbers */

void read_dcmatrix2 ( int n, int m, dcmplx a[n][m], FILE *ifp );
/* reads an n-by-m matrix a of double floating-point complex numbers from an input file */

void random_dcmatrix ( int n, int m, dcmplx a[n][m] );
/* creates a random n-by-m matrix */

void random_dcmatrix0 ( int n, int m, dcmplx a[n][m] );
/* creates a random n-by-m matrix with the imaginary part of each entry is 0 */

void copy_dcmatrix ( int n, int m, dcmplx a[n][m], dcmplx b[n][m] );
/* copies the n-by-m matrix a to b */

void multiply_dcmatrix ( int n, int m, int l, dcmplx a[n][m], dcmplx b[m][l], dcmplx c[n][l] );
/* multiply matrices a and b, put the result in c matrix */

void add_dcmatrix ( int n, int m, dcmplx a[n][m], dcmplx b[n][m], dcmplx c[n][m] );
/* add matrices a and b, save the result in c matrix */

void sub_dcmatrix ( int n, int m, dcmplx a[n][m], dcmplx b[n][m], dcmplx c[n][m]);
/* subtracts matrices b from a, save the result in c matrix */

void free_dcmatrix ( int n, dcmplx ** a);
/* free the memory assigned to a */

int lufac ( int n, dcmplx a[n][n], int ipvt[n] );
/* computes an LU-factorization of an n-by-n matrix a,
   the matrix a on return contains the factors l and u,
   while the pivoting information is stored in ipvt;
   the normal return value of lufac is -1, otherwise the
   value on return is the column containing a zero pivot */

void lusolve ( int n, dcmplx a[n][n], int ipvt[n], dcmplx b[n] );
/* solves the system a*x = b using the factors computed by lufac,
   the matrix a and array ipvt are the output of lufac, the array
   b contains the solution on return */

dcmplx determinant( int n, dcmplx (*a)[n] );
/* returns the determinant of matrix a */

void I_dcmatrix(int n, dcmplx a[n][n]);
/* returns an n dimension identity matrix */
