/* file "pieri_sols.h" contains the prototype for transfering a polynomial matrix from
   ada to c format */
#include "poly_dcmplx.h"

void skip(FILE *ifp);
/* skips the dimensions information */

void ada2c_poly_matrix( int n, int m, POLY c[n][m], int l, int start_a, int a[l], int q, int start_b, double b[q] );
/* read the polynomial of each matrix and save all the information in matrix c */

int pieri_sols( int m, int p, int q, int nbsols, int da, int a[da], int db, double b[db],
                 int npt, double pt[npt], int npl, double pl[npl],
                 char *filename );
/* read an integer array and a double array from ada subroutine, convert to the input matrices (M1 and M2) of
   realization algorithm */ 
