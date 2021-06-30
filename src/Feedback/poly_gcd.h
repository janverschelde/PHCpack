/* file poly_gcd.h contains the prototypes for the extended gcd method */

#include "dcmplx.h"
#include "poly_dcmplx.h"

POLY* ExtPolyGcd1 ( POLY a, POLY b );
/* returns the coefficients of the extended gcd with an array of polynomial
   k*a+l*b=gcd; the items in the returned array are k, l, gcd, b/gcd, a/gcd
   respectively   */

dcmplx **ExtPolyGcd 
 ( int n, dcmplx a[n], int m, dcmplx b[m], int *dgcd, int *dk, int *dl,
   int *dbd, int *dad );
/* returns the coefficients of the extended gcd method with an array 
 * of complex pointer */

void rootsGCD
 ( int n, dcmplx a[n], int m, dcmplx b[m], int *l,
   dcmplx ra[n-1], dcmplx rb[m-1] );
/* get the common roots of polynomial a and b, then save them as 
   the first l elements in the ra and rb */

POLY get_gcd ( POLY a, POLY b );
/* returns the gcd of the polynomials a and b */

void extended_gcd
 ( int n, dcmplx ra[n], int m, dcmplx rb[m], int c, 
   dcmplx k[m-c], dcmplx l[n-c] );
/* gets the coefficient of the extended gcd method and
 * saves them in the k and l */

dcmplx *get_poly ( int n, dcmplx * root );
/* given an array of roots of a polynomial, returns a polynomial 
 * with the leading coefficient is one */

void free_gcd_coeff ( dcmplx ** gcd_coeff );
/* free the memory allocated for ExtPolyGcd function */

void free_gcd_coeff1(POLY * result);
/* free the memory allocated for ExtPolyGcd1 function */

int group_points ( int n, double tol, dcmplx x[n], dcmplx f[n], int m[n] );
/* makes the multiple points adjacent each other;
 * if there exists multiple points, returns 1; otherwise, returns 0 */

dcmplx **rational_derivative
 ( int n, dcmplx *num, int m, dcmplx *denom, int *d_num, int *d_denom );
/* find the derivative of the rational polynomial */

void hermite_derivative
 ( int root_num, dcmplx root[root_num], int p_num,
   dcmplx x[p_num], dcmplx f[p_num], int m[p_num] );
/* this is a subroutine for Hermite interpolation algorithm, 
 * if some multiple points are found, say x[i] and x[i+1] have the same value,
 * f[i] will save the function value at the point x[i], f[i+1] will be used
 * to save the first order derivative of the function, and so on. */
