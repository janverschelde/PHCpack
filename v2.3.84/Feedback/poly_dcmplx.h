#ifndef POLY_DCMPLX_H
#define POLY_DCMPLX_H

#include "dcmplx.h"

typedef struct poly POLY;

struct poly
{
   int d;                /* degree of the polynomial */
   dcmplx *p;            /* coefficient of the polynomial */
};


void Read_Poly1 ( int n, dcmplx a[] );
/* read one polynomial from input */

void Read_Poly ( int n, dcmplx a[], int m, dcmplx b[] );
/* read two polynomials from input */

POLY random_poly (int max_deg);
/* get a random polynomial with the maximum degree is max_deg */

POLY random_poly1 (int deg);
/* get a random polynomial with the degree specified */

void Print_Poly ( int k, dcmplx * c );
/* print a polynomial without "\n" at the end */

POLY neg_poly ( POLY a);
/* return a polynomial with the negative coefficient of a  */

dcmplx *add_poly ( int n, dcmplx a[], int m, dcmplx b[], int *size_sum );
/* returns the coefficient array of the sum of the polynomials a and b,
   the size of this array is in size_sum */

POLY add_poly1 ( POLY a, POLY b);
/* return polynomial a+b */

dcmplx *min_poly ( int n, dcmplx a[], int m, dcmplx b[], int *size_sum );
/* returns the coefficient array of the polynomials a - b,
   the size of this array is in size_sum */

POLY min_poly1( POLY a, POLY b);
/* return polynomial a-b */

dcmplx* mul_poly ( int n, dcmplx a[], int m, dcmplx b[], int *size_sum );
/* returns the coefficient array of the polynomials a*b,
   the size of this array is in size_sum */

POLY mul_poly1 ( POLY a, POLY b);
/* return polynomial a*b */

POLY mul_dcmplx_poly ( dcmplx a, POLY b );
/* return the product of the complex number a and polynomial b */

dcmplx** div_poly ( int n, dcmplx a[], int m, dcmplx b[],
			    int* dq, int* dr);
/* */

dcmplx* div_poly2 ( int n, dcmplx a[], int m, dcmplx b[], int *dq );
/* Can only be used when the remainder is zero */

int degree ( dcmplx *a, int d );

int iszero ( int n,  dcmplx a[] );

void Test_Div ( int n, int m );

/*dcmplx dabs (dcmplx x); */ 

dcmplx * assign(int n, dcmplx a[]);
/* returns a polynomial which is same as a */

POLY assign_poly(POLY a);
/* returns a polynomial which is same as a*/

void negative(int n, dcmplx a[]);

int equal_poly(int n, dcmplx a[], int m, dcmplx b[]);
/* returns 1 if a and b are equal, otherwise return 0 */

void divide_by_number( int n, dcmplx a[], dcmplx num);
/*  divide all the coefficient of a  by a complex number */

void mult_by_number( int n, dcmplx a[], dcmplx num);
/*  multiply all the coefficient of a by a complex number */

#endif











