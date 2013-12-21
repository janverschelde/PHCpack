/* file "dcmplx.h" defines complex numbers of doubles */
#include <stdio.h>

#ifndef DCMPLX_H
#define DCMPLX_H

typedef struct dcmplx dcmplx;

struct dcmplx
{
   double re;   /* real part of the complex number */
   double im;   /* imaginary part of the complex number */
};


static dcmplx zero = { 0.0, 0.0 };
static dcmplx one = { 1.0, 0.0 };

dcmplx create1 ( double r );
/* returns a complex number with real part r and zero imaginary part */

dcmplx create2 ( double r, double i );
/* returns a complex number with real part r and imaginary part i */

dcmplx random_dcmplx1 ( void );
/* returns random complex number of modulus one */

dcmplx polar ( double r, double a );
/* returns the complex number which has r*cos(a) and r*sin(a) 
   as real and imaginary parts */

void read_dcmplx ( dcmplx *z );
/* reads a complex number from standard input */

void read_dcmplx0 ( dcmplx *z, FILE *ifp );
/* reads a complex number with zero imaginary part from an input file with real format */ 

void read_dcmplx1 (dcmplx *z, FILE *ifp );
/* reads a complex number from an input file */

void write_dcmplx ( dcmplx z );
/* writes a complex number to the screen */

void write_dcmplx1 ( dcmplx z, FILE *ofp );
/* writes a complex number to a file */

void writeln_dcmplx ( dcmplx z );
/* writes a complex number to the screen, followed by newline symbol */

void writeln_dcmplx1 ( dcmplx z, FILE *ofp );
/* writes a complex number to a file, followed by newline symbol */

double dabs ( double x );
/* returns the absolute value of x */

double dcabs ( dcmplx z );
/* returns abs(z.re) + abs(z.im) */

int equal_dcmplx ( dcmplx z1, dcmplx z2, double tol );
/* returns 1 if | z1.re - z2.re | <= tol and | z1.im - z2.im | <= tol,
   returns 0 otherwise */

double modulus ( dcmplx z );
/* returns the modulus of the complex number z */

dcmplx conjugate ( dcmplx z );
/* returns the complex conjugate of z */

dcmplx min_dcmplx ( dcmplx z1 );              /* returns -z1 */
dcmplx add_dcmplx ( dcmplx z1, dcmplx z2 );   /* returns z1 + z2 */
dcmplx sub_dcmplx ( dcmplx z1, dcmplx z2 );   /* returns z1 - z2 */
dcmplx mul_dcmplx ( dcmplx z1, dcmplx z2 );   /* returns z1 * z2 */
dcmplx div_dcmplx ( dcmplx z1, dcmplx z2 );   /* returns z1 / z2 */

dcmplx add_double ( dcmplx z1, double z2 );   /* returns z1 + z2 */
dcmplx sub_double ( dcmplx z1, double z2 );   /* returns z1 - z2 */
dcmplx mul_double ( dcmplx z1, double z2 );   /* returns z1 * z2 */
dcmplx div_double ( dcmplx z1, double z2 );   /* returns z1 / z2 */

void swap (dcmplx* a, dcmplx* b);
/* swap the values of the two elements */

#endif
