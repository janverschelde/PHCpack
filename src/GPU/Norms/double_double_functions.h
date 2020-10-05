// This file double_double_functions.h defines the arithmetical operations 
// for double double numbers, defined by high and low doubles.

#ifndef __double_double_functions_h__
#define __double_double_functions_h__

/************************** additions ********************************/

double dd_quick_two_sum ( double a, double b, double* err );
/*
 * DESCRIPTION :
 *   Assuming |a| >= |b|, returns a+b and in err the error.
 *
 * ON ENTRY :
 *   a,b      two doubles: |a| >= |b|.
 *
 * ON RETURN :
 *   s        returned sum of a and b.
 *   err      error value, b - (s - a). */

double dd_two_sum ( double a, double b, double* err );
/*
 * DESCRIPTION :
 *   Computes fl(a+b) and err(a+b).
 *
 * ON ENTRY :
 *   a,b      two doubles.
 *
 * ON RETURN :
 *   s        approximation for the sum of a and b is returned;
 *   err      error of a + b. */

void dd_add
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double* c_hi, double* c_lo );
/*
 * DESCRIPTION : c = a + b.
 *   Adds two double doubles in a (a_hi, a_lo) and b (b_hi, b_lo)
 *   to make the double double c (c_hi, c_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b.
 *
 * ON RETURN :
 *   c_hi     high part of the double double c;
 *   c_lo     low part of the double double c. */

double dd_quick_two_diff ( double a, double b, double *err );
/*
 * DESCRIPTION :
 *   Assuming |a| >= |b|, returns a-b and in err the error.
 *
 * ON ENTRY :
 *   a,b      two doubles: |a| >= |b|.
 *
 * ON RETURN :
 *   s        returned a minus b.
 *   err      error value, (a - s) - b. */

double dd_two_diff ( double a, double b, double *err );
/*
 * DESCRIPTION :
 *   Computes fl(a-b) and err(a-b).
 *
 * ON ENTRY :
 *   a,b      two doubles.
 *
 * ON RETURN :
 *   s        approximation for the difference of a with b is returned;
 *   err      error of a - b. */

void dd_sub
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double* c_hi, double* c_lo );
/*
 * DESCRIPTION : c = a - b.
 *   Subtracts the double double in b (b_hi, b_lo) 
 *   from the double double in a (a_hi, a_lo)
 *   and places the result in the double double c (c_hi, c_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b.
 *
 * ON RETURN :
 *   c_hi     high part of the double double c;
 *   c_lo     low part of the double double c. */

void dd_sub_dd_d
 ( double a_hi, double a_lo, double b, double* c_hi, double* c_lo );
/*
 * DESCRIPTION : c = a - b.
 *   Subtracts the double b from the double double in a (a_hi, a_lo)
 *   and places the result in the double double c (c_hi, c_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b        some double.
 *
 * ON RETURN :
 *   c_hi     high part of the double double c;
 *   c_lo     low part of the double double c. */

void dd_inc ( double *a_hi, double *a_lo, double b_hi, double b_lo );
/*
 * DESCRIPTION : a = a + b.
 *   Inplace increment of the double double a (a_hi, a_lo)
 *   with the double double in b (b_hi, b_lo) 
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b.
 *
 * ON RETURN :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a. */

void dd_inc_d ( double *a_hi, double *a_lo, double b );
/*
 * DESCRIPTION : a = a + b.
 *   Inplace increment of the double double a (a_hi, a_lo)
 *   with the double b.
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b        some double.
 *
 * ON RETURN :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a. */

void dd_dec ( double *a_hi, double *a_lo, double b_hi, double b_lo );
/*
 * DESCRIPTION : a = a - b.
 *   Inplace decrement of the double double a (a_hi, a_lo)
 *   with the double double in b (b_hi, b_lo) 
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b.
 *
 * ON RETURN :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a. */

void dd_dec_d ( double *a_hi, double *a_lo, double b );
/*
 * DESCRIPTION : a = a - b.
 *   Inplace decrement of the double double a (a_hi, a_lo)
 *   with the double b.
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b        some double.
 *
 * ON RETURN :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a. */

/************************ multiplications ********************************/

void dd_split ( double a, double *hi, double *lo );
/*
 * DESCRIPTION :
 *   Computes high and low word of a.
 *
 * ON ENTRY :
 *   a        some double float.
 *
 * ON RETURN :
 *   hi       high word of a;
 *   lo       low word of a. */ 

double dd_two_prod ( double a, double b, double *err );
/*
 * DESCRIPTION :
 *   Computes fl(a*b) and err(a*b).
 *
 * ON ENTRY :
 *   a,b      two doubles.
 *
 * ON RETURN :
 *   p        returned approximation for a*b;
 *   err      error on the approximated product. */

double dd_two_sqr ( double a, double *err );
/*
 * DESCRIPTION :
 *   Computes fl(a*a) and err(a*a) faster than two_prod. */

void dd_mul
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double* c_hi, double* c_lo );
/*
 * DESCRIPTION : c = a * b.
 *   Multiplies two double doubles in a (a_hi, a_lo) and b (b_hi, b_lo)
 *   to make the double double c (c_hi, c_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b.
 *
 * ON RETURN :
 *   c_hi     high part of the double double c;
 *   c_lo     low part of the double double c. */

void dd_sqr ( double a_hi, double a_lo, double *b_hi, double *b_lo );
/*
 * DESCRIPTION :
 *   Returns in the double double b (b_hi, b_lo) 
 *   the square of the double double a (a_hi, a_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a.
 *
 * ON RETURN :
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b. */

void dd_mul_d_dd
 ( double a, double b_hi, double b_lo, double *c_hi, double *c_lo );
/*
 * DESCRIPTION : c = a * b.
 *   Multiplies the double a with the double double b (b_hi, b_lo)
 *   to make the double double c (c_hi, c_lo).
 *
 * ON ENTRY :
 *   a        some double;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b;
 *
 * ON RETURN :
 *   c_hi     high part of the double double c;
 *   c_lo     low part of the double double c. */

/*************************** divisions ***************************/

void dd_div
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double* c_hi, double* c_lo );
/*
 * DESCRIPTION : c = a / b.
 *   Divides the double doubles in a (a_hi, a_lo) by b (b_hi, b_lo)
 *   to make the double double c (c_hi, c_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a;
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b.
 *
 * ON RETURN :
 *   c_hi     high part of the double double c;
 *   c_lo     low part of the double double c. */

/*************************** sqrt ***************************/

void dd_sqrt ( double a_hi, double a_lo, double *b_hi, double *b_lo );
/*
 * DESCRIPTION :
 *   Returns in the double double b (b_hi, b_lo) 
 *   the square root of the double double a (a_hi, a_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a.
 *
 * ON RETURN :
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b. */

#endif
