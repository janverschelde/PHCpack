/* This file double_double_functions.h defines the arithmetical operations 
   for double double numbers, defined by high and low doubles.

The algorithms are from the QD library of Hida, Li, and Bailey,
with the modification that a double double is not stored as an array
of two doubles, but plainly by two doubles: a high and a low double.
All functions have the prefix ddf_ so the dd_ functions can be applied
for comparison and for input and output. */

#ifndef __double_double_functions_h__
#define __double_double_functions_h__

/************************** additions ********************************/

double ddf_quick_two_sum ( double a, double b, double *err );
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

double ddf_two_sum ( double a, double b, double *err );
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

void ddf_add
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo );
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

double ddf_quick_two_diff ( double a, double b, double *err );
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

double ddf_two_diff ( double a, double b, double *err );
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

void ddf_minus ( double *a_hi, double *a_lo );
/*
 * DESCRIPTION : a = -a, unary minus,
 *   Flips the sign of both high and low parts of the double double a.
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_hi     low part of the double double a.
 *
 * ON RETURN :
 *   a_hi     high part of the double double -a;
 *   a_hi     low part of the double double -a. */

void ddf_sub
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo );
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

void ddf_sub_dd_d
 ( double a_hi, double a_lo, double b, double *c_hi, double *c_lo );
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

/********** incrementers, decrementers, and multipliers ****************/

void ddf_inc ( double *a_hi, double *a_lo, double b_hi, double b_lo );
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

void ddf_inc_d ( double *a_hi, double *a_lo, double b );
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

void ddf_dec ( double *a_hi, double *a_lo, double b_hi, double b_lo );
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

void ddf_dec_d ( double *a_hi, double *a_lo, double b );
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

void ddf_mlt ( double *a_hi, double *a_lo, double b_hi, double b_lo );
/*
 * DESCRIPTION : a = a * b.
 *   Inplace multiplication of the double double a (a_hi, a_lo)
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

void ddf_mlt_d ( double *a_hi, double *a_lo, double b );
/*
 * DESCRIPTION : a = a * b.
 *   Inplace multiplication of the double double a (a_hi, a_lo)
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

void ddf_split ( double a, double *hi, double *lo );
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

double ddf_two_prod ( double a, double b, double *err );
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

double ddf_two_sqr ( double a, double *err );
/*
 * DESCRIPTION :
 *   Computes fl(a*a) and err(a*a) faster than two_prod. */

void ddf_mul
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo );
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

void ddf_sqr ( double a_hi, double a_lo, double *b_hi, double *b_lo );
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

void ddf_mul_d_dd
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

void ddf_div
 ( double a_hi, double a_lo, double b_hi, double b_lo,
   double *c_hi, double *c_lo );
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

void ddf_sqrt ( double a_hi, double a_lo, double *b_hi, double *b_lo );
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

void ddf_abs ( double a_hi, double a_lo, double *b_hi, double *b_lo );
/*
 * DESCRIPTION :
 *   Returns in the double double b (b_hi, b_lo) the absolute value
 *   of the double double a (a_hi, a_lo).
 *
 * ON ENTRY :
 *   a_hi     high part of the double double a;
 *   a_lo     low part of the double double a.
 *
 * ON RETURN :
 *   b_hi     high part of the double double b;
 *   b_lo     low part of the double double b. */

#endif
