/* file: quad_double.h */

/* This file contains the header files for a standalone, self-contained
   collection of routines to provide quad-double arithmetic, based on
   the QD-2.3.9 software library. */

/* A quad double is represented by an array of four doubles. */

#ifndef _quad_double_h
#define _quad_double_h

/* Part I: basic functions from qd_inline.h */

/********************* renormalizations ***************************/

void qd_quick_renorm ( double *c0, double *c1, double *c2,
                       double *c3, double *c4 );
/*
 * DESCRIPTION :
 *   Does a quick renormalization of the number
 *   represented by the five given doubles. */

void qd_renorm4 ( double *c0, double *c1, double *c2, double *c3 );
/*
 * DESCRIPTION :
 *   Does a renormalization of the number
 *   represented by the four given doubles. */

void qd_renorm5 ( double *c0, double *c1, double *c2,
                  double *c3, double *c4 );
/*
 * DESCRIPTION :
 *   Does a renormalization of the number
 *   represented by the five given doubles. */

/************************* additions *****************************/

void qd_three_sum ( double *a, double *b, double *c );
/*
 * DESCRIPTION :
 *   Generalization of the dd_two_sum. */

void qd_three_sum2 ( double *a, double *b, double *c );
/*
 * DESCRIPTION :
 *   Less accurate version of qd_three_sum. */

/* quad double = quad double + double */
void qd_add_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds the quad double a to the double in b 
 *   to make the quad double c. */

/* quad double = quad double + double double */
void qd_add_dd ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds the quad double a to the double double in b 
 *   to make the quad double c. */

double qd_quick_three_accum ( double *a, double *b, double c );
/*
 * DESCRIPTION :
 *   Adds c to the double double pair (a,b).
 *   If the result does not fit in two doubles, 
 *   then the sum is returned and (a,b) contains the remainder.
 *   Otherwise, the returned value is zero
 *   and on return (a,b) stores the sum. */

/* quad double = quad double + quad double */
void qd_add ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds the quad double a to the quad double b
 *   and stores the result in the quad double c. */

/******* constructor, copy, abs, type casts, and unary minus ********/

void qd_real ( double a, double *b );
/*
 * DESCRIPTION :
 *   Copies the double a to the first word of the quad double b,
 *   setting all other words in b to zero. */

void qd_copy ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Copies the four doubles from a to b. */

void qd_abs ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns in b the absolute value of the quad double in a. */

int qd_to_int ( const double *a );
/*
 * DESCRIPTION :
 *   Casts the first word a[0] to an integer. */

double qd_to_double ( const double *a );
/*
 * DESCRIPTION :
 *   Returns a[0], the first word of the quad double a. */

void qd_floor ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns the largest integer to the quad double a
 *   in the quad double b. */

/* quad double = - quad double */
void qd_minus ( double *a );
/*
 * DESCRIPTION : a = -a, or in words:
 *   Flips the sign of all components of a. */

/*********************** subtractions **************************/

/* quad double = quad double - quad double */
void qd_sub ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Substracts the quad double b from the quad double a
 *   and stores the result in the quad double c. */

/* quad double = quad double - double */
void qd_sub_qd_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the double b from the quad double a
 *   and stores the result in the quad double c. */

/* quad double = double - quad double */
void qd_sub_d_qd ( double a, const double *b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the double double b from the double a
 *   and stores the result in the quad double c. */

/* quad double = quad double - double double */
void qd_sub_qd_dd ( const double *a, double *b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the double double b from the quad double a
 *   and stores the result in the quad double c. */

/* quad double = double double - quad double */
void qd_sub_dd_qd ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the quad double b from the double double a
 *   and stores the result in the quad double c. */

/************************* multiplications *************************/

/* quad double = quad double * double, where double is a power of 2 */
void qd_mul_pwr2 ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a*b, as a plain multiplication of all
 *   components of a with b, because b is a power of 2. */

/* quad double = quad double * double */
void qd_mul_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a*b, or in words:
 *   Multiplies the quad double a with the double b
 *   and stores the result in the quad double c. */

/* quad double = quad double * double double */
void qd_mul_dd ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a*b, or in words:
 *   Multiplies the quad double a with the double double b
 *   and stores the result in the quad double c. */

/* quad double = quad double * quad double */
void qd_mul ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a*b, or in words:
 *   Multiplies the quad double a with the quad double b
 *   and stores the result in the quad double c. */

/********************** comparisons *****************************/

int qd_is_zero ( const double *a );
/*
 * DESCRIPTION : a == 0, or in words:
 *   Returns 1 if a equals zero, returns 0 otherwise. */ 

int qd_is_one ( const double *a );
/*
 * DESCRIPTION : a == 1, or in words:
 *   Returns 1 if a equals one, returns 0 otherwise. */ 

int qd_is_positive ( const double *a );
/*
 * DESCRIPTION : a > 0, or in words:
 *   Returns 1 if a is positive, returns 0 otherwise. */

int qd_is_negative ( const double *a );
/*
 * DESCRIPTION : a < 0, or in words:
 *   Returns 1 if a is negative, returns 0 otherwise. */

/* quad double == double */
int qd_eq_qd_d ( const double *a, double b );
/*
 * DESCRIPTION : a == b, or in words:
 *   Returns 1 if the quad double a equals the double b,
 *   returns 0 otherwise. */

/* quad double == double double */
int qd_eq_qd_dd ( const double *a, const double *b );
/*
 * DESCRIPTION : a == b, or in words:
 *   Returns 1 if the quad double a equals the double double b,
 *   returns 0 otherwise. */

/* quad double == quad double */
int qd_eq ( const double *a, const double *b );
/*
 * DESCRIPTION : a == b, or in words:
 *   Returns 1 if the quad double a equals the quad double b,
 *   returns 0 otherwise. */

/* quad double < double */
int qd_lt_qd_d ( const double *a, double b );
/*
 * DESCRIPTION : a < b, or in words:
 *   Returns 1 if the quad double a is less than the double b,
 *   returns 0 otherwise. */

/* quad double < double double */
int qd_lt_qd_dd ( const double *a, const double *b );
/*
 * DESCRIPTION : a < b, or in words:
 *   Returns 1 if the quad double a is less than the double double b,
 *   returns 0 otherwise. */

/* quad double < quad double */
int qd_lt ( const double *a, const double *b );
/*
 * DESCRIPTION : a < b, or in words:
 *   Returns 1 if the quad double a is less than the quad double b,
 *   returns 0 otherwise. */

/* quad double <= double */
int qd_leq_qd_d ( const double *a, double b );
/*
 * DESCRIPTION : a <= b, or in words:
 *   Returns 1 if the quad double a is less than or equal to the double b,
 *   returns 0 otherwise. */

/* quad double <= double double */
int qd_leq_qd_dd ( const double *a, const double *b );
/*
 * DESCRIPTION : a <= b, or in words:
 *   Returns 1 if the quad double a is less than or equal to
 *   the double double b, returns 0 otherwise. */

/* quad double <= quad double */
int qd_leq ( const double *a, const double *b );
/*
 * DESCRIPTION : a <= b, or in words:
 *   Returns 1 if the quad double a is less than or equal to 
 *   the quad double b, returns 0 otherwise. */

/* quad double > double */
int qd_gt_qd_d ( const double *a, double b );
/*
 * DESCRIPTION : a > b, or in words:
 *   Returns 1 if the quad double a is greater than the double b,
 *   returns 0 otherwise. */

/* quad double > double double */
int qd_gt_qd_dd ( const double *a, const double *b );
/*
 * DESCRIPTION : a > b, or in words:
 *   Returns 1 if the quad double a is greater than the double double b,
 *   returns 0 otherwise. */

/* quad double > quad double */
int qd_gt ( const double *a, const double *b );
/*
 * DESCRIPTION : a > b, or in words:
 *   Returns 1 if the quad double a is greater than the quad double b,
 *   returns 0 otherwise. */

/* quad double >= double */
int qd_geq_qd_d ( const double *a, double b );
/*
 * DESCRIPTION : a >= b, or in words:
 *   Returns 1 if the quad double a is greater than or equal to the double b,
 *   returns 0 otherwise. */

/* quad double >= double double */
int qd_geq_qd_dd ( const double *a, const double *b );
/*
 * DESCRIPTION : a >= b, or in words:
 *   Returns 1 if the quad double a is greater than or equal to
 *   the double double b, returns 0 otherwise. */

/* quad double >= quad double */
int qd_geq ( const double *a, const double *b );
/*
 * DESCRIPTION : a >= b, or in words:
 *   Returns 1 if the quad double a is greater than or equal to
 *   the quad double b, returns 0 otherwise. */

/********** Part II : operations from qd_real.cpp ****************/

/************************ divisions ******************************/

/* quad double = quad double / double */
void qd_div_qd_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a/b, or in words:
 *   Divides the quad double a by the double b
 *   and stores the result in the quad double c. */

/* quad double = quad double / double double */
void qd_div_qd_dd ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a/b, or in words:
 *   Divides the quad double a by the double double b
 *   and stores the result in the quad double c. */

/* quad double = quad double / quad double */
void qd_div ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a/b, or in words:
 *   Divides the quad double a by the quad double b
 *   and stores the result in the quad double c. */

/* quad double = double / quad double */
void qd_div_d_qd ( double a, const double *b, double *c );
/*
 * DESCRIPTION : c = a/b, or in words:
 *   Divides the double a by the quad double b
 *   and stores the result in the quad double c. */

/* quad double = double double / quad double */
void qd_div_dd_qd ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a/b, or in words:
 *   Divides the double double a by the quad double b
 *   and stores the result in the quad double c. */

/******************* squaring and power *************************/

void qd_sqr ( const double *a, double *b );
/*
 * DESCRIPTION : b = a^2, or in words:
 *   Takes the square of the double double a
 *   and stores the result in b. */

void qd_npwr ( const double *a, int n, double *b );
/*
 * DESCRIPTION : b = a^n, or in words:
 *   Returns in b the n-th power of the quad double in a. */

void qd_ldexp ( double *a, int n );
/*
 * DESCRIPTION : a = a*2^n, or in words:
 *   Multiplies the quad double in a with 2^n. */

/*********************** exp and log ******************************/

void qd_exp ( const double *a, double *b );
/*
 * DESCRIPTION : b = exp(a), or in words:
 *   Computes the exponential of the quad double in a
 *   and stores the result in b.  Overflow returns -1. */

void qd_log ( const double *a, double *b );
/*
 * DESCRIPTION : b = log(a), or in words:
 *   Computes the natural logarithm of the quad double in a
 *   and stores the result in b.  If a <= 0, then b is -1. */

void qd_log10 ( const double *a, double *b );
/*
 * DESCRIPTION : b = log10(a), or in words:
 *   Computes the decimal logarithm of the quad double in a
 *   and stores the result in b.  If a <= 0, then b is -1. */

/********** PART III : input/output operations *****************/

int qd_read ( const char *s, double *a );
/*
 * DESCRIPTION :
 *   Reads a quad double from the string s, the result is stored in a.
 *   If all goes well, 0 is returned, otherwise -1 is returned. */

void qd_to_digits ( const double *a, char *s, int *expn, int precision );
/*
 * DESCRIPTION :
 *   Converts the quad double in a to decimal format,
 *   written to the string s, called by qd_to_string below.
 *   See the specification of to_string for more information. */

void qd_to_string ( const double *a, char *s, int precision, int width,
                    int fixed, int showpos, int uppercase, char fill,
                    int *endstring );
/*
 * DESCRIPTION :
 *   Writes the quad double a to the string in s.
 *
 * ON ENTRY :
 *   a          a quad double;
 *   s          a string large enough to hold at least d+8 characters,
 *              where d is the number of significant digits to write;
 *   precision  equals the number of significant digits to write;
 *   width      if large enough, the string s will be padded with
 *              characters equals to fill;
 *   fixed      if 1, then fixed format is used,
 *              otherwise, a is displayed in scientific format;
 *   showpos    if 1, then the + sign will be shown for positive numbers;
 *   uppercase  if 1, then E, otherwise e is used for the exponent;
 *   fill       character used for padding the string in case width
 *              is large enough.
 *
 * ON RETURN :
 *   s          string that represents the quad double a;
 *   endstring  marks the end of the string: s[endstring] == '\0'. */

void qd_write ( const double *a, int precision );
/*
 * DESCRIPTION :
 *   Writes the quad double to standard output, with given precision,
 *   and using default values for the qd_to_string. */

#endif /* _quad_double_h */
