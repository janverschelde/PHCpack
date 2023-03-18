/* file: double_double.h */

/* This file contains the header files for a standalone, self-contained
   collection of routines to provide double-double arithmetic,
   based on the QD-2.3.9 software library. */

#ifndef __double_double_h__
#define __double_double_h__

/* Part I: basic functions from inline.h */

double dd_quick_two_sum ( double a, double b, double *err );
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

double dd_two_sum ( double a, double b, double *err );
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

double dd_nint_d ( double d );
/*
 * DESCRIPTION :
 *   Returns nearest integer to d. */

double dd_aint ( double d );
/*
 * DESCRIPTION :
 *   Returns the truncated integer. */

/* Part II: arithmetical operations from dd_inline.h */

/* A double double x is represented as arrays of two doubles,
   the high word of x is in x[0] and the low word in x[1]. */

/*************************** additions ***************************/

/* double double = double double + double double */
void dd_add ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds two double doubles in a and b to make the double double c. */

/* double double = double + double double */
void dd_add_d_dd ( double a, const double *b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds the double a to the double double b to make the double double c. */

/* double_double = double double + double */
void dd_add_dd_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds the double double a to the double b to make the double double c. */

/* double double = double + double */
void dd_add_d_d ( double a, double b, double *c );
/*
 * DESCRIPTION : c = a + b, or in words:
 *   Adds the double a to the double b to make the double double c. */

/* double double += double double */
void dd_inc ( double *a, const double *b );
/*
 * DESCRIPTION : a += b, or in words:
 *   Inplace addition (increment) of the double double a
 *   with the double double b. */

/* double double += double */
void dd_inc_d ( double *a, double b );
/*
 * DESCRIPTION : a += b, or in words:
 *   Inplace addition (increment) of the double double a with the double b. */

/*************************** subtractions ***************************/

/* double double = double double - double double */
void dd_sub ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the double double in b from the double double in a
 *   and places the result in the double double c. */

/* double double = double - double double */
void dd_sub_d_dd ( double a, const double *b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the double double b from the double a 
 *   and places the result in the double double c. */

/* double double = double double - double */
void dd_sub_dd_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts from the double double a the double b 
 *   and places the result in the double double c. */

/* double double = double - double */
void dd_sub_d_d ( double a, double b, double *c );
/*
 * DESCRIPTION : c = a - b, or in words:
 *   Subtracts the double b from the double a to make the double double c. */

/* double double -= double double */
void dd_dec ( double *a, const double *b );
/*
 * DESCRIPTION : a -= b, or in words:
 *   Inplace difference (decrement) of the double double a
 *   with the double double b. */

/* double double -= double */
void dd_dec_d ( double *a, double b );
/*
 * DESCRIPTION : a -= b, or in words:
 *   Inplace subtraction (decrement) of the double b 
 *   from the double double a. */

/******************************* unary minus **************************/

/* double double = - double double */
void dd_minus ( double *a );
/*
 * DESCRIPTION : a = -a, or in words:
 *   Flips the sign of both low and high words of a. */

/*************************** multiplications ***************************/

/* double double = double double * double double */
void dd_mul ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a * b, or in words:
 *   Multiplies two double doubles in a and b to make the double double c. */

/* double double = double * double double */
void dd_mul_d_dd ( double a, const double *b, double *c );
/*
 * DESCRIPTION : c = a * b, or in words:
 *   Multiplies the double a to the double double b 
 *   to make the double double c. */

/* double double = double double * double */
void dd_mul_dd_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a * b, or in words:
 *   Multiplies the double double a to the double b 
 *   to make the double double c. */

/* double double = double * double */
void dd_mul_d_d ( double a, double b, double *c );
/*
 * DESCRIPTION : c = a * b, or in words:
 *   Multiplies the double a to the double b to make the double double c. */

/* double double *= double double */
void dd_mlt ( double *a, const double *b );
/*
 * DESCRIPTION : a *= b, or in words:
 *   Inplace multiplication of the double double a
 *   with the double double b. */

/* double double *= double */
void dd_mlt_d ( double *a, double b );
/*
 * DESCRIPTION : a *= b, or in words:
 *   Inplace multiplication of the double double a with the double b. */

/* double double = double double * double, where double is a power of 2. */
void dd_mul_pwr2 ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a*b as a plain multiplication
 *   of high and low word of a, because b is a power of 2. */

/*************************** divisions ***************************/

/* double double = double double / double double */
void dd_div ( const double *a, const double *b, double *c );
/*
 * DESCRIPTION : c = a / b, or in words:
 *   Divides the double doubles a by the double double b
 *   and places the result in the double double c. */

/* double double = double / double double */
void dd_div_d_dd ( double a, const double *b, double *c );
/*
 * DESCRIPTION : c = a / b, or in words:
 *   Divides the double a by the double double b 
 *   and places the result in the double double c. */

/* double double = double double / double */
void dd_div_dd_d ( const double *a, double b, double *c );
/*
 * DESCRIPTION : c = a / b, or in words:
 *   Divides the double double a by the double b 
 *   and places the result in the double double c. */

/* double double = double / double */
void dd_div_d_d ( double a, double b, double *c );
/*
 * DESCRIPTION : c = a * b, or in words:
 *   Divides the double a by the double b to make the double double c. */

/************************* comparisons *******************************/

int dd_is_zero ( const double *a );
/*
 * DESCRIPTION : a == 0, or in words:
 *   Returns 1 if a equals zero, returns 0 otherwise. */

int dd_is_one ( const double *a );
/*
 * DESCRIPTION : a == 1, or in words:
 *   Returns 1 if a equals one, returns 0 otherwise. */

int dd_is_positive ( const double *a );
/*
 * DESCRIPTION : a > 0, or in words:
 *   Returns 1 if a is positive, returns 0 otherwise. */

int dd_is_negative ( const double *a );
/*
 * DESCRIPTION : a < 0, or in words:
 *   Returns 1 if a is negative, returns 0 otherwise. */

/* double double == double */
int dd_eq_dd_d ( const double *a, double b );
/*
 * DESCRIPTION : a == b, or in words:
 *   Returns 1 if the double double a equals the double b,
 *   returns 0 otherwise. */

/* double double == double double */
int dd_eq ( const double *a, const double *b );
/*
 * DESCRIPTION : a == b, or in words:
 *   Returns 1 if the double double a equals the double double b,
 *   returns 0 otherwise. */

/* double == double double */
int dd_eq_d_dd ( double a, const double *b );
/*
 * DESCRIPTION : a == b, or in words:
 *   Returns 1 if the double a equals the double double b,
 *   returns 0 otherwise. */

/* double double != double */
int dd_neq_dd_d ( const double *a, double b );
/*
 * DESCRIPTION : a != b, or in words:
 *   Returns 1 if the double double a does not equal the double b,
 *   returns 0 otherwise. */

/* double double != double double */
int dd_neq ( const double *a, const double *b );
/*
 * DESCRIPTION : a != b, or in words:
 *   Returns 1 if the double double a does not equal the double double b,
 *   returns 0 otherwise. */

/* double != double double */
int dd_neq_d_dd ( double a, const double *b );
/*
 * DESCRIPTION : a != b, or in words:
 *   Returns 1 if the double a does not equal the double b,
 *   returns 0 otherwise. */

/* double double > double */
int dd_gt_dd_d ( const double *a, double b );
/*
 * DESCRIPTION : a > b, or in words:
 *   Returns 1 if the double double a is greater than the double b,
 *   returns 0 otherwise. */

/* double double > double double */
int dd_gt ( const double *a, double *b );
/*
 * DESCRIPTION : a > b, or in words:
 *   Returns 1 if the double double a is greater than the double double b,
 *   returns 0 otherwise. */

/* double > double double */
int dd_gt_d_dd ( double a, const double *b );
/*
 * DESCRIPTION : a > b, or in words:
 *   Returns 1 if the double a is greater than the double double b,
 *   returns 0 otherwise. */

/* double double >= double */
int dd_geq_dd_d ( const double *a, double b );
/*
 * DESCRIPTION : a >= b, or in words:
 *   Returns 1 if the double double a is greater than or equal to
 *   the double b, returns 0 otherwise. */

/* double double >= double double */
int dd_geq ( const double *a, const double *b ); 
/*
 * DESCRIPTION : a >= b, or in words:
 *   Returns 1 if the double double a is greater than or equal to
 *   the double double b, returns 0 otherwise. */

/* double >= double double */
int dd_geq_d_dd ( double a, const double *b );
/*
 * DESCRIPTION : a >= b, or in words:
 *   Returns 1 if the double a is greater than or equal to
 *   the double double b, returns 0 otherwise. */

/* double double < double */
int dd_lt_dd_d ( const double *a, double b );
/*
 * DESCRIPTION : a < b, or in words:
 *   Returns 1 if the double double a is less than the double b,
 *   returns 0 otherwise. */

/* double double < double double */
int dd_lt ( const double *a, const double *b );
/*
 * DESCRIPTION : a < b, or in words:
 *   Returns 1 if the double double a is less than the double double b,
 *   returns 0 otherwise. */

/* double < double double */
int dd_lt_d_dd ( double a, const double *b );
/*
 * DESCRIPTION : a < b, or in words:
 *   Returns 1 if the double a is less than the double double b,
 *   returns 0 otherwise. */

/* double double <= double */
int dd_leq_dd_d ( const double *a, double b );
/*
 * DESCRIPTION : a <= b, or in words:
 *   Returns 1 if the double double a is less than or equal to
 *   the double b, returns 0 otherwise. */

/* double double <= double double */
int dd_leq ( const double *a, const double *b );
/*
 * DESCRIPTION : a <= b, or in words:
 *   Returns 1 if the double double a is less than or equal to
 *   the double double b, returns 0 otherwise. */ 

/* double <= double double */
int dd_leq_d_dd ( double a, const double *b );
/*
 * DESCRIPTION : a <= b, or in words:
 *   Returns 1 if the double a is less than or equal to 
 *   the double double b, returns 0 otherwise. */

/******************** squaring and power ***************************/

void dd_sqr ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns in the double double b the square of the double double a. */

void dd_sqr_d ( double a, double *b );
/*
 * DESCRIPTION :
 *   Returns in the double double b the square of the double a. */

void dd_npwr ( const double *a, int n, double *b );
/*
 * DESCRIPTION :
 *   Returns in b the n-th power of the double double a. */

void dd_ldexp ( double *a, int n );
/*
 * DESCRIPTION :
 *   Multiplies the double double in a with 2^n. */

/*********************** exp and log ******************************/

void dd_exp ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns in b = exp(a).  Overflow returns -1.0. */

void dd_log ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns b = log(a).  If a <=0, b equals -1.0. */

void dd_log10 ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns b = log10(a).  If a <=0, b equals -1.0. */

/**************** copy, type casts, abs, floor, nint *****************/

void dd_copy ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Copies the content of the double double a to the double double b. */

int dd_to_int ( const double *a );
/*
 * DESCRIPTION :
 *   Casts the high word of a in a[0] to an integer. */

double dd_to_double ( const double *a );
/*
 * DESCRIPTION :
 *   Returns a[0], the high word of a. */

void dd_abs ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns in b the absolute value of the double double in a. */

void dd_floor ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns the largest integer to a as the double double b. */

void dd_nint ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns the nearest integer to a as the double double b. */

/************************** sqrt, sin and cos ***************************/  

void dd_sqrt ( const double* a, double* b );
/*
 * DESCRIPTION :
 *   Returns in the double double b 
 *   the square root of the double double a. */

void dd_sin_taylor ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns in b = sin(a) using Taylor series, assuming |a| <= pi/32. */ 

void dd_cos_taylor ( const double *a, double *b );
/*
 * DESCRIPTION :
 *   Returns in b = cos(a) using Taylor series. */ 

void dd_sincos_taylor ( const double *a, double *sin_a, double *cos_a );
/*
 * DESCRIPTION :
 *   Returns sin_a = sin(a), cos_a = cos(a) for |a| <= pi/32. */

void dd_reduce_modulo_2pi
 ( const double *x, double *t, int *j, int *k, int *abs_k, int* fail );
/*
 * DESCRIPTION :
 *   Reduces x modulo 2*pi, modulo pi/2, and then modulo pi/16.
 *   If this reduction does not work, then an error message is printed
 *   and fail is true on return.
 *
 * ON ENTRY :
 *   x        some double double number.
 *
 * ON RETURN :
 *   t        what remains of x;
 *   j        result after reduction modulo pi/2;
 *   k        result after reduction modulo pi/16;
 *   abs_k    absolute value of k;
 *   fail     true if reduction fails. */

void dd_sin ( const double *a, double *sin_a );
/*
 * DESCRIPTION :
 *   Returns sin_a = sin(a). */

void dd_cos ( const double *a, double *cos_a );
/*
 * DESCRIPTION :
 *   Returns sin_a = cos(a). */

/********************* input/output operations *************************/

int dd_read ( const char *s, double *a );
/*
 * DESCRIPTION :
 *   Reads a double double from the string s, the result is stored in a.
 *   If all goes well, 0 is returned, otherwise -1 is returned. */

void dd_to_digits ( const double *a, char *s, int *expn, int precision );
/*
 * DESCRIPTION :
 *   Converts the double double in a to decimal format,
 *   written to the string s, called by dd_to_string below.
 *   See the specification of to_string for more information. */

void dd_to_string ( const double *a, char *s, int precision, int width,
                    int fixed, int showpos, int uppercase, char fill,
                    int *endstring );
/*
 * DESCRIPTION :
 *   Writes the double double a to the string in s.
 *
 * ON ENTRY :
 *   a          a double double;
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
 *   s          string that represents the double double a;
 *   endstring  marks the end of the string: s[endstring] == '\0'. */

void dd_write ( const double *a, int precision );
/*
 * DESCRIPTION :
 *   Writes the double double to standard output, with given precision,
 *   and using default values for the dd_to_string. */

#endif /* _double_double_h */
