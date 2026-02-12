/* Collection of functions to split doubles. */

#ifndef __SPLITTING_DOUBLES_H__
#define __SPLITTING_DOUBLES_H__

#define uint64 unsigned long long int

void write_52bits ( int k, uint64 nbr );
/*
 * Writes the bits of the number nbr.
 * Initially, call with k = 52. */

void write_52double ( double nbr );
/*
 * Writes the 52 bits of the fraction of the number nbr,
 * followed by the exponent. */

uint64 last_bits ( int k, uint64 nbr );
/*
 * Returns the last k bits of the number nbr. */

void quarter_bits
 ( uint64 nbr, uint64 *b0, uint64 *b1, uint64 *b2, uint64 *b3, int vrblvl=0 );
/*
 * Splits the 52 bits in the number nbr in four parts,
 * each with 13 bits taken from the parts of nbr.
 * If the verbose level vrblvl > 0, then intermediate results are shown. */

void octo_split_bits
 ( uint64 nbr, uint64 *b0, uint64 *b1, uint64 *b2, uint64 *b3,
   uint64 *b4, uint64 *b5, uint64 *b6, uint64 *b7, int vrblvl=0 );
/*
 * Splits the 52 bits in the number nbr in eight parts,
 * the first four parts b0, b1, b2, b3 have each 7 bits taken out of nbr,
 * the next four parts b4, b5, b6, b7 have each 6 bits taken out of nbr.
 * If the verbose level vrblvl > 0, then intermediate results are shown. */

int leading_zeros ( uint64 nbr, int idxpwr, int vrblvl=0 );
/*
 * Returns the number of leading zeros in the number nbr,
 * relative to 2**idxpwr.
 * If vrblvl > 0, then the progression of the count is shown. */

double first_half ( double x, int vrblvl=0 );
/*
 * Returns the first 26 bits of x.
 * Assumes that x > 0. */

void half_split ( double x, double *x0, double *x1, int vrblvl=0 );
/*
 * Assuming x > 0, splits x into two doubles x0 and x1,
 * selecting the first 26 bits of the fraction of x to go into x0,
 * and the rest to go into x1.
 * If vrblvl > 0, then intermediate results are shown. */

void quarter_split
 ( double x, double *x0, double *x1, double *x2, double *x3, int vrblvl=0 );
/*
 * Assuming x > 0, splits x into four doubles x0, x1, x2, x3,
 * selecting the 13 bits of the fraction of x to go into x0,
 * then next 13 into x1, the next 13 into x2, and the rest in x3.
 * If vrblvl > 0, then intermediate results are shown. */

void octo_split
 ( double x, double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7, int vrblvl=0 );
/*
 * Assuming x > 0, splits x into eight doubles x0, x1, .., x7,
 * distributing the fraction of x into x0, x1, x2, x3, 7 bits each,
 * and then into x4, x5, x6, x7, 6 bits each.
 * If vrblvl > 0, then intermediate results are shown. */

bool is_balanced ( double x, double y, int threshold=13, int vrblvl=0 );
/*
 * Returns true if the difference in exponents between two consecutive
 * quarters x and y equals the threshold, or is less.
 * If vrblvl > 0, then the exponents are written. */

bool is_quarter_balanced ( double x, double y, int vrblvl=0 );
/*
 * Returns true if the difference in exponents between two consecutive
 * quarters x and y is 13, or less.
 * If vrblvl > 0, then the exponents are written. */

void balance ( double *x, double *y, int threshold=14, int vrblvl=0 );
/*
 * Makes x and y balanced by reduction of one bit of x and
 * addition of one bit to y, asumming both x and y are positive,
 * using the threshold exponent. */

void quarter_balance ( double *x, double *y, int vrblvl=0 );
/*
 * Makes x and y balanced by reduction of one bit of x and
 * addition of one bit to y, asumming both x and y are positive. */

void balance_quarters
 ( double *x0, double *x1, double *x2, double *x3, int vrblvl=0 );
/*
 * Balances the four quarters x0, x1, x2, x3 of a double.
 * If vrblvl > 0, output is written. */

void octo_balance
 ( double *x0, double *x1, double *x2, double *x3, 
   double *x4, double *x5, double *x6, double *x7, int vrblvl=0 );
/*
 * Balances the eight parts x0, x1, .., x7 of a double.
 * If vrblvl > 0, output is written. */

void make_exponent_zero ( double *x, double *pow2fac, int vrblvl=0 );
/*
 * Multiplies x with a factor pow2fac to make its exponent zero. */

#endif
