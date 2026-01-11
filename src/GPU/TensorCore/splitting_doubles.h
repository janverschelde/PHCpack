/* Collection of functions to split doubles. */

#ifndef __SPLITTING_DOUBLES_H__
#define __SPLITTING_DOUBLES_H__

#define uint64 unsigned long long int

void write_52bits ( int k, uint64 nbr );
/*
 * Writes the bits of the number nbr.
 * Initially, call with k = 52. */

uint64 last_bits ( int k, uint64 nbr );
/*
 * Returns the last k bits of the number nbr. */

double first_half ( double x, int vrblvl );
/*
 * Returns the first 26 bits of x. */

void half_split ( double x, double *x0, double *x1, int vrblvl );
/*
 * Splits x into two doubles x0 and x1,
 * selecting the first 26 bits of x0. */

int test ( void );
/*
 * Generates a random number, splits, and then checks if
 * adding the parts gives the original number. */

#endif
