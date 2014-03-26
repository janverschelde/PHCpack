// Conversions between real multiprecision CPU and GPU types
// and parsing command line arguments.

#ifndef UTIL
#define UTIL

#include <qd/qd_real.h>
#include <cstdlib>
#include "complex.h"
#include "complexH.h"

// copying the content of a to b

void qd2gqd ( qd_real *a, gqd_real *b );
void gqd2qd ( gqd_real *a, qd_real *b );

void qd2gqd ( dd_real *a, gdd_real *b ); 
void gqd2qd ( gdd_real *a, dd_real *b ); 

void qd2gqd ( double *a, double *b );
void gqd2qd ( double *a, double *b );

int print_modes ( void );
// Prints the possible execution modes to screen.

int parse_arguments
 ( int argc, char *argv[], int *BS, int *dim, int *r, int *mode );
/*
   Parses the argc arguments on the command line in argv[].
   Returns 0 if okay, otherwise the function returns 1.
   If okay, then on return, BS contains the block size, dim the dimension,
   r the number of repeated runs, and mode the execution mode.  */

double random_double ( void ); // Returns a random double in [0,1].

#endif
