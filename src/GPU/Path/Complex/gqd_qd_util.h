// Conversions between real multiprecision CPU and GPU types
// and parsing command line arguments.

#ifndef __GQD_QD_UTIL__
#define __GQD_QD_UTIL__

#include "complex.h"
#include "complexH.h"
#include <iostream>

// copying the content of a to b

void qd2gqd ( double *a, double *b );
void gqd2qd ( double *a, double *b );

void qd2gqd ( qd_real *a, gqd_real *b );
void gqd2qd ( gqd_real *a, qd_real *b );

void qd2gqd ( dd_real *a, gdd_real *b ); 
void gqd2qd ( gdd_real *a, dd_real *b );

void comp1_gqd2qd(gd_complex* a, complexH<double>* b);
void comp1_gqd2qd(gdd_complex* a, complexH<dd_real>* b);
void comp1_gqd2qd(gqd_complex* a, complexH<qd_real>* b);

void comp1_qd2gqd(complexH<double>* a, gd_complex* b);
void comp1_qd2gqd(complexH<dd_real>* a, gdd_complex* b);
void comp1_qd2gqd(complexH<qd_real>* a, gqd_complex* b);


#endif /* __GQD_QD_UTIL__ */
