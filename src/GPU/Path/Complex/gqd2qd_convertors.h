// conversions between real multiprecision CPU and GPU types
// and the complex data types

#ifndef __GQD2QD_CONVERTORS__
#define __GQD2QD_CONVERTORS__

#include <qd/qd_real.h>
#include "complexH.h"
#include "complexD.h"

// copying the content of a to b

void real_d2gd ( double *a, double *b ); 
// double a on host to double b on device
void real_gd2d ( double *a, double *b );
// double a on device to double b on host

void real_dd2gdd ( dd_real *a, gdd_real *b );
// double double a on host to double double b on device
void real_gdd2dd ( gdd_real *a, dd_real *b );
// double double a on device to double double b on host

void real_qd2gqd ( qd_real *a, gqd_real *b ); 
// quad double a on host to quad double b on device
void real_gqd2qd ( gqd_real *a, qd_real *b );
// quad double a on device to quad double b on host

void complex_d2gd ( complexH<double> *a, complexD<double> *b ); 
// complex double a on host to complex double b on device
void complex_gd2d ( complexD<double> *a, complexH<double> *b );
// complex double a on device to complex double b on host

void complex_dd2gdd ( complexH<dd_real> *a, complexD<gdd_real> *b );
// complex double double a on host to complex double double b on device
void complex_gdd2dd ( complexD<gdd_real> *a, complexH<dd_real> *b );
// complex double double a on device to complex double double b on host

void complex_qd2gqd ( complexH<qd_real> *a, complexD<gqd_real> *b ); 
// complex quad double a on host to complex quad double b on device
void complex_gqd2qd ( complexD<gqd_real> *a, complexH<qd_real> *b );
// complex quad double a on device to complex quad double b on host

#endif /* __GQD2QD_CONVERTORS__ */
