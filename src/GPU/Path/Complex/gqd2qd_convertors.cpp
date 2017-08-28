// Defines the functions with prototypes in gqd2qd_convertors.h. 

#include "gqd2qd_convertors.h"

void real_d2gd ( double *a, double *b )
{
   *b = *a;
}

void real_gd2d ( double *a, double *b )
{
   *b = *a;
}

void real_dd2gdd ( dd_real *a, gdd_real *b )
{
   b->x = a->x[0];
   b->y = a->x[1];
}

void real_gdd2dd ( gdd_real *a, dd_real *b )
{
   b->x[0] = a->x;
   b->x[1] = a->y;
}


void real_qd2gqd ( qd_real *a, gqd_real *b ) 
{
   b->x = a->x[0];
   b->y = a->x[1];
   b->z = a->x[2];
   b->w = a->x[3];
}

void real_gqd2qd ( gqd_real *a, qd_real *b )
{
   b->x[0] = a->x;
   b->x[1] = a->y;
   b->x[2] = a->z;
   b->x[3] = a->w;
}

void complex_d2gd ( complexH<double> *a, complexD<double> *b )
{
   b->real = a->real;
   b->imag = a->imag;
}

void complex_gd2d ( complexD<double> *a, complexH<double> *b )
{
   b->real = a->real;
   b->imag = a->imag;
}

void complex_dd2gdd ( complexH<dd_real> *a, complexD<gdd_real> *b )
{
   real_dd2gdd(&(a->real),&(b->real));
   real_dd2gdd(&(a->imag),&(b->imag));
}

void complex_gdd2dd ( complexD<gdd_real> *a, complexH<dd_real> *b )
{
   real_gdd2dd(&(a->real),&(b->real));
   real_gdd2dd(&(a->imag),&(b->imag));
}

void complex_qd2gqd ( complexH<qd_real> *a, complexD<gqd_real> *b )
{
   real_qd2gqd(&(a->real),&(b->real));
   real_qd2gqd(&(a->imag),&(b->imag));
}

void complex_gqd2qd ( complexD<gqd_real> *a, complexH<qd_real> *b )
{
   real_gqd2qd(&(a->real),&(b->real));
   real_gqd2qd(&(a->imag),&(b->imag));
}
