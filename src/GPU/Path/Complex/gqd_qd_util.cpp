// Conversions between real multiprecision CPU and GPU types
// and parsing command line arguments as defined in gqd_qd_util.h.

#include "gqd_qd_util.h"

using namespace std;

void qd2gqd ( qd_real *a, gqd_real *b )
{
   b->x = a->x[0];
   b->y = a->x[1];
   b->z = a->x[2];
   b->w = a->x[3];
}

void gqd2qd ( gqd_real *a, qd_real *b )
{
   b->x[0] = a->x;
   b->x[1] = a->y;
   b->x[2] = a->z;
   b->x[3] = a->w;
}

void qd2gqd ( dd_real *a, gdd_real *b )
{
   b->x = a->x[0];
   b->y = a->x[1];
}

void gqd2qd ( gdd_real *a, dd_real *b )
{
   b->x[0] = a->x;
   b->x[1] = a->y;
}

void qd2gqd ( double *a, double *b )
{
   *b = *a;
}

void gqd2qd ( double *a, double *b )
{
   *b = *a;
}

void comp1_gqd2qd(gd_complex* a, complexH<double>* b)
{
   b->real = a->real;
   b->imag = a->imag;
}

void comp1_gqd2qd(gdd_complex* a, complexH<dd_real>* b)
{
   gqd2qd(&(a->real),&(b->real));
   gqd2qd(&(a->imag),&(b->imag));
}

void comp1_gqd2qd(gqd_complex* a, complexH<qd_real>* b)
{
   gqd2qd(&(a->real),&(b->real));
   gqd2qd(&(a->imag),&(b->imag));
}

void comp1_qd2gqd(complexH<double>* a, gd_complex* b)
{
	b->real = a->real;
	b->imag = a->imag;
}

void comp1_qd2gqd(complexH<dd_real>* a, gdd_complex* b)
{
   qd2gqd(&(a->real),&(b->real));
   qd2gqd(&(a->imag),&(b->imag));
}

void comp1_qd2gqd(complexH<qd_real>* a, gqd_complex* b)
{
   qd2gqd(&(a->real),&(b->real));
   qd2gqd(&(a->imag),&(b->imag));
}

