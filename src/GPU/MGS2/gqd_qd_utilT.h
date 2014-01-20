// Conversions between complex templated multiprecision CPU and GPU types
// Generating complex random vectors of given precision 

#ifndef __UTILT__
#define __UTILT__

#include <qd/qd_real.h>
#include <cstdlib>
#include "complex.h"
#include "complexH.h"
#include "gqd_qd_util.h"

template<class T, class T1>
void comp1_gqd2qd(complex<T>* a, complexH<T1>* b)
{
   gqd2qd(&(a->real),&(b->real));
   gqd2qd(&(a->imag),&(b->imag));
}

template<class T, class T1>
void comp1_qd2gqd(complexH<T1>* a, complex<T>* b)
{
   qd2gqd(&(a->real),&(b->real));
   qd2gqd(&(a->imag),&(b->imag));
}

template<class T, class T1>
void comp1_gqdArr2qdArr(complex<T>* a, complexH<T1>* b, int dim)
{
   int i;

   for(i=0; i<dim; i++) comp1_gqd2qd(a+i,b+i);
}

template<class T, class T1>
void random_point ( int dim, complex<T> *x_h, complexH<T1> *x )
// Generates a random complex point of length dim,
// stored twice in the arrays x and x_h.
{
   for(int i=0; i<dim; i++)
   {
      double temp = random_double()*2*M_PI;
      x_h[i].initH(cos(temp),sin(temp));
      x[i].init(cos(temp),sin(temp));
   }
}

#endif
