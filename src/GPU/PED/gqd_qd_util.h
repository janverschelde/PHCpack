// utilities for the polynomial evaluation and differentiation

#ifndef UTIL
#define UTIL

#include <qd/qd_real.h>
#include "complexD.h"
#include "complexH.h"
#include "cstdlib"

void qd2gqd ( qd_real *a, gqd_real *b );
void gqd2qd ( gqd_real *a, qd_real *b );
void qd2gqd ( dd_real *a, gdd_real *b );
void gqd2qd ( gdd_real *a, dd_real *b );
void qd2gqd ( double *a, double *b );
void gqd2qd ( double *a, double *b );

template<class realD, class realH>
void comp1_gqd2qd ( complexD<realD>* a, complexH<realH>* b )
{
   gqd2qd(&(a->re),&(b->re));
   gqd2qd(&(a->im),&(b->im));
}

template<class realD, class realH>
void comp1_gqdArr2qdArr ( complexD<realD>* a, complexH<realH>* b, int dim )
{
   int i;

   for(i=0;i<dim;i++) comp1_gqd2qd(a+i,b+i);
}

int parse_arguments
 ( int argc, char *argv[],
   int *BS, int *dim, int *NM, int *NV, int *deg, int *r, int *mode );
/*
 * Parses the command line arguments.
 * Returns 0 if okay, otherwise the function returns 1.
 *
 * ON ENTRY :
 *   argc     the number of command line arguments;
 *   argv     the command line arguments.
 *
 * ON RETURN :
 *   BS       number of threads in a block, the block size;
 *   dim      dimension of the problem;
 *   NM       number of monomials;
 *   NV       number of variables;
 *   deg      highest degree of the variables;
 *   r        frequency of the runs;
 *   mode     mode of execution, 0, 1, or 2. */

int parse_args
 ( int argc, char *argv[], int *prc,
   int *BS, int *dim, int *NM, int *NV, int *deg, int *r, int *mode );
/*
 * Parses the command line arguments.
 * Returns 0 if okay, otherwise the function returns 1.
 *
 * ON ENTRY :
 *   argc     the number of command line arguments;
 *   argv     the command line arguments.
 *
 * ON RETURN :
 *   prc      the working precision, 0 for d, 1 for dd, 2 for qd;
 *   BS       number of threads in a block, the block size;
 *   dim      dimension of the problem;
 *   NM       number of monomials;
 *   NV       number of variables;
 *   deg      highest degree of the variables;
 *   r        frequency of the runs;
 *   mode     mode of execution, 0, 1, or 2. */

void random_positions ( int dim, int NM, int NV, int *p_int, char *p_char );
/*
 * Generates random positions for the NM monomials in NV variables, for
 * a system of dimension dim, stored twice in the arrays p_int and p_char. */

void random_exponents ( int NM, int NV, int d, int *e_int, char *e_char );
/*
 * Generates random exponents for the NM monomials in NV variables,
 * of degree at most d, stored twice in the arrays e_int and e_char. */

double random_double ( void ) // Returns a random double in [0,1].
{
   return rand()/((double) RAND_MAX);
}

template<class realD, class realH>
void random_point ( int dim, complexD<realD> *x_d, complexH<realH> *x_h )
/*
 * Generates a random complex point of length dim,
 * stored twice in the arrays x_d and x_h. */
{
   for(int i=0; i<dim; i++)
   {
      double temp = random_double()*2*M_PI;
      x_d[i].initH(cos(temp),sin(temp));
      x_h[i].init(cos(temp),sin(temp));
   }
}

template<class realD, class realH>
void random_coefficients ( int n, complexD<realD> *c_h, complexH<realH> *c )
/*
 * Generates n random coefficients
 * stored twice in the arrays c_h and c. */
{
   for(int i=0; i<n; i++)
   {
      double temp = random_double()*2*M_PI;
      c_h[i].initH(cos(temp),sin(temp));
      c[i].init(cos(temp),sin(temp));
   }
}

template<class realD, class realH>
void generate_system
 ( int dim, int NM, int NV, int deg, int *p_int, char *p_char,
   int *e_int, char *e_char, complexD<realD> *c_d, complexH<realH> *c_h )
/*
 * Generates a random system of dimension dim.
 * On input are the number of monomials NM, number of variables NV,
 * largest degree deg.  On returns are positions in p_int, p_char,
 * exponents e_int, e_char, and coefficients in c_d, c_h. */
{
   random_positions(dim,NM,NV,p_int,p_char);
   random_exponents(NM,NV,deg,e_int,e_char);
   random_coefficients(NM,c_d,c_h);
}

template<class realD, class realH>
void error_on_factors
 ( int NM, complexD<realD> *factors_h, complexH<realH> *factors_s )
/*
 * Compares the factors computes on the host with those computed by
 * the CPU, for a number of monomials equal to NM. */
{
   complexH<realH> factors_error(0.0,0.0);   
   complexH<realH> temp;

   for(int i=0; i<NM; i++)
   {
      comp1_gqd2qd(&factors_h[i], &temp);
      factors_error=factors_error+(factors_s[i]-temp)*(factors_s[i]-temp);
      cout << temp;
      cout << factors_s[i]; 
   }
   cout << "factors_error=" << factors_error; 
}

template<class realD, class realH>
void error_on_derivatives
 ( int dim, complexH<realH> **derivatives, complexD<realD> *values_h )
/*
 * Compares the derivatives of the system of dimension dim
 * computed on the GPU with those computed by the CPU. */
{
   complexH<realH> poly_error;
   complexH<realH> temp;

   for(int i=0;i<dim;i++)
      for(int j=0; j<dim; j++)
      {
         comp1_gqd2qd(&values_h[dim+dim*i+j], &temp);
         poly_error=poly_error
            +(derivatives[j][i]-temp)*(derivatives[j][i]-temp);
      }
   cout << "poly_error=" << poly_error;
}

template<class realH>
void write_system
 ( int dim, int NM, int NV, complexH<realH> *c, int *myp, int *e )
/*
 * Writes the system of dimension dim, with number of monomials NM
 * number of variables NV with coefficients in c, positions of nonzero
 * exponents in p, and values of the exponents in e. */
{
   cout << "          dimension : " << dim << endl;
   cout << "number of monomials : " << NM << endl;
   cout << "number of variables : " << NV << endl;
   cout << "   the coefficients : " << endl;
   cout << scientific << setprecision(8);
   for(int i=0; i<NM; i++)
   {
      cout << "c[" << i << "] : ";
      cout << "(" <<  c[i] << ")";
      for(int j=0; j<NV; j++)
         cout << "*x[" << myp[NV*i+j] << "]^" << e[NV*i+j];
      cout << endl;
   }
}

#endif
