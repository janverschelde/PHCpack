#ifndef UTIL
#define UTIL

#include<qd/qd_real.h>
#include"complexD.h"
#include"complexH.h"
#include"cstdlib"

void qd2gqd( qd_real *a, gqd_real *b );
void gqd2qd( gqd_real *a, qd_real *b );
void qd2gqd ( dd_real *a, gdd_real *b );
void gqd2qd ( gdd_real *a, dd_real *b );
void qd2gqd ( double *a, double *b );
void gqd2qd ( double *a, double *b );

template<class T, class T1>
void comp1_gqd2qd(complexD<T>* a, complexH<T1>* b)
{
   gqd2qd(&(a->real),&(b->real));
   gqd2qd(&(a->imag),&(b->imag));
}

template<class T, class T1>
void comp1_gqdArr2qdArr(complexD<T>* a, complexH<T1>* b, int dim)
{
   int i;

   for(i=0;i<dim;i++) comp1_gqd2qd(a+i,b+i);
}

int parse_arguments
 ( int argc, char *argv[],
   int *BS, int *dim, int *NM, int *NV, int *deg, int *r, int *mode );

void random_positions ( int dim, int NM, int NV, int *p_int, char *p_char );
// Generates random positions for the NM monomials in NV variables,
// for a system of dimension dim, stored twice in the arrays p_int and p_char.

void random_exponents ( int NM, int NV, int d, int *e_int, char *e_char );
// Generates random exponents for the NM monomials in NV variables,
// of degree at most d, stored twice in the arrays e_int and e_char.

double random_double ( void )
// Returns a random double in [0,1].
{
   return rand()/((double) RAND_MAX);
}

template<class T, class T1>
void random_point ( int dim, complexD<T> *x_h, complexH<T1> *x )
// Generates a random complex point of length dim,
// stored twice in the arrays x and x_h.
{
/*   for(int i=0; i<dim; i++)
   {
      double temp = random_double()*2*M_PI;
      x[i].re = cos(temp);   x[i].im = sin(temp);
      x_h[i].re = cos(temp); x_h[i].im = sin(temp);
   } */
   for (int i=0;i<dim;i++)
     {
        double temp = random_double()*2*M_PI;
        // cout << "cos=" << cos(temp) << endl;
        // cout << "sin=" << sin(temp) << endl;
        x_h[i].initH(cos(temp),sin(temp));
        x[i].init(cos(temp),sin(temp));
        // cout << "cos=" << cos(temp) << endl;
        // cout << "sin=" << sin(temp) << endl;
        // cout << x[i];
     }
}

template<class T, class T1>
void random_coefficients ( int n, complexD<T> *c_h, complexH<T1> *c )
// Generates n random coefficients
// stored twice in the arrays c_h and c.
{
   for(int i=0; i<n; i++)
   {
      double temp = random_double()*2*M_PI;
      c_h[i].initH(cos(temp),sin(temp));
      c[i].init(cos(temp),sin(temp));
      //c_h[i].re = cos(temp); c_h[i].im = sin(temp);
      //c[i].re = cos(temp);   c[i].im = sin(temp);
   }
}

template<class T, class T1>
void generate_system
 ( int dim, int NM, int NV, int deg, int *p_int, char *p_char,
   int *e_int, char *e_char, complexD<T> *c_h, complexH<T1> *c )
// Generates a random system of dimension dim.
// On input are the number of monomials NM, number of variables NV,
// largest degree deg.  On returns are positions in p_int, p_char,
// exponents e_int, e_char, and coefficients in c_h, c.
{
   random_positions(dim,NM,NV,p_int,p_char);
   random_exponents(NM,NV,deg,e_int,e_char);
   random_coefficients(NM,c_h,c);
}

template<class T, class T1>
void error_on_factors ( int NM, complexD<T> *factors_h, complexH<T1> *factors_s )
// Compares the factors computes on the host with those computed by
// the CPU, for a number of monomials equal to NM.
{

   complexH<T1> factors_error(0.0,0.0);   
   complexH<T1> temp;

   for(int i=0;i<NM;i++)
      {
        comp1_gqd2qd(&factors_h[i], &temp);
        factors_error=factors_error+(factors_s[i]-temp)*(factors_s[i]-temp);
        cout << temp;
        cout << factors_s[i]; 
      }

  // cout << scientific << setprecision(sizeof(T)*2);
   cout << "factors_error=" << factors_error; 
/*
   float factors_error_inf_re = 0.0;
   float factors_error_inf_im = 0.0;
   int factors_re_err_ind_i = 0;
   int factors_im_err_ind_i = 0;

   for(int i=0; i<NM; i++)
   {
      if(abs(factors_h[i].re-factors_s[i].re) > factors_error_inf_re)
      {
         factors_error_inf_re = abs(factors_h[i].re-factors_s[i].re);
         factors_re_err_ind_i = i;
      }
      if(abs(factors_h[i].im-factors_s[i].im) > factors_error_inf_im)
      {
         factors_error_inf_im = abs(factors_h[i].im-factors_s[i].im);
         factors_im_err_ind_i = i;
      }
   }
   cout << scientific << setprecision(8);
   cout << " factors_error_re = " << factors_error_inf_re << " "
        << factors_re_err_ind_i << endl;
   cout << " factors_error_im = " << factors_error_inf_im << " "
        << factors_im_err_ind_i << endl;

*/
}

template<class T, class T1>
void error_on_derivatives
 ( int dim, complexH<T1> **derivatives, complexD<T> *values_h )
// compares the derivatives of the system of dimension dim
// computed on the GPU with those computed by the CPU
{

   complexH<T1> poly_error;
   complexH<T1> temp;

   for(int i=0;i<dim;i++)
        for(int j=0; j<dim; j++)

      {
        comp1_gqd2qd(&values_h[dim+dim*i+j], &temp);
        poly_error=poly_error+(derivatives[j][i]-temp)*(derivatives[j][i]-temp);
      }

   cout << "poly_error=" << poly_error;
}

#endif
