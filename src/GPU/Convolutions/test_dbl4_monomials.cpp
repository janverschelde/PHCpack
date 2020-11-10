/* Tests the evaluation and differentiation of a monomial
 * in quad double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_numbers.h"
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "random4_series.h"
#include "dbl4_monomials_host.h"
// #include "dbl4_monomials_kernels.h"

using namespace std;

bool sorted_insert ( int n, int *data );
/*
 * DESCRIPTION :
 *   Inserts data[n] in the sequence of n sorted numbers in data.
 *   Returns true if data[n] was already inserted, that is:
 *   there is an index k less than n, for which data[k] == data[n].
 *   Returns false if data[n] is not a duplicate number. */

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cffhihi has space allocated for deg+1 doubles;
 *   cfflohi has space allocated for deg+1 doubles;
 *   cffhilo has space allocated for deg+1 doubles;
 *   cfflolo has space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffhihi stores the deg+1 highest doubles of the coefficients;
 *   cfflohi stores the deg+1 second highest doubles of the coefficients;
 *   cffhilo stores the deg+1 second lowest doubles of the coefficients;
 *   cfflolo stores the deg+1 lowest doubles of the coefficients. */

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi,
   double *cffimhilo, double *cffimlolo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nvr       number of variables with positive power in the monomial;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   exp       space allocated for nvr integers;
 *   cffrehihi has space allocated for deg+1 doubles,
 *             for the highest doubles of the real parts of the coefficients;
 *   cffrelohi has space allocated for deg+1 doubles, for the second
 *             highest doubles of the real parts of the coefficients;
 *   cffrehilo has space allocated for deg+1 doubles, for the second
 *             lowest doubles of the real parts of the coefficients;
 *   cffrelolo has space allocated for deg+1 doubles, for the lowest
 *             doubles of the real parts of the coefficients;
 *   cffimhihi has space allocated for deg+1 doubles, for the highest
 *             doubles of the imaginary parts of the coefficients;
 *   cffimlohi has space allocated for deg+1 doubles, for the second
 *             highest doubles of the imaginary parts of the coefficients;
 *   cffimhilo has space allocated for deg+1 doubles, for the second
 *             lowest doubles of the imaginary parts of the coefficients.
 *   cffimlolo has space allocated for deg+1 doubles, for the lowest
 *             doubles of the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx       nvr integers in the range from 0 to dim-1,
 *             idx(k) is the index of the k-th variable in the monomial;
 *   exp       nvr positive integers with the powers of the variables,
 *             exp(k) is the power of the variable with index idx(k);
 *   cffrehihi holds the deg+1 highest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffrelohi holds the deg+1 second highest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffrehilo holds the deg+1 second lowest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffrelolo holds the deg+1 lowest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffimhihi holds the deg+1 highest doubles of the imaginary parts
 *             of the coefficients of the power series;
 *   cffimlohi holds the deg+1 second highest doubles of the imaginary parts
 *             of the coefficients of the power series;
 *   cffimhilo holds the deg+1 second lowest doubles of the imaginary parts
 *             of the coefficients of the power series;
 *   cffimlolo holds deg+1 lowest doubles of the imaginary parts
 *             of the coefficients of the power series. */

void common_factors ( int nvr, int *exp, int *nbrfac, int *expfac );
/*
 * DESCRIPTION :
 *   Extracts all exponents strictly larger than one.
 *  
 * ON ENTRY :
 *   nvr     number of variables in exp, exp[k] >= 1,
 *           for all k from 0 to nvr-1;
 *   exp     exponents of a monomial;
 *   expfac  space for nvr integers.
 *
 * ON RETURN :
 *   nbrfac  number of exponents in exp strictly larger than one;
 *   expfac  exponents of the common factor,
 *           if exp[k] > 1, then expfac[k] = exp[k]-1.  */

void make_real_input
 ( int dim, int deg, double **datahihi, double **datalohi,
                     double **datahilo, double **datalolo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   datahihi has space allocated for dim arrays of deg+1 doubles;
 *   datalohi has space allocated for dim arrays of deg+1 doubles;
 *   datahilo has space allocated for dim arrays of deg+1 doubles;
 *   datalolo has space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datahihi holds the highest doubles of the input series,
 *            datahihi[i][j] is the highest double of the j-th coefficient 
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalohi holds the second highest doubles of the input series,
 *            datalohi[i][j] is the 2nd highest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datahilo holds the second lowest doubles of the input series,
 *            datalolo[i][j] is the 2nd lowest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalolo holds the lowest doubles of the input series,
 *            datalolo[i][j] is the lowest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg. */

//void make_complex_input
// ( int dim, int deg,
//   double **datarehi, double **datarelo,
//   double **dataimhi, double **dataimlo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   datarehi has space allocated for the high doubles of the real parts
 *            of dim series of degree deg.
 *   datarelo has space allocated for the low doubles of the real parts
 *            of dim series of degree deg;
 *   dataimhi has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg;
 *   dataimlo has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg.
 *
 * ON RETURN :
 *   datarehi contains the high doubles of the real parts of the data,
 *   datarehi contains the low doubles of the real parts of the data,
 *   dataimhi contains the high doubles of the imaginary parts of the data,
 *   dataimlo contains the low doubles of the imaginary parts of the data,
 *            data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

int test_real ( int dim, int nvr, int pwr, int deg );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series. */

int test_complex ( int dim, int nvr, int pwr, int deg );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series. */

int main ( void )
{
   int dim,nvr,pwr,deg;

   cout << "Give the dimension : "; cin >> dim;
   cout << "Give the number of variables, <= "; cout << dim;
   cout << " : "; cin >> nvr;
   cout << "Give the largest power of each variable : "; cin >> pwr;
   cout << "Give the degree of the series : "; cin >> deg;

   cout << endl << "Testing for real input data ... " << endl;
   test_real(dim,nvr,pwr,deg);
   // cout << endl << "Testing for complex input data ..." << endl;
   // test_complex(dim,nvr,pwr,deg);

   return 0;
}

bool sorted_insert ( int n, int *data )
{
   if(n == 0)
      return false;
   else
   {
      int nbr = data[n];
      int idx = n;

      for(int i=0; i<n; i++)
         if(data[i] >= nbr)
         {
            idx = i; break;
         }

      if(idx == n)
         return false;                  // sequence is already sorted
      else
      {
         if(data[idx] == nbr)           // found duplicate number
            return true;
         else
         {
            for(int i=n; i>idx; i--)
               data[i] = data[i-1];     // shift the numbers

            data[idx] = nbr;            // insert number
            return false;
         }
      }
   }
}

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo )
{
   bool fail;

   if(nvr > dim)
   {
      cout << "ERROR: nvr = " << nvr << " > " << dim << " dim" << endl;
      return true;
   }
   else
   {
      for(int i=0; i<=deg; i++)
         random_quad_double(&cffhihi[i],&cfflohi[i],&cffhilo[i],&cfflolo[i]);

      for(int i=0; i<nvr; i++)
      {
         exp[i] = 1 + (rand() % pwr);
         do
         {
            idx[i] = rand() % dim;
            fail = sorted_insert(i,idx);
         }
         while(fail);
      }
      return false;
   }
}

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo )
{
   bool fail;

   if(nvr > dim)
   {
      cout << "ERROR: nvr = " << nvr << " > " << dim << " dim" << endl;
      return true;
   }
   else
   {
      double rndhihi,rndlohi,rndhilo,rndlolo;
      double sinhihi,sinlohi,sinhilo,sinlolo;

      for(int i=0; i<=deg; i++)
      {
         random_quad_double
            (&rndhihi,&rndlohi,&rndhilo,&rndlolo);           // random cos

         cffrehihi[i] = rndhihi; cffrelohi[i] = rndlohi;     // cos(angle)
         cffrehilo[i] = rndhilo; cffrelolo[i] = rndlolo;
         qdf_sqrt(rndhihi,rndlohi,rndhilo,rndlolo,
                  &sinhihi,&sinlohi,&sinhilo,&sinlolo);      // cos^(angle)
         qdf_minus(&sinhihi,&sinlohi,&sinhilo,&sinlolo);     // -cos^(angle)
         qdf_inc_d(&sinhihi,&sinlohi,&sinhilo,&sinlolo,1.0); // 1-cos^2(angle)
         qdf_sqrt(sinhihi,sinlohi,sinhilo,sinlolo,
                  &cffimhihi[i],&cffimlohi[i],
                  &cffimhilo[i],&cffimlolo[i]);              // sin is sqrt
      }
      for(int i=0; i<nvr; i++)
      {
         exp[i] = 1 + (rand() % pwr);
         do
         {
            idx[i] = rand() % dim;
            fail = sorted_insert(i,idx);
         }
         while(fail);
      }
      return false;
   }
}

void common_factors ( int nvr, int *exp, int *nbrfac, int *expfac )
{
   *nbrfac = 0;

   for(int i=0; i<nvr; i++)
   {
      if(exp[i] <= 1)
         expfac[i] = 0;
      else
      {
         expfac[i] = exp[i] - 1;
         *nbrfac = *nbrfac + 1;
      }
   }
}

void make_real_input
 ( int dim, int deg, double **datahihi, double **datalohi,
                     double **datahilo, double **datalolo )
{
   double rndhihi,rndlohi,rndhilo,rndlolo;
   double* pluxhihi = new double[deg+1];
   double* pluxlohi = new double[deg+1];
   double* pluxhilo = new double[deg+1];
   double* pluxlolo = new double[deg+1];
   double* minxhihi = new double[deg+1];
   double* minxlohi = new double[deg+1];
   double* minxhilo = new double[deg+1];
   double* minxlolo = new double[deg+1];

   for(int i=0; i<dim; i++)
   {
      random_dbl4_exponentials
         (deg,&rndhihi,&rndlohi,&rndhilo,&rndlolo,
              pluxhihi,pluxlohi,pluxhilo,pluxlolo,
              minxhihi,minxlohi,minxhilo,minxlolo);

      for(int j=0; j<=deg; j++)
      {
         datahihi[i][j] = pluxhihi[j];
         datalohi[i][j] = pluxlohi[j];
         datalolo[i][j] = pluxhilo[j];
         datalolo[i][j] = pluxlolo[j];
      }
   }
   free(pluxhihi); free(pluxlohi); free(pluxhilo); free(pluxlolo);
   free(minxhihi); free(minxlohi); free(minxhilo); free(minxlolo); 
}

/*
void make_complex_input
 ( int dim, int deg, double **datarehi, double **datarelo,
   double **dataimhi, double **dataimlo )
{
   double rndrehi,rndrelo,rndimhi,rndimlo;
   double* pluxrehi = new double[deg+1];
   double* pluxrelo = new double[deg+1];
   double* pluximhi = new double[deg+1];
   double* pluximlo = new double[deg+1];
   double* minxrehi = new double[deg+1];
   double* minxrelo = new double[deg+1];
   double* minximhi = new double[deg+1];
   double* minximlo = new double[deg+1];

   for(int i=0; i<dim; i++)
   {
      random_cmplx4_exponentials(deg,
         &rndrehi,&rndrelo,&rndimhi,&rndimlo,
         pluxrehi,pluxrelo,pluximhi,pluximlo,
         minxrehi,minxrelo,minximhi,minximlo);

      for(int j=0; j<=deg; j++)
      {
         datarehi[i][j] = pluxrehi[j]; datarelo[i][j] = pluxrelo[j];
         dataimhi[i][j] = pluximhi[j]; dataimlo[i][j] = pluximlo[j];
      }
   }
   free(pluxrehi); free(pluxrelo); free(pluximhi); free(pluximlo);
   free(minxrehi); free(minxrelo); free(minximhi); free(minximlo); 
}
*/

int test_real ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffhihi = new double[deg+1]; // highest coefficient doubles
   double *cfflohi = new double[deg+1]; // second highest coefficient doubles
   double *cffhilo = new double[deg+1]; // second lowest coefficient doubles
   double *cfflolo = new double[deg+1]; // lowest coefficient doubles

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputhihi = new double*[dim];
   for(int i=0; i<dim; i++) inputhihi[i] = new double[deg+1];
   double **inputlohi = new double*[dim];
   for(int i=0; i<dim; i++) inputlohi[i] = new double[deg+1];
   double **inputhilo = new double*[dim];
   for(int i=0; i<dim; i++) inputhilo[i] = new double[deg+1];
   double **inputlolo = new double*[dim];
   for(int i=0; i<dim; i++) inputlolo[i] = new double[deg+1];
   double **outputhihi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihi_h[i] = new double[deg+1];
   double **outputlohi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohi_h[i] = new double[deg+1];
   double **outputhilo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilo_h[i] = new double[deg+1];
   double **outputlolo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlolo_h[i] = new double[deg+1];
   double **outputhihi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihi_d[i] = new double[deg+1];
   double **outputlohi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohi_d[i] = new double[deg+1];
   double **outputhilo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilo_d[i] = new double[deg+1];
   double **outputlolo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlolo_d[i] = new double[deg+1];

   srand(time(NULL));

   bool fail = make_real_monomial(dim,nvr,pwr,deg,idx,exp,
                                  cffhihi,cfflohi,cffhilo,cfflolo);

   if(!fail)
   {
      cout << "Generated a random monomial :" << endl;
      cout << "   the indices :";
      for(int i=0; i<nvr; i++) cout << " " << idx[i];
      cout << endl;

      cout << " the exponents :";
      for(int i=0; i<nvr; i++) cout << " " << exp[i];
      cout << endl;

      common_factors(nvr,exp,&nbrfac,expfac);

      cout << "common factors :";
      for(int i=0; i<nvr; i++) cout << " " << expfac[i];
      cout << endl;
      cout << "number of common factors : " << nbrfac << endl;

      cout << scientific << setprecision(16);
      cout << "the coefficients :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << " " << cffhihi[i] << " " << cfflohi[i] << endl;;
         cout << " " << cffhilo[i] << " " << cfflolo[i] << endl;;
      }
   }
   make_real_input(dim,deg,inputhihi,inputlohi,inputhilo,inputlolo);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl;
         cout << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
      }
   }
   CPU_dbl4_evaldiff(dim,nvr,deg,idx,cffhihi,cfflohi,cffhilo,cfflolo,
                     inputhihi,inputlohi,inputhilo,inputlolo,
                     outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << outputhihi_h[dim][i] << "  " << outputlohi_h[dim][i] << endl;
      cout << outputhilo_h[dim][i] << "  " << outputlolo_h[dim][i] << endl;
   }

//   double errsum = 0.0;
//   double errtot = 0.0;
 
   if(nvr > 2)
   {
/*
      GPU_dbl4_evaldiff(deg+1,dim,nvr,deg,idx,cffhi,cfflo,inputhi,inputlo,
                        outputhi_d,outputlo_d);
      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputhi_d[dim][i] << "  " << outputlo_d[dim][i] << endl;
         errsum = abs(outputhi_h[dim][i] - outputhi_d[dim][i])
                + abs(outputlo_h[dim][i] - outputlo_d[dim][i]);
      }
      cout << "Sum of errors : " << errsum << endl; errtot += errsum;
 */
   }
   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputhihi_h[idx[k]][i] << "  "
              << outputlohi_h[idx[k]][i] << endl;
         cout << outputhilo_h[idx[k]][i] << "  "
              << outputlolo_h[idx[k]][i] << endl;
      }
/*
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         for(int i=0; i<=deg; i++)
         {
            cout << outputhi_d[idx[k]][i] << "  "
                 << outputlo_d[idx[k]][i] << endl;
            errsum = abs(outputhi_h[idx[k]][i] - outputhi_d[idx[k]][i])
                   + abs(outputlo_h[idx[k]][i] - outputlo_d[idx[k]][i]);
            cout << "Sum of errors : " << errsum << endl; errtot += errsum;
         }
      }
   }
   cout << "Total sum of all errors : " << errtot << endl;
 */
   }
   return 0;
}

/*

int test_complex ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];             // indices of variables in the monomial
   int *exp = new int[nvr];             // exponents of the variables
   int *expfac = new int[nvr];          // exponents of common factor
   int nbrfac;                          // number of common factors
   double *cffrehi = new double[deg+1]; // high real parts of coefficients
   double *cffrelo = new double[deg+1]; // high real parts of coefficients
   double *cffimhi = new double[deg+1]; // low imaginary parts of coefficients
   double *cffimlo = new double[deg+1]; // low imaginary parts of coefficients

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrehi = new double*[dim];
   double **inputrelo = new double*[dim];
   double **inputimhi = new double*[dim];
   double **inputimlo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputrehi[i] = new double[deg+1];
      inputrelo[i] = new double[deg+1];
      inputimhi[i] = new double[deg+1];
      inputimlo[i] = new double[deg+1];
   }
   double **outputrehi_h = new double*[dim+1];
   double **outputrelo_h = new double*[dim+1];
   double **outputimhi_h = new double*[dim+1];
   double **outputimlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_h[i] = new double[deg+1];
      outputrelo_h[i] = new double[deg+1];
      outputimhi_h[i] = new double[deg+1];
      outputimlo_h[i] = new double[deg+1];
   }
   double **outputrehi_d = new double*[dim+1];
   double **outputrelo_d = new double*[dim+1];
   double **outputimhi_d = new double*[dim+1];
   double **outputimlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_d[i] = new double[deg+1];
      outputrelo_d[i] = new double[deg+1];
      outputimhi_d[i] = new double[deg+1];
      outputimlo_d[i] = new double[deg+1];
   }
   srand(time(NULL));

   bool fail = make_complex_monomial
     (dim,nvr,pwr,deg,idx,exp,cffrehi,cffrelo,cffimhi,cffimlo);

   if(!fail)
   {
      cout << "Generated a random monomial :" << endl;
      cout << "   the indices :";
      for(int i=0; i<nvr; i++) cout << " " << idx[i];
      cout << endl;

      cout << " the exponents :";
      for(int i=0; i<nvr; i++) cout << " " << exp[i];
      cout << endl;

      common_factors(nvr,exp,&nbrfac,expfac);

      cout << "common factors :";
      for(int i=0; i<nvr; i++) cout << " " << expfac[i];
      cout << endl;
      cout << "number of common factors : " << nbrfac << endl;

      cout << scientific << setprecision(16);
      cout << "the coefficients :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << cffrehi[i] << "  " << cffrelo[i] << endl;
         cout << cffimhi[i] << "  " << cffimlo[i] << endl;
      }
   }
   make_complex_input(dim,deg,inputrehi,inputrelo,inputimhi,inputimlo);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl;
         cout << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   CPU_cmplx4_evaldiff
      (dim,nvr,deg,idx,cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi,inputrelo,inputimhi,inputimlo,
       outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << outputrehi_h[dim][i] << "  " << outputrelo_h[dim][i] << endl;
      cout << outputimhi_h[dim][i] << "  " << outputimlo_h[dim][i] << endl;
   }

   double errsum = 0.0;
   double errtot = 0.0;

   if(nvr > 2)
   {
      GPU_cmplx4_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrehi,cffrelo,cffimhi,cffimlo,
         inputrehi,inputrelo,inputimhi,inputimlo,
         outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d);

      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++) 
      {
         cout << outputrehi_d[dim][i] << "  " << outputrelo_d[dim][i] << endl;
         cout << outputimhi_d[dim][i] << "  " << outputimlo_d[dim][i] << endl;
         errsum = abs(outputrehi_h[dim][i] - outputrehi_d[dim][i])
                + abs(outputrelo_h[dim][i] - outputrelo_d[dim][i])
                + abs(outputimhi_h[dim][i] - outputimhi_d[dim][i])
                + abs(outputimlo_h[dim][i] - outputimlo_d[dim][i]);
      }
      cout << "The sum of errors : " << errsum << endl; errtot += errsum;
   }

   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrehi_h[idx[k]][i] << "  "
              << outputrelo_h[idx[k]][i] << endl;
         cout << outputimhi_h[idx[k]][i] << "  "
              << outputimlo_h[idx[k]][i] << endl;
      }
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         for(int i=0; i<=deg; i++)
         {
            cout << outputrehi_d[idx[k]][i] << "  "
                 << outputrelo_d[idx[k]][i] << endl;
            cout << outputimhi_d[idx[k]][i] << "  "
                 << outputimlo_d[idx[k]][i] << endl;
            errsum = abs(outputrehi_h[idx[k]][i] - outputrehi_d[idx[k]][i])
                   + abs(outputrelo_h[idx[k]][i] - outputrelo_d[idx[k]][i])
                   + abs(outputimhi_h[idx[k]][i] - outputimhi_d[idx[k]][i])
                   + abs(outputimlo_h[idx[k]][i] - outputimlo_d[idx[k]][i]);
         }
         cout << "The sum of errors : " << errsum << endl; errtot += errsum;
      }
   }
   cout << "Total sum of all errors : " << errtot << endl;

   return 0;
}
*/
