/* Tests the evaluation and differentiation of a monomial
 * in triple double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_numbers.h"
#include "triple_double_functions.h"
#include "random3_vectors.h"
#include "random3_series.h"
#include "dbl3_monomials_host.h"
#include "dbl3_monomials_kernels.h"

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
   double *cffhi, double *cffmi, double *cfflo );
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
 *   cffhi   space allocated for deg+1 doubles;
 *   cffmi   space allocated for deg+1 doubles;
 *   cfflo   space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffhi   deg+1 doubles with the high doubles of the coefficients;
 *   cffmi   deg+1 doubles with the middle doubles of the coefficients;
 *   cfflo   deg+1 doubles with the low doubles of the coefficients. */

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cffrehi has space allocated for deg+1 doubles,
 *           for the high doubles of the real parts of the coefficients;
 *   cffremi has space allocated for deg+1 doubles,
 *           for the middle doubles of the real parts of the coefficients;
 *   cffrelo has space allocated for deg+1 doubles,
 *           for the low doubles of the real parts of the coefficients;
 *   cffimhi has space allocated for deg+1 doubles,
 *           for the high doubles of the imaginary parts of the coefficients;
 *   cffimmi has space allocated for deg+1 doubles,
 *           for the middle doubles of the imaginary parts of the coefficients;
 *   cffimlo has space allocated for deg+1 doubles,
 *           for the low doubles of the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffrehi holds deg+1 doubles with the high doubles of the real parts
 *           of the coefficients of the power series;
 *   cffremi holds deg+1 doubles with the middle doubles of the real parts
 *           of the coefficients of the power series;
 *   cffrelo holds deg+1 doubles with the low doubles of the real parts
 *           of the coefficients of the power series;
 *   cffimhi holds deg+1 doubles with the high doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimmi holds deg+1 doubles with the middle doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimlo holds deg+1 doubles with the low doubles of the imaginary parts
 *           of the coefficients of the power series. */

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
 ( int dim, int deg, double **datahi, double **datami, double **datalo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   datahi   space allocated for dim arrays of deg+1 doubles;
 *   datami   space allocated for dim arrays of deg+1 doubles;
 *   datalo   space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datahi   datahi[i][j] is the high double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datami   datahi[i][j] is the middle double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalo   datalo[i][j] is the low double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg. */

void make_complex_input
 ( int dim, int deg, double **datarehi, double **dataremi, double **datarelo,
   double **dataimhi, double **dataimmi, double **dataimlo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   datarehi has space allocated for the high doubles of the real parts
 *            of dim series of degree deg.
 *   dataremi has space allocated for the middle doubles of the real parts
 *            of dim series of degree deg.
 *   datarelo has space allocated for the low doubles of the real parts
 *            of dim series of degree deg;
 *   dataimhi has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg;
 *   dataimmi has space allocated for the middle doubles of the imaginary parts
 *            of dim series of degree deg;
 *   dataimlo has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg.
 *
 * ON RETURN :
 *   datarehi contains the high doubles of the real parts of the data,
 *   dataremi contains the middle doubles of the real parts of the data,
 *   datarehi contains the low doubles of the real parts of the data,
 *   dataimhi contains the high doubles of the imaginary parts of the data,
 *   dataimmi contains the middle doubles of the imaginary parts of the data,
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

  // cout << endl << "Testing for real input data ... " << endl;
  // test_real(dim,nvr,pwr,deg);
   cout << endl << "Testing for complex input data ..." << endl;
   test_complex(dim,nvr,pwr,deg);

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
   double *cffhi, double *cffmi, double *cfflo )
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
         random_triple_double(&cffhi[i],&cffmi[i],&cfflo[i]);

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
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo )
{
   bool fail;

   if(nvr > dim)
   {
      cout << "ERROR: nvr = " << nvr << " > " << dim << " dim" << endl;
      return true;
   }
   else
   {
      double rndhi,rndmi,rndlo,sinhi,sinmi,sinlo;

      for(int i=0; i<=deg; i++)
      {
         random_triple_double(&rndhi,&rndmi,&rndlo);    // random cos
                                                        // cos(angle)
         cffrehi[i] = rndhi; cffremi[i] = rndmi;
         cffrelo[i] = rndlo;
         tdf_sqrt(rndhi,rndmi,rndlo,&sinhi,&sinmi,&sinlo);  // cos^(angle)
         tdf_minus(&sinhi,&sinmi,&sinlo);                   // -cos^(angle)
         tdf_inc_d(&sinhi,&sinmi,&sinlo,1.0);               // 1-cos^2(angle)
         // sin is sqrt
         tdf_sqrt(sinhi,sinmi,sinlo,&cffimhi[i],&cffimmi[i],&cffimlo[i]);
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
 ( int dim, int deg, double **datahi, double **datami, double **datalo )
{
   double rndhi,rndmi,rndlo;
   double* pluxhi = new double[deg+1];
   double* pluxmi = new double[deg+1];
   double* pluxlo = new double[deg+1];
   double* minxhi = new double[deg+1];
   double* minxmi = new double[deg+1];
   double* minxlo = new double[deg+1];

   for(int i=0; i<dim; i++)
   {
      random_dbl3_exponentials
         (deg,&rndhi,&rndmi,&rndlo,pluxhi,pluxmi,pluxlo,minxhi,minxmi,minxlo);
      for(int j=0; j<=deg; j++)
      {
         datahi[i][j] = pluxhi[j];
         datami[i][j] = pluxmi[j];
         datalo[i][j] = pluxlo[j];
      }
   }
   free(pluxhi); free(pluxmi); free(pluxlo);
   free(minxhi); free(minxmi); free(minxlo); 

//      for(int j=0; j<=deg; j++)
//         random_triple_double(&datahi[i][j],&datami[i][j],&datalo[i][j]);

}

void make_complex_input
 ( int dim, int deg, double **datarehi, double **dataremi, double **datarelo,
   double **dataimhi, double **dataimmi, double **dataimlo )
{
   double rndrehi,rndremi,rndrelo,rndimhi,rndimmi,rndimlo;
   double* pluxrehi = new double[deg+1];
   double* pluxremi = new double[deg+1];
   double* pluxrelo = new double[deg+1];
   double* pluximhi = new double[deg+1];
   double* pluximmi = new double[deg+1];
   double* pluximlo = new double[deg+1];
   double* minxrehi = new double[deg+1];
   double* minxremi = new double[deg+1];
   double* minxrelo = new double[deg+1];
   double* minximhi = new double[deg+1];
   double* minximmi = new double[deg+1];
   double* minximlo = new double[deg+1];

   for(int i=0; i<dim; i++)
   {
      random_cmplx3_exponentials(deg,
         &rndrehi,&rndremi,&rndrelo,&rndimhi,&rndimmi,&rndimlo,
         pluxrehi,pluxremi,pluxrelo,pluximhi,pluximmi,pluximlo,
         minxrehi,minxremi,minxrelo,minximhi,minximmi,minximlo);

      for(int j=0; j<=deg; j++)
      {
         datarehi[i][j] = pluxrehi[j]; dataremi[i][j] = pluxremi[j];
         datarelo[i][j] = pluxrelo[j];
         dataimhi[i][j] = pluximhi[j]; dataimmi[i][j] = pluximmi[j];
         dataimlo[i][j] = pluximlo[j];
      }
   }
   free(pluxrehi); free(pluxremi); free(pluxrelo);
   free(pluximhi); free(pluximmi); free(pluximlo);
   free(minxrehi); free(minxremi); free(minxrelo);
   free(minximhi); free(minximmi); free(minximlo); 
}

int test_real ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffhi = new double[deg+1]; // high doubles of series coefficient
   double *cffmi = new double[deg+1]; // middle doubles of series coefficient
   double *cfflo = new double[deg+1]; // low doubles of series coefficient

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputhi = new double*[dim];
   for(int i=0; i<dim; i++) inputhi[i] = new double[deg+1];
   double **inputmi = new double*[dim];
   for(int i=0; i<dim; i++) inputmi[i] = new double[deg+1];
   double **inputlo = new double*[dim];
   for(int i=0; i<dim; i++) inputlo[i] = new double[deg+1];
   double **outputhi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhi_h[i] = new double[deg+1];
   double **outputhi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhi_d[i] = new double[deg+1];
   double **outputmi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputmi_h[i] = new double[deg+1];
   double **outputmi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputmi_d[i] = new double[deg+1];
   double **outputlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlo_h[i] = new double[deg+1];
   double **outputlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlo_d[i] = new double[deg+1];

   srand(time(NULL));

   bool fail = make_real_monomial(dim,nvr,pwr,deg,idx,exp,cffhi,cffmi,cfflo);

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
         cout << cffhi[i] << "  " << cffmi[i] << "  " << cfflo[i] << endl;
   }

   make_real_input(dim,deg,inputhi,inputmi,inputlo);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
         cout << inputhi[i][j] << "  " << inputmi[i][j]
                               << "  " << inputlo[i][j] << endl;
   }
   CPU_dbl3_evaldiff(dim,nvr,deg,idx,cffhi,cffmi,cfflo,
                     inputhi,inputmi,inputlo,
                     outputhi_h,outputmi_h,outputlo_h);
   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
      cout << outputhi_h[dim][i] << "  " << outputmi_h[dim][i] 
                                 << "  " << outputlo_h[dim][i] << endl;

   double errsum = 0.0;
   double errtot = 0.0;

   if(nvr > 2)
   {
      GPU_dbl3_evaldiff(deg+1,dim,nvr,deg,idx,cffhi,cffmi,cfflo,inputhi,
                        inputmi,inputlo,outputhi_d,outputmi_d,outputlo_d);
      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputhi_d[dim][i] << "  "
              << outputmi_d[dim][i] << "  "
              << outputlo_d[dim][i] << endl;
         errsum = abs(outputhi_h[dim][i] - outputhi_d[dim][i])
                + abs(outputmi_h[dim][i] - outputmi_d[dim][i])
                + abs(outputlo_h[dim][i] - outputlo_d[dim][i]);
      }
      cout << "Sum of errors : " << errsum << endl; errtot += errsum;
   }

   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputhi_h[idx[k]][i] << "  "
              << outputmi_h[idx[k]][i] << "  "
              << outputlo_h[idx[k]][i] << endl;

      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         for(int i=0; i<=deg; i++)
         {
            cout << outputhi_d[idx[k]][i] << "  "
                 << outputmi_d[idx[k]][i] << "  "
                 << outputlo_d[idx[k]][i] << endl;
            errsum = abs(outputhi_h[idx[k]][i] - outputhi_d[idx[k]][i])
                   + abs(outputmi_h[idx[k]][i] - outputmi_d[idx[k]][i])
                   + abs(outputlo_h[idx[k]][i] - outputlo_d[idx[k]][i]);
         }
         cout << "Sum of errors : " << errsum << endl; errtot += errsum;
      }
      cout << "Total sum of all errors : " << errtot << endl;
   }

   return 0;
}

int test_complex ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];       // indices of variables in the monomial
   int *exp = new int[nvr];       // exponents of the variables
   int *expfac = new int[nvr];    // exponents of common factor
   int nbrfac;                    // number of common factors
   double *cffrehi = new double[deg+1]; // high real parts of coefficients
   double *cffremi = new double[deg+1]; // middle real parts of coefficients
   double *cffrelo = new double[deg+1]; // low real parts of coefficients
   double *cffimhi = new double[deg+1]; // high imaginary coefficients
   double *cffimmi = new double[deg+1]; // middle imaginary coefficients
   double *cffimlo = new double[deg+1]; // low imaginary coefficients

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrehi = new double*[dim];
   double **inputremi = new double*[dim];
   double **inputrelo = new double*[dim];
   double **inputimhi = new double*[dim];
   double **inputimmi = new double*[dim];
   double **inputimlo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputrehi[i] = new double[deg+1]; inputremi[i] = new double[deg+1];
      inputrelo[i] = new double[deg+1];
      inputimhi[i] = new double[deg+1]; inputimmi[i] = new double[deg+1];
      inputimlo[i] = new double[deg+1];
   }
   double **outputrehi_h = new double*[dim+1];
   double **outputremi_h = new double*[dim+1];
   double **outputrelo_h = new double*[dim+1];
   double **outputimhi_h = new double*[dim+1];
   double **outputimmi_h = new double*[dim+1];
   double **outputimlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_h[i] = new double[deg+1];
      outputremi_h[i] = new double[deg+1];
      outputrelo_h[i] = new double[deg+1];
      outputimhi_h[i] = new double[deg+1];
      outputimmi_h[i] = new double[deg+1];
      outputimlo_h[i] = new double[deg+1];
   }
   double **outputrehi_d = new double*[dim+1];
   double **outputremi_d = new double*[dim+1];
   double **outputrelo_d = new double*[dim+1];
   double **outputimhi_d = new double*[dim+1];
   double **outputimmi_d = new double*[dim+1];
   double **outputimlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_d[i] = new double[deg+1];
      outputremi_d[i] = new double[deg+1];
      outputrelo_d[i] = new double[deg+1];
      outputimhi_d[i] = new double[deg+1];
      outputimmi_d[i] = new double[deg+1];
      outputimlo_d[i] = new double[deg+1];
   }

   srand(time(NULL));

   bool fail = make_complex_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo);

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
         cout << cffrehi[i] << "  " << cffremi[i] << "  " << cffrelo[i]
              << endl;
         cout << cffimhi[i] << "  " << cffimmi[i] << "  " << cffimlo[i]
              << endl;
      }
   }
   make_complex_input
      (dim,deg,inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << inputrehi[i][j] << "  " << inputremi[i][j]
                                 << "  " << inputrelo[i][j] << endl;
         cout << inputimhi[i][j] << "  " << inputimmi[i][j]
                                 << "  " << inputimlo[i][j] << endl;
      }
   }
   CPU_cmplx3_evaldiff(dim,nvr,deg,idx,
      cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
      inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
      outputrehi_h,outputremi_h,outputrelo_h,
      outputimhi_h,outputimmi_h,outputimlo_h);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << outputrehi_h[dim][i] << "  " << outputremi_h[dim][i]
                                   << "  " << outputrelo_h[dim][i] << endl;
      cout << outputimhi_h[dim][i] << "  " << outputimmi_h[dim][i]
                                   << "  " << outputimlo_h[dim][i] << endl;
   }

   double errsum = 0.0;
   double errtot = 0.0;

   if(nvr > 2)
   {
      GPU_cmplx3_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
         inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
         outputrehi_d,outputremi_d,outputrelo_d,
         outputimhi_d,outputimmi_d,outputimlo_d);

      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++) 
      {
         cout << outputrehi_d[dim][i] << "  " << outputremi_d[dim][i]
                                      << "  " << outputrelo_d[dim][i] << endl;
         cout << outputimhi_d[dim][i] << "  " << outputimmi_d[dim][i] 
                                      << "  " << outputimlo_d[dim][i] << endl;
         errsum = abs(outputrehi_h[dim][i] - outputrehi_d[dim][i])
                + abs(outputremi_h[dim][i] - outputremi_d[dim][i])
                + abs(outputrelo_h[dim][i] - outputrelo_d[dim][i])
                + abs(outputimhi_h[dim][i] - outputimhi_d[dim][i])
                + abs(outputimmi_h[dim][i] - outputimmi_d[dim][i])
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
              << outputremi_h[idx[k]][i] << "  "
              << outputrelo_h[idx[k]][i] << endl;
         cout << outputimhi_h[idx[k]][i] << "  "
              << outputimmi_h[idx[k]][i] << "  "
              << outputimlo_h[idx[k]][i] << endl;
      }
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         for(int i=0; i<=deg; i++)
         {
            cout << outputrehi_d[idx[k]][i] << "  "
                 << outputremi_d[idx[k]][i] << "  "
                 << outputrelo_d[idx[k]][i] << endl;
            cout << outputimhi_d[idx[k]][i] << "  "
                 << outputimmi_d[idx[k]][i] << "  "
                 << outputimlo_d[idx[k]][i] << endl;
            errsum = abs(outputrehi_h[idx[k]][i] - outputrehi_d[idx[k]][i])
                   + abs(outputremi_h[idx[k]][i] - outputremi_d[idx[k]][i])
                   + abs(outputrelo_h[idx[k]][i] - outputrelo_d[idx[k]][i])
                   + abs(outputimhi_h[idx[k]][i] - outputimhi_d[idx[k]][i])
                   + abs(outputimmi_h[idx[k]][i] - outputimmi_d[idx[k]][i])
                   + abs(outputimlo_h[idx[k]][i] - outputimlo_d[idx[k]][i]);
         }
         cout << "The sum of errors : " << errsum << endl; errtot += errsum;
      }
   }
   cout << "Total sum of all errors : " << errtot << endl;

   return 0;
}
