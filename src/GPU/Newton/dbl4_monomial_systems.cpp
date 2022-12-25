// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "random4_series.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomial_systems.h"

using namespace std;

void make_real4_exponentials
 ( int dim, int  deg,
   double **shihi, double **slohi, double **shilo, double **slolo )
{
   double rndhihi,rndlohi,rndhilo,rndlolo;
   double acchihi,acclohi,acchilo,acclolo;

   const double fac = 64.0;      // 1024.0;
   const double inc = 63.0/64.0; // 1023.0/1024.0;

   for(int i=0; i<dim; i++)
   {
      // rnd is in [-1, +1]
      random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo); 
      qdf_div(rndhihi,rndlohi,rndhilo,rndlolo,fac,0.0,0.0,0.0,
              &acchihi,&acclohi,&acchilo,&acclolo); // acc in [-1/fac, +1/fac]
      
      if(rndhihi < 0) // if -1/fac <= rnd       < 0
      {               // then   -1 <= rnd - inc < -inc
         qdf_sub(acchihi,acclohi,acchilo,acclolo,inc,0.0,0.0,0.0,
                 &rndhihi,&rndlohi,&rndhilo,&rndlolo);
      }
      else            // if    0  <= rnd       <= 1/fac
      {               // then inc <= rnd + inc <= 1
         qdf_add(acchihi,acclohi,acchilo,acclolo,inc,0.0,0.0,0.0,
                 &rndhihi,&rndlohi,&rndhilo,&rndlolo);
      }
      dbl4_exponential(deg,rndhihi,rndlohi,rndhilo,rndlolo,
                       shihi[i],slohi[i],shilo[i],slolo[i]);
   }
}

void make_complex4_exponentials
 ( int dim, int deg,
   double **srehihi, double **srelohi, double **srehilo, double **srelolo,
   double **simhihi, double **simlohi, double **simhilo, double **simlolo )
{
   double xrehihi,xrelohi,xrehilo,xrelolo;
   double ximhihi,ximlohi,ximhilo,ximlolo;
   double yhihi,ylohi,yhilo,ylolo;

   for(int i=0; i<dim; i++)
   {                                            // cosine of some angle
      random_quad_double(&xrehihi,&xrelohi,&xrehilo,&xrelolo);
 
      qdf_sqr(xrehihi,xrelohi,xrehilo,xrelolo,
              &yhihi,&ylohi,&yhilo,&ylolo);        // y = cos^2
      qdf_minus(&yhihi,&ylohi,&yhilo,&ylolo);      // y = -cos^2
      qdf_inc_d(&yhihi,&ylohi,&yhilo,&ylolo,1.0);  // y = 1 - cos^2
      qdf_sqrt(yhihi,ylohi,yhilo,ylolo,
               &ximhihi,&ximlohi,&ximhilo,&ximlolo); // sin is sqrt(1-cos^2)

      cmplx4_exponential
         (deg,xrehihi,xrelohi,xrehilo,xrelolo,
              ximhihi,ximlohi,ximhilo,ximlolo,
          srehihi[i],srelohi[i],srehilo[i],srelolo[i],
          simhihi[i],simlohi[i],simhilo[i],simlolo[i]);
   }
}

void evaluate_real4_monomials
 ( int dim, int deg, int **rowsA,
   double **shihi, double **slohi, double **shilo, double **slolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo )
{
   const int degp1 = deg+1;

   double *acchihi = new double[degp1]; // accumulates product
   double *acclohi = new double[degp1];
   double *acchilo = new double[degp1];
   double *acclolo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_dbl4_product
                  (deg,shihi[j],slohi[j],shilo[j],slolo[j],
                   rhshihi[i],rhslohi[i],rhshilo[i],rhslolo[i],
                   acchihi,acclohi,acchilo,acclolo);

               for(int L=0; L<degp1; L++)
               {
                  rhshihi[i][L] = acchihi[L];
                  rhslohi[i][L] = acclohi[L];
                  rhshilo[i][L] = acchilo[L];
                  rhslolo[i][L] = acclolo[L];
               }
            }
         }
      }
   }
   free(acchihi); free(acclohi); free(acchilo); free(acclolo);
}

void evaluate_complex4_monomials
 ( int dim, int deg, int **rowsA,
   double **srehihi, double **srelohi, double **srehilo, double **srelolo,
   double **simhihi, double **simlohi, double **simhilo, double **simlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo )
{
   const int degp1 = deg+1;

   double *accrehihi = new double[degp1]; // accumulates product
   double *accrelohi = new double[degp1];
   double *accrehilo = new double[degp1];
   double *accrelolo = new double[degp1];
   double *accimhihi = new double[degp1];
   double *accimlohi = new double[degp1];
   double *accimhilo = new double[degp1];
   double *accimlolo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx4_product
                  (deg,srehihi[j],srelohi[j],srehilo[j],srelolo[j],
                       simhihi[j],simlohi[j],simhilo[j],simlolo[j],
                   rhsrehihi[i],rhsrelohi[i],rhsrehilo[i],rhsrelolo[i],
                   rhsimhihi[i],rhsimlohi[i],rhsimhilo[i],rhsimlolo[i],
                   accrehihi,accrelohi,accrehilo,accrelolo,
                   accimhihi,accimlohi,accimhilo,accimlolo);

               for(int L=0; L<degp1; L++)
               {
                  rhsrehihi[i][L] = accrehihi[L];
                  rhsrelohi[i][L] = accrelohi[L];
                  rhsrehilo[i][L] = accrehilo[L];
                  rhsrelolo[i][L] = accrelolo[L];
                  rhsimhihi[i][L] = accimhihi[L];
                  rhsimlohi[i][L] = accimlohi[L];
                  rhsimhilo[i][L] = accimhilo[L];
                  rhsimlolo[i][L] = accimlolo[L];
               }
            }
         }
      }
   }
   free(accrehihi); free(accrelohi); free(accrehilo); free(accrelolo);
   free(accimhihi); free(accimlohi); free(accimhilo); free(accimlolo);
}

void make_real4_coefficients
 ( int nbrcol, int dim,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo )
{
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         random_quad_double
            (&cffhihi[i][j][0],&cfflohi[i][j][0],
             &cffhilo[i][j][0],&cfflolo[i][j][0]);
}

void make_complex4_coefficients
 ( int nbrcol, int dim,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo )
{
   double xrehihi,xrelohi,ximhihi,ximlohi,yhihi,ylohi;
   double xrehilo,xrelolo,ximhilo,ximlolo,yhilo,ylolo;

   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
      {
         random_quad_double(&xrehihi,&xrelohi,&xrehilo,&xrelolo);
         // cosine of some angle

         qdf_sqr(xrehihi,xrelohi,xrehilo,xrelolo,
                 &yhihi,&ylohi,&yhilo,&ylolo);        // y = cos^2
         qdf_minus(&yhihi,&ylohi,&yhilo,&ylolo);      // y = -cos^2
         qdf_inc_d(&yhihi,&ylohi,&yhilo,&ylolo,1.0);  // y = 1 - cos^2
         qdf_sqrt(yhihi,ylohi,yhilo,ylolo,
                  &ximhihi,&ximlohi,&ximhilo,&ximlolo); // sin is sqrt(1-cos^2)

         cffrehihi[i][j][0] = xrehihi;
         cffrelohi[i][j][0] = xrelohi;
         cffrehilo[i][j][0] = xrehilo;
         cffrelolo[i][j][0] = xrelolo;
         cffimhihi[i][j][0] = ximhihi;
         cffimlohi[i][j][0] = ximlohi;
         cffimhilo[i][j][0] = ximhilo;
         cffimlolo[i][j][0] = ximlolo;
      }
}

void evaluate_real4_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   int vrblvl )
{
   const int degp1 = deg+1;

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "evaluating at the series x ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "x[" << i << "][" << j << "] : "
                 << xhihi[i][j] << "  " << xlohi[i][j] << endl << "  "
                 << xhilo[i][j] << "  " << xlolo[i][j] << endl;
   }
   double **prdrhshihi = new double*[dim];
   double **prdrhslohi = new double*[dim];
   double **prdrhshilo = new double*[dim];
   double **prdrhslolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhshihi[i] = new double[degp1];
      prdrhslohi[i] = new double[degp1];
      prdrhshilo[i] = new double[degp1];
      prdrhslolo[i] = new double[degp1];

      for(int k=0; k<degp1; k++) // initialize sum to 0
      {
         rhshihi[i][k] = 0.0; rhslohi[i][k] = 0.0;
         rhshilo[i][k] = 0.0; rhslolo[i][k] = 0.0;
      }
   }
   rhshihi[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++)   // initialize product to coefficient
      {
         prdrhshihi[i][0] = cffhihi[col][i][0];
         prdrhslohi[i][0] = cfflohi[col][i][0];
         prdrhshilo[i][0] = cffhilo[col][i][0];
         prdrhslolo[i][0] = cfflolo[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhshihi[i][k] = 0.0; prdrhslohi[i][k] = 0.0;
            prdrhshilo[i][k] = 0.0; prdrhslolo[i][k] = 0.0;
         }
      }
      if(vrblvl > 1)
         cout << "Evaluating at column " << col << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         for(int k=0; k<dim; k++) rowsA[i][k] = 0;
         for(int k=0; k<nvr[col][i]; k++) rowsA[i][idx[col][i][k]] = 1;
         if(vrblvl > 1)
         {
            for(int k=0; k<dim; k++) cout << " " << rowsA[i][k];
            cout << endl;
         }
      }
      evaluate_real4_monomials
         (dim,deg,rowsA,xhihi,xlohi,xhilo,xlolo,
          prdrhshihi,prdrhslohi,prdrhshilo,prdrhslolo);

      if(vrblvl > 1)
      {
         cout << scientific << setprecision(16);
         for(int i=0; i<dim; i++)
         {
            cout << "value at dimension " << i << " :" << endl;
            for(int j=0; j<degp1; j++)
               cout << "prdrhs[" << i << "][" << j << "] : "
                    << prdrhshihi[i][j] << "  "
                    << prdrhslohi[i][j] << endl << "  "
                    << prdrhshilo[i][j] << "  "
                    << prdrhslolo[i][j] << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               // rhs[i][k] += prdrhs[i][k];
               qdf_inc(&rhshihi[i][k],&rhslohi[i][k],
                       &rhshilo[i][k],&rhslolo[i][k],
                       prdrhshihi[i][k],prdrhslohi[i][k],
                       prdrhshilo[i][k],prdrhslolo[i][k]);
            }
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "the evaluated series ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << rhshihi[i][j] << "  " << rhslohi[i][j] << endl << "  "
                 << rhshilo[i][j] << "  " << rhslolo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhshihi[i]); free(prdrhslohi[i]);
      free(prdrhshilo[i]); free(prdrhslolo[i]);
   }
   free(prdrhshihi); free(prdrhslohi);
   free(prdrhshilo); free(prdrhslolo);
}

void evaluate_complex4_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **xrehihi, double **xrelohi, double **xrehilo, double **xrelolo,
   double **ximhihi, double **ximlohi, double **ximhilo, double **ximlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo, int vrblvl )
{
   const int degp1 = deg+1;

   double **prdrhsrehihi = new double*[dim];
   double **prdrhsrelohi = new double*[dim];
   double **prdrhsrehilo = new double*[dim];
   double **prdrhsrelolo = new double*[dim];
   double **prdrhsimhihi = new double*[dim];
   double **prdrhsimlohi = new double*[dim];
   double **prdrhsimhilo = new double*[dim];
   double **prdrhsimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhsrehihi[i] = new double[degp1];
      prdrhsrelohi[i] = new double[degp1];
      prdrhsrehilo[i] = new double[degp1];
      prdrhsrelolo[i] = new double[degp1];
      prdrhsimhihi[i] = new double[degp1];
      prdrhsimlohi[i] = new double[degp1];
      prdrhsimhilo[i] = new double[degp1];
      prdrhsimlolo[i] = new double[degp1];

      for(int k=0; k<degp1; k++)  // initialize sum to zero
      {
         rhsrehihi[i][k] = 0.0; rhsimhihi[i][k] = 0.0;
         rhsrelohi[i][k] = 0.0; rhsimlohi[i][k] = 0.0;
         rhsrehilo[i][k] = 0.0; rhsimhilo[i][k] = 0.0;
         rhsrelolo[i][k] = 0.0; rhsimlolo[i][k] = 0.0;
      }
   }
   rhsrehihi[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++)     // initialize product to coefficient
      {
         prdrhsrehihi[i][0] = cffrehihi[col][i][0];
         prdrhsrelohi[i][0] = cffrelohi[col][i][0];
         prdrhsrehilo[i][0] = cffrehilo[col][i][0];
         prdrhsrelolo[i][0] = cffrelolo[col][i][0];
         prdrhsimhihi[i][0] = cffimhihi[col][i][0];
         prdrhsimlohi[i][0] = cffimlohi[col][i][0];
         prdrhsimhilo[i][0] = cffimhilo[col][i][0];
         prdrhsimlolo[i][0] = cffimlolo[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhsrehihi[i][k] = 0.0; prdrhsimhihi[i][k] = 0.0;
            prdrhsrelohi[i][k] = 0.0; prdrhsimlohi[i][k] = 0.0;
            prdrhsrehilo[i][k] = 0.0; prdrhsimhilo[i][k] = 0.0;
            prdrhsrelolo[i][k] = 0.0; prdrhsimlolo[i][k] = 0.0;
         }
      }
      if(vrblvl > 1)
         cout << "Evaluating at column " << col << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         for(int k=0; k<dim; k++) rowsA[i][k] = 0;
         for(int k=0; k<nvr[col][i]; k++) rowsA[i][idx[col][i][k]] = 1;
         if(vrblvl > 1)
         {
            for(int k=0; k<dim; k++) cout << " " << rowsA[i][k];
            cout << endl;
         }
      }
      evaluate_complex4_monomials
         (dim,deg,rowsA,
          xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
          prdrhsrehihi,prdrhsrelohi,prdrhsrehilo,prdrhsrelolo,
          prdrhsimhihi,prdrhsimlohi,prdrhsimhilo,prdrhsimlolo);

      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               // rhsre[i][k] += prdrhsre[i][k];
               qdf_inc(&rhsrehihi[i][k],&rhsrelohi[i][k],
                       &rhsrehilo[i][k],&rhsrelolo[i][k],
                       prdrhsrehihi[i][k],prdrhsrelohi[i][k],
                       prdrhsrehilo[i][k],prdrhsrelolo[i][k]);
               // rhsim[i][k] += prdrhsim[i][k];
               qdf_inc(&rhsimhihi[i][k],&rhsimlohi[i][k],
                       &rhsimhilo[i][k],&rhsimlolo[i][k],
                       prdrhsimhihi[i][k],prdrhsimlohi[i][k],
                       prdrhsimhilo[i][k],prdrhsimlolo[i][k]);
            }
      }
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhsrehihi[i]); free(prdrhsimhihi[i]);
      free(prdrhsrelohi[i]); free(prdrhsimlohi[i]);
      free(prdrhsrehilo[i]); free(prdrhsimhilo[i]);
      free(prdrhsrelolo[i]); free(prdrhsimlolo[i]);
   }
   free(prdrhsrehihi); free(prdrhsimhihi);
   free(prdrhsrelohi); free(prdrhsimlohi);
   free(prdrhsrehilo); free(prdrhsimhilo);
   free(prdrhsrelolo); free(prdrhsimlolo);
}
