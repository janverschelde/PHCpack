// The file random5_monomials.cpp defines functions specified 
// in random5_monomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_monomials.h"
#include "penta_double_functions.h"
#include "random5_vectors.h"
#include "random5_series.h"
#include "random5_monomials.h"

bool make_real5_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk )
{
   bool fail;

   if(nvr > dim)
   {
      std::cout << "ERROR: nvr = " << nvr << " > " << dim << " dim"
                << std::endl;

      return true;
   }
   else
   {
      for(int i=0; i<=deg; i++)
         random_penta_double(&cfftb[i],&cffix[i],&cffmi[i],
                             &cffrg[i],&cffpk[i]);

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

bool make_complex5_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk )
{
   bool fail;

   if(nvr > dim)
   {
      std::cout << "ERROR: nvr = " << nvr << " > " << dim << " dim"
                << std::endl;

      return true;
   }
   else
   {
      double rndtb,rndix,rndmi,rndrg,rndpk;
      double sintb,sinix,sinmi,sinrg,sinpk;

      for(int i=0; i<=deg; i++)
      {
         random_penta_double
            (&rndtb,&rndix,&rndmi,&rndrg,&rndpk);       // random cos
                                                        // cos(angle)
         cffretb[i] = rndtb; cffreix[i] = rndix; cffremi[i] = rndmi;
         cffrerg[i] = rndrg; cffrepk[i] = rndpk;
         pdf_sqrt(rndtb,rndix,rndmi,rndrg,rndpk,
                  &sintb,&sinix,&sinmi,&sinrg,&sinpk);      // cos^(angle)
         pdf_minus(&sintb,&sinix,&sinmi,&sinrg,&sinpk);     // -cos^(angle)
         pdf_inc_d(&sintb,&sinix,&sinmi,&sinrg,&sinpk,1.0); // 1-cos^2(angle)
         // sin is sqrt
         pdf_sqrt(sintb,sinix,sinmi,sinrg,sinpk,
                  &cffimtb[i],&cffimix[i],&cffimmi[i],
                  &cffimrg[i],&cffimpk[i]);
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

void make_real5_input
 ( int dim, int deg, double **datatb, double **dataix,
   double **datami, double **datarg, double **datapk )
{
   double rndtb,rndix,rndmi,rndrg,rndpk;
   double* pluxtb = new double[deg+1];
   double* pluxix = new double[deg+1];
   double* pluxmi = new double[deg+1];
   double* pluxrg = new double[deg+1];
   double* pluxpk = new double[deg+1];
   double* minxtb = new double[deg+1];
   double* minxix = new double[deg+1];
   double* minxmi = new double[deg+1];
   double* minxrg = new double[deg+1];
   double* minxpk = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_dbl5_exponentials
         (deg,&rndtb,&rndix,&rndmi,&rndrg,&rndpk,
          pluxtb,pluxix,pluxmi,pluxrg,pluxpk,
          minxtb,minxix,minxmi,minxrg,minxpk);

      for(int j=0; j<=deg; j++)
      {
         datatb[i][j] = pluxtb[j]; dataix[i][j] = pluxix[j];
         datami[i][j] = pluxmi[j]; datarg[i][j] = pluxrg[j];
         datapk[i][j] = pluxpk[j];
         datatb[i+1][j] = minxtb[j]; dataix[i+1][j] = minxix[j];
         datami[i+1][j] = minxmi[j]; datarg[i+1][j] = minxrg[j];
         datapk[i+1][j] = minxpk[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datatb[dim-1][0] = 1.0; dataix[dim-1][0] = 0.0; datami[dim-1][0] = 0.0;
      datarg[dim-1][0] = 0.0; datapk[dim-1][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         datatb[dim-1][j] = 0.0; dataix[dim-1][j] = 0.0;
         datami[dim-1][j] = 0.0; datarg[dim-1][j] = 0.0;
         datapk[dim-1][j] = 0.0;
      }
   }
   free(pluxtb); free(pluxix); free(pluxmi); free(pluxrg); free(pluxpk);
   free(minxtb); free(minxix); free(minxmi); free(minxrg); free(minxpk); 
}

void make_complex5_input
 ( int dim, int deg,
   double **dataretb, double **datareix, double **dataremi,
   double **datarerg, double **datarepk,
   double **dataimtb, double **dataimix, double **dataimmi,
   double **dataimrg, double **dataimpk )
{
}
