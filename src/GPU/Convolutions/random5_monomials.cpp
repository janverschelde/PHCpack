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
   double rndretb,rndreix,rndremi,rndrerg,rndrepk;
   double rndimtb,rndimix,rndimmi,rndimrg,rndimpk;
   double* pluxretb = new double[deg+1];
   double* pluxreix = new double[deg+1];
   double* pluxremi = new double[deg+1];
   double* pluxrerg = new double[deg+1];
   double* pluxrepk = new double[deg+1];
   double* pluximtb = new double[deg+1];
   double* pluximix = new double[deg+1];
   double* pluximmi = new double[deg+1];
   double* pluximrg = new double[deg+1];
   double* pluximpk = new double[deg+1];
   double* minxretb = new double[deg+1];
   double* minxreix = new double[deg+1];
   double* minxremi = new double[deg+1];
   double* minxrerg = new double[deg+1];
   double* minxrepk = new double[deg+1];
   double* minximtb = new double[deg+1];
   double* minximix = new double[deg+1];
   double* minximmi = new double[deg+1];
   double* minximrg = new double[deg+1];
   double* minximpk = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_cmplx5_exponentials(deg,
         &rndretb,&rndreix,&rndremi,&rndrerg,&rndrepk,
         &rndimtb,&rndimix,&rndimmi,&rndimrg,&rndimpk,
         pluxretb,pluxreix,pluxremi,pluxrerg,pluxrepk,
         pluximtb,pluximix,pluximmi,pluximrg,pluximpk,
         minxretb,minxreix,minxremi,minxrerg,minxrepk,
         minximtb,minximix,minximmi,minximrg,minximpk);

      for(int j=0; j<=deg; j++)
      {
         dataretb[i][j] = pluxretb[j]; dataimtb[i][j] = pluximtb[j];
         datareix[i][j] = pluxreix[j]; dataimix[i][j] = pluximix[j];
         dataremi[i][j] = pluxremi[j]; dataimmi[i][j] = pluximmi[j];
         datarerg[i][j] = pluxrerg[j]; dataimrg[i][j] = pluximrg[j];
         datarepk[i][j] = pluxrepk[j]; dataimpk[i][j] = pluximpk[j];
         dataretb[i+1][j] = minxretb[j]; dataimtb[i+1][j] = minximtb[j];
         datareix[i+1][j] = minxreix[j]; dataimix[i+1][j] = minximix[j];
         dataremi[i+1][j] = minxremi[j]; dataimmi[i+1][j] = minximmi[j];
         datarerg[i+1][j] = minxrerg[j]; dataimrg[i+1][j] = minximrg[j];
         datarepk[i+1][j] = minxrepk[j]; dataimpk[i+1][j] = minximpk[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      dataretb[dim-1][0] = 1.0; dataimtb[dim-1][0] = 0.0;
      datareix[dim-1][0] = 0.0; dataimix[dim-1][0] = 0.0;
      dataremi[dim-1][0] = 0.0; dataimmi[dim-1][0] = 0.0;
      datarerg[dim-1][0] = 0.0; dataimrg[dim-1][0] = 0.0;
      datarepk[dim-1][0] = 0.0; dataimpk[dim-1][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         dataretb[dim-1][j] = 0.0; dataimtb[dim-1][j] = 0.0;
         datareix[dim-1][j] = 0.0; dataimix[dim-1][j] = 0.0;
         dataremi[dim-1][j] = 0.0; dataimmi[dim-1][j] = 0.0;
         datarerg[dim-1][j] = 0.0; dataimrg[dim-1][j] = 0.0;
         datarepk[dim-1][j] = 0.0; dataimpk[dim-1][j] = 0.0;
      }
   }
   free(pluxretb); free(pluxreix); free(pluxremi);
   free(pluxrerg); free(pluxrepk);
   free(pluximtb); free(pluximix); free(pluximmi);
   free(pluximrg); free(pluximpk);
   free(minxretb); free(minxreix); free(minxremi);
   free(minxrerg); free(minxrepk);
   free(minximtb); free(minximix); free(minximmi);
   free(minximrg); free(minximpk); 
}
