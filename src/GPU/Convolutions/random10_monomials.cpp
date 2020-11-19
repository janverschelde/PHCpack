// The file random10_monomials.cpp defines functions specified 
// in random10_monomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_monomials.h"
#include "deca_double_functions.h"
#include "random10_vectors.h"
#include "random10_series.h"
#include "random10_monomials.h"

bool make_real10_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrtb, double *cffrix, double *cffrmi,
   double *cffrrg, double *cffrpk,
   double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk )
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
         random_deca_double
            (&cffrtb[i],&cffrix[i],&cffrmi[i],&cffrrg[i],&cffrpk[i],
             &cffltb[i],&cfflix[i],&cfflmi[i],&cfflrg[i],&cfflpk[i]);

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

bool make_complex10_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrertb, double *cffrerix, double *cffrermi,
   double *cffrerrg, double *cffrerpk,
   double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi,
   double *cffimrrg, double *cffimrpk,
   double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk )
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
      double rndrtb,rndrix,rndrmi,rndrrg,rndrpk;
      double rndltb,rndlix,rndlmi,rndlrg,rndlpk;
      double sinrtb,sinrix,sinrmi,sinrrg,sinrpk;
      double sinltb,sinlix,sinlmi,sinlrg,sinlpk;

      for(int i=0; i<=deg; i++)
      {
         random_deca_double
            (&rndrtb,&rndrix,&rndrmi,&rndrrg,&rndrpk,
             &rndltb,&rndlix,&rndlmi,&rndlrg,&rndlpk);  // random cos

         cffrertb[i] = rndrtb; cffrerix[i] = rndrix; cffrermi[i] = rndrmi;
         cffrerrg[i] = rndrrg; cffrerpk[i] = rndrpk;
         cffreltb[i] = rndltb; cffrelix[i] = rndlix; cffrelmi[i] = rndlmi;
         cffrelrg[i] = rndlrg; cffrelpk[i] = rndlpk;

         daf_sqrt(rndrtb,rndrix,rndrmi,rndrrg,rndrpk,
                  rndltb,rndlix,rndlmi,rndlrg,rndlpk,
                  &sinrtb,&sinrix,&sinrmi,&sinrrg,&sinrpk,
                  &sinltb,&sinlix,&sinlmi,&sinlrg,&sinlpk); // cos^(angle)
         daf_minus(&sinrtb,&sinrix,&sinrmi,&sinrrg,&sinrpk,
                   &sinltb,&sinlix,&sinlmi,&sinlrg,&sinlpk); // -cos^(angle)
         daf_inc_d(&sinrtb,&sinrix,&sinrmi,&sinrrg,&sinrpk,
                   &sinltb,&sinlix,&sinlmi,&sinlrg,&sinlpk,1.0);
         // sin is sqrt of 1-cos^2(angle)
         daf_sqrt(sinrtb,sinrix,sinrmi,sinrrg,sinrpk,
                  sinltb,sinlix,sinlmi,sinlrg,sinlpk,
                  &cffimrtb[i],&cffimrix[i],&cffimrmi[i],
                  &cffimrrg[i],&cffimrpk[i],
                  &cffimltb[i],&cffimlix[i],&cffimlmi[i],
                  &cffimlrg[i],&cffimlpk[i]);
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

void make_real10_input
 ( int dim, int deg,
   double **datartb, double **datarix, double **datarmi,
   double **datarrg, double **datarpk,
   double **dataltb, double **datalix, double **datalmi,
   double **datalrg, double **datalpk )
{
   double rndrtb,rndrix,rndrmi,rndrrg,rndrpk;
   double rndltb,rndlix,rndlmi,rndlrg,rndlpk;
   double* pluxrtb = new double[deg+1];
   double* pluxrix = new double[deg+1];
   double* pluxrmi = new double[deg+1];
   double* pluxrrg = new double[deg+1];
   double* pluxrpk = new double[deg+1];
   double* pluxltb = new double[deg+1];
   double* pluxlix = new double[deg+1];
   double* pluxlmi = new double[deg+1];
   double* pluxlrg = new double[deg+1];
   double* pluxlpk = new double[deg+1];
   double* minxrtb = new double[deg+1];
   double* minxrix = new double[deg+1];
   double* minxrmi = new double[deg+1];
   double* minxrrg = new double[deg+1];
   double* minxrpk = new double[deg+1];
   double* minxltb = new double[deg+1];
   double* minxlix = new double[deg+1];
   double* minxlmi = new double[deg+1];
   double* minxlrg = new double[deg+1];
   double* minxlpk = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_dbl10_exponentials
         (deg,&rndrtb,&rndrix,&rndrmi,&rndrrg,&rndrpk,
              &rndltb,&rndlix,&rndlmi,&rndlrg,&rndlpk,
          pluxrtb,pluxrix,pluxrmi,pluxrrg,pluxrpk,
          pluxltb,pluxlix,pluxlmi,pluxlrg,pluxlpk,
          minxrtb,minxrix,minxrmi,minxrrg,minxrpk,
          minxltb,minxlix,minxlmi,minxlrg,minxlpk);

      for(int j=0; j<=deg; j++)
      {
         datartb[i][j] = pluxrtb[j]; datarix[i][j] = pluxrix[j];
         datarmi[i][j] = pluxrmi[j]; datarrg[i][j] = pluxrrg[j];
         datarpk[i][j] = pluxrpk[j];
         dataltb[i][j] = pluxltb[j]; datalix[i][j] = pluxlix[j];
         datalmi[i][j] = pluxlmi[j]; datalrg[i][j] = pluxlrg[j];
         datalpk[i][j] = pluxlpk[j];
         datartb[i+1][j] = minxrtb[j]; datarix[i+1][j] = minxrix[j];
         datarmi[i+1][j] = minxrmi[j]; datarrg[i+1][j] = minxrrg[j];
         datarpk[i+1][j] = minxrpk[j];
         dataltb[i+1][j] = minxltb[j]; datalix[i+1][j] = minxlix[j];
         datalmi[i+1][j] = minxlmi[j]; datalrg[i+1][j] = minxlrg[j];
         datalpk[i+1][j] = minxlpk[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datartb[dim-1][0] = 1.0;
      datarix[dim-1][0] = 0.0; datarmi[dim-1][0] = 0.0;
      datarrg[dim-1][0] = 0.0; datarpk[dim-1][0] = 0.0;
      dataltb[dim-1][0] = 0.0;
      datalix[dim-1][0] = 0.0; datalmi[dim-1][0] = 0.0;
      datalrg[dim-1][0] = 0.0; datalpk[dim-1][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         datartb[dim-1][j] = 0.0; datarix[dim-1][j] = 0.0;
         datarmi[dim-1][j] = 0.0; datarrg[dim-1][j] = 0.0;
         datarpk[dim-1][j] = 0.0;
         dataltb[dim-1][j] = 0.0; datalix[dim-1][j] = 0.0;
         datalmi[dim-1][j] = 0.0; datalrg[dim-1][j] = 0.0;
         datalpk[dim-1][j] = 0.0;
      }
      int middle = dim/2; // move the last input series to the middle
      double tmp;

      for(int j=0; j<=deg; j++)
      {
         tmp = datartb[dim-1][j]; datartb[dim-1][j] = datartb[middle][j];
         datartb[middle][j] = tmp;
         tmp = datarix[dim-1][j]; datarix[dim-1][j] = datarix[middle][j];
         datarix[middle][j] = tmp;
         tmp = datarmi[dim-1][j]; datarmi[dim-1][j] = datarmi[middle][j];
         datarmi[middle][j] = tmp;
         tmp = datarrg[dim-1][j]; datarrg[dim-1][j] = datarrg[middle][j];
         datarrg[middle][j] = tmp;
         tmp = datarpk[dim-1][j]; datarpk[dim-1][j] = datarpk[middle][j];
         datarpk[middle][j] = tmp;
         tmp = dataltb[dim-1][j]; dataltb[dim-1][j] = dataltb[middle][j];
         dataltb[middle][j] = tmp;
         tmp = datalix[dim-1][j]; datalix[dim-1][j] = datalix[middle][j];
         datalix[middle][j] = tmp;
         tmp = datalmi[dim-1][j]; datalmi[dim-1][j] = datalmi[middle][j];
         datalmi[middle][j] = tmp;
         tmp = datalrg[dim-1][j]; datalrg[dim-1][j] = datalrg[middle][j];
         datalrg[middle][j] = tmp;
         tmp = datalpk[dim-1][j]; datalpk[dim-1][j] = datalpk[middle][j];
         datalpk[middle][j] = tmp;
      }
   }
   free(pluxrtb); free(pluxrix); free(pluxrmi); free(pluxrrg); free(pluxrpk);
   free(pluxltb); free(pluxlix); free(pluxlmi); free(pluxlrg); free(pluxlpk);
   free(minxrtb); free(minxrix); free(minxrmi); free(minxrrg); free(minxrpk); 
   free(minxltb); free(minxlix); free(minxlmi); free(minxlrg); free(minxlpk); 
}

void make_complex10_input
 ( int dim, int deg,
   double **datarertb, double **datarerix, double **datarermi,
   double **datarerrg, double **datarerpk,
   double **datareltb, double **datarelix, double **datarelmi,
   double **datarelrg, double **datarelpk,
   double **dataimrtb, double **dataimrix, double **dataimrmi,
   double **dataimrrg, double **dataimrpk,
   double **dataimltb, double **dataimlix, double **dataimlmi,
   double **dataimlrg, double **dataimlpk )
{
   double rndrertb,rndrerix,rndrermi,rndrerrg,rndrerpk;
   double rndreltb,rndrelix,rndrelmi,rndrelrg,rndrelpk;
   double rndimrtb,rndimrix,rndimrmi,rndimrrg,rndimrpk;
   double rndimltb,rndimlix,rndimlmi,rndimlrg,rndimlpk;
   double* pluxrertb = new double[deg+1];
   double* pluxrerix = new double[deg+1];
   double* pluxrermi = new double[deg+1];
   double* pluxrerrg = new double[deg+1];
   double* pluxrerpk = new double[deg+1];
   double* pluxreltb = new double[deg+1];
   double* pluxrelix = new double[deg+1];
   double* pluxrelmi = new double[deg+1];
   double* pluxrelrg = new double[deg+1];
   double* pluxrelpk = new double[deg+1];
   double* pluximrtb = new double[deg+1];
   double* pluximrix = new double[deg+1];
   double* pluximrmi = new double[deg+1];
   double* pluximrrg = new double[deg+1];
   double* pluximrpk = new double[deg+1];
   double* pluximltb = new double[deg+1];
   double* pluximlix = new double[deg+1];
   double* pluximlmi = new double[deg+1];
   double* pluximlrg = new double[deg+1];
   double* pluximlpk = new double[deg+1];
   double* minxrertb = new double[deg+1];
   double* minxrerix = new double[deg+1];
   double* minxrermi = new double[deg+1];
   double* minxrerrg = new double[deg+1];
   double* minxrerpk = new double[deg+1];
   double* minxreltb = new double[deg+1];
   double* minxrelix = new double[deg+1];
   double* minxrelmi = new double[deg+1];
   double* minxrelrg = new double[deg+1];
   double* minxrelpk = new double[deg+1];
   double* minximrtb = new double[deg+1];
   double* minximrix = new double[deg+1];
   double* minximrmi = new double[deg+1];
   double* minximrrg = new double[deg+1];
   double* minximrpk = new double[deg+1];
   double* minximltb = new double[deg+1];
   double* minximlix = new double[deg+1];
   double* minximlmi = new double[deg+1];
   double* minximlrg = new double[deg+1];
   double* minximlpk = new double[deg+1];

   for(int i=0; i<dim-1; i=i+2)
   {
      random_cmplx10_exponentials(deg,
         &rndrertb,&rndrerix,&rndrermi,&rndrerrg,&rndrerpk,
         &rndreltb,&rndrelix,&rndrelmi,&rndrelrg,&rndrelpk,
         &rndimrtb,&rndimrix,&rndimrmi,&rndimrrg,&rndimrpk,
         &rndimltb,&rndimlix,&rndimlmi,&rndimlrg,&rndimlpk,
         pluxrertb,pluxrerix,pluxrermi,pluxrerrg,pluxrerpk,
         pluxreltb,pluxrelix,pluxrelmi,pluxrelrg,pluxrelpk,
         pluximrtb,pluximrix,pluximrmi,pluximrrg,pluximrpk,
         pluximltb,pluximlix,pluximlmi,pluximlrg,pluximlpk,
         minxrertb,minxrerix,minxrermi,minxrerrg,minxrerpk,
         minxreltb,minxrelix,minxrelmi,minxrelrg,minxrelpk,
         minximrtb,minximrix,minximrmi,minximrrg,minximrpk,
         minximltb,minximlix,minximlmi,minximlrg,minximlpk);

      for(int j=0; j<=deg; j++)
      {
         datarertb[i][j] = pluxrertb[j]; dataimrtb[i][j] = pluximrtb[j];
         datarerix[i][j] = pluxrerix[j]; dataimrix[i][j] = pluximrix[j];
         datarermi[i][j] = pluxrermi[j]; dataimrmi[i][j] = pluximrmi[j];
         datarerrg[i][j] = pluxrerrg[j]; dataimrrg[i][j] = pluximrrg[j];
         datarerpk[i][j] = pluxrerpk[j]; dataimrpk[i][j] = pluximrpk[j];
         datareltb[i][j] = pluxreltb[j]; dataimltb[i][j] = pluximltb[j];
         datarelix[i][j] = pluxrelix[j]; dataimlix[i][j] = pluximlix[j];
         datarelmi[i][j] = pluxrelmi[j]; dataimlmi[i][j] = pluximlmi[j];
         datarelrg[i][j] = pluxrelrg[j]; dataimlrg[i][j] = pluximlrg[j];
         datarelpk[i][j] = pluxrelpk[j]; dataimlpk[i][j] = pluximlpk[j];
         datarertb[i+1][j] = minxrertb[j]; dataimrtb[i+1][j] = minximrtb[j];
         datarerix[i+1][j] = minxrerix[j]; dataimrix[i+1][j] = minximrix[j];
         datarermi[i+1][j] = minxrermi[j]; dataimrmi[i+1][j] = minximrmi[j];
         datarerrg[i+1][j] = minxrerrg[j]; dataimrrg[i+1][j] = minximrrg[j];
         datarerpk[i+1][j] = minxrerpk[j]; dataimrpk[i+1][j] = minximrpk[j];
         datareltb[i+1][j] = minxreltb[j]; dataimltb[i+1][j] = minximltb[j];
         datarelix[i+1][j] = minxrelix[j]; dataimlix[i+1][j] = minximlix[j];
         datarelmi[i+1][j] = minxrelmi[j]; dataimlmi[i+1][j] = minximlmi[j];
         datarelrg[i+1][j] = minxrelrg[j]; dataimlrg[i+1][j] = minximlrg[j];
         datarelpk[i+1][j] = minxrelpk[j]; dataimlpk[i+1][j] = minximlpk[j];
      }
   }
   if(dim % 2 == 1) // in odd case, set the last input series to one
   {
      datarertb[dim-1][0] = 1.0; dataimrtb[dim-1][0] = 0.0;
      datarerix[dim-1][0] = 0.0; dataimrix[dim-1][0] = 0.0;
      datarermi[dim-1][0] = 0.0; dataimrmi[dim-1][0] = 0.0;
      datarerrg[dim-1][0] = 0.0; dataimrrg[dim-1][0] = 0.0;
      datarerpk[dim-1][0] = 0.0; dataimrpk[dim-1][0] = 0.0;
      datareltb[dim-1][0] = 0.0; dataimltb[dim-1][0] = 0.0;
      datarelix[dim-1][0] = 0.0; dataimlix[dim-1][0] = 0.0;
      datarelmi[dim-1][0] = 0.0; dataimlmi[dim-1][0] = 0.0;
      datarelrg[dim-1][0] = 0.0; dataimlrg[dim-1][0] = 0.0;
      datarelpk[dim-1][0] = 0.0; dataimlpk[dim-1][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         datarertb[dim-1][j] = 0.0; dataimrtb[dim-1][j] = 0.0;
         datarerix[dim-1][j] = 0.0; dataimrix[dim-1][j] = 0.0;
         datarermi[dim-1][j] = 0.0; dataimrmi[dim-1][j] = 0.0;
         datarerrg[dim-1][j] = 0.0; dataimrrg[dim-1][j] = 0.0;
         datarerpk[dim-1][j] = 0.0; dataimrpk[dim-1][j] = 0.0;
         datareltb[dim-1][j] = 0.0; dataimltb[dim-1][j] = 0.0;
         datarelix[dim-1][j] = 0.0; dataimlix[dim-1][j] = 0.0;
         datarelmi[dim-1][j] = 0.0; dataimlmi[dim-1][j] = 0.0;
         datarelrg[dim-1][j] = 0.0; dataimlrg[dim-1][j] = 0.0;
         datarelpk[dim-1][j] = 0.0; dataimlpk[dim-1][j] = 0.0;
      }
   }
   free(pluxrertb); free(pluxrerix); free(pluxrermi);
   free(pluxrerrg); free(pluxrerpk);
   free(pluxreltb); free(pluxrelix); free(pluxrelmi);
   free(pluxrelrg); free(pluxrelpk);
   free(pluximrtb); free(pluximrix); free(pluximrmi);
   free(pluximrrg); free(pluximrpk);
   free(pluximltb); free(pluximlix); free(pluximlmi);
   free(pluximlrg); free(pluximlpk);
   free(minxrertb); free(minxrerix); free(minxrermi);
   free(minxrerrg); free(minxrerpk);
   free(minxreltb); free(minxrelix); free(minxrelmi);
   free(minxrelrg); free(minxrelpk);
   free(minximrtb); free(minximrix); free(minximrmi);
   free(minximrrg); free(minximrpk); 
   free(minximltb); free(minximlix); free(minximlmi);
   free(minximlrg); free(minximlpk); 
}
