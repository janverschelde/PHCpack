// The file dbl4_systems_host.cpp defines the functions with prototypes in
// the file dbl4_systems_host.h.

#include <iostream>
#include <iomanip>
#include "quad_double_functions.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"

using namespace std;

void CPU_dbl4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_dbl4_product
                     (deg,inputhihi[idxvar],inputlohi[idxvar],
                          inputhilo[idxvar],inputlolo[idxvar],
                      cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
                      acchihi,   acclohi,   acchilo,   acclolo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffhihi[i][L] = acchihi[L];
                     cfflohi[i][L] = acclohi[L];
                     cffhilo[i][L] = acchilo[L];
                     cfflolo[i][L] = acclolo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << cffhihi[i][j] << "  " << cfflohi[i][j] << endl;
            cout << cffhilo[i][j] << "  " << cfflolo[i][j] << endl;
         }
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl;
            cout << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
         }
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl4_evaldiff
         (dim,nvr[i],deg,idx[i],cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
           inputhihi,    inputlohi,    inputhilo,    inputlolo,
          outputhihi[i],outputlohi[i],outputhilo[i],outputlolo[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << outputhihi[i][dim][j] << "  "
                 << outputlohi[i][dim][j] << endl;
            cout << outputhilo[i][dim][j] << "  "
                 << outputlolo[i][dim][j] << endl;
         }
      }
   }

   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];
               double tmphihi,tmplohi,tmphilo,tmplolo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // output[i][idxvar][k] = factor*output[i][idxvar][k];
                  qdf_mul_qd_d
                     (outputhihi[i][idxvar][k],outputlohi[i][idxvar][k],
                      outputhilo[i][idxvar][k],outputlolo[i][idxvar][k],
                      factor,&tmphihi,&tmplohi,&tmphilo,&tmplolo);
                  outputhihi[i][idxvar][k] = tmphihi;
                  outputlohi[i][idxvar][k] = tmplohi;
                  outputhilo[i][idxvar][k] = tmphilo;
                  outputlolo[i][idxvar][k] = tmplolo;
               }
            }
         }
      }
   }
}

void CPU_cmplx4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi, double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi, double *accimhilo, double *accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi, 
   double **inputimhilo, double **inputimlolo, 
   double ***outputrehihi, double ***outputrelohi, 
   double ***outputrehilo, double ***outputrelolo, 
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo, int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_cmplx4_product
                     (deg,inputrehihi[idxvar],inputrelohi[idxvar],
                          inputrehilo[idxvar],inputrelolo[idxvar],
                          inputimhihi[idxvar],inputimlohi[idxvar],
                          inputimhilo[idxvar],inputimlolo[idxvar],
                      cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
                      cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
                      accrehihi,accrelohi,accrehilo,accrelolo,
                      accimhihi,accimlohi,accimhilo,accimlolo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffrehihi[i][L] = accrehihi[L];
                     cffrelohi[i][L] = accrelohi[L];
                     cffrehilo[i][L] = accrehilo[L];
                     cffrelolo[i][L] = accrelolo[L];
                     cffimhihi[i][L] = accimhihi[L];
                     cffimlohi[i][L] = accimlohi[L];
                     cffimhilo[i][L] = accimhilo[L];
                     cffimlolo[i][L] = accimlolo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffrehihi[i][j] << "  " << cffrelohi[i][j] << endl << "  "
                 << cffrehilo[i][j] << "  " << cffrelolo[i][j] << endl << "  "
                 << cffimhihi[i][j] << "  " << cffimlohi[i][j] << endl << "  "
                 << cffimhilo[i][j] << "  " << cffimlolo[i][j] << endl;
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehihi[i][j] << "  " << inputrelohi[i][j] << endl
                 << "  "
                 << inputrehilo[i][j] << "  " << inputrelolo[i][j] << endl
                 << "  "
                 << inputimhihi[i][j] << "  " << inputimlohi[i][j] << endl
                 << "  "
                 << inputimhilo[i][j] << "  " << inputimlolo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_cmplx4_evaldiff
         (dim,nvr[i],deg,idx[i],
          cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
          cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
          inputrehihi,inputrelohi,inputrehilo,inputrelolo,
          inputimhihi,inputimlohi,inputimhilo,inputimlolo,
          outputrehihi[i],outputrelohi[i],outputrehilo[i],outputrelolo[i],
          outputimhihi[i],outputimlohi[i],outputimhilo[i],outputimlolo[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehihi[i][dim][j] << "  "
                 << outputrelohi[i][dim][j] << endl << "  "
                 << outputrehilo[i][dim][j] << "  "
                 << outputrelolo[i][dim][j] << endl << "  "
                 << outputimhihi[i][dim][j] << "  "
                 << outputimlohi[i][dim][j] << endl << "  "
                 << outputimhilo[i][dim][j] << "  "
                 << outputimlolo[i][dim][j] << endl;
      }
   }
   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];
               double acchihi,acclohi,acchilo,acclolo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // outputre[i][idxvar][k] = factor*outputre[i][idxvar][k];
                  qdf_mul_qd_d
                     (outputrehihi[i][idxvar][k],outputrelohi[i][idxvar][k],
                      outputrehilo[i][idxvar][k],outputrelolo[i][idxvar][k],
                      factor,&acchihi,&acclohi,&acchilo,&acclolo);
                  outputrehihi[i][idxvar][k] = acchihi;
                  outputrelohi[i][idxvar][k] = acclohi;
                  outputrehilo[i][idxvar][k] = acchilo;
                  outputrelolo[i][idxvar][k] = acclolo;
                  // outputim[i][idxvar][k] = factor*outputim[i][idxvar][k];
                  qdf_mul_qd_d
                     (outputimhihi[i][idxvar][k],outputimlohi[i][idxvar][k],
                      outputimhilo[i][idxvar][k],outputimlolo[i][idxvar][k],
                      factor,&acchihi,&acclohi,&acchilo,&acclolo);
                  outputimhihi[i][idxvar][k] = acchihi;
                  outputimlohi[i][idxvar][k] = acclohi;
                  outputimhilo[i][idxvar][k] = acchilo;
                  outputimlolo[i][idxvar][k] = acclolo;
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "after multiplication with the factors ..." << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehihi[i][dim][j] << "  "
                 << outputrelohi[i][dim][j] << endl << "  "
                 << outputrehilo[i][dim][j] << "  "
                 << outputrelolo[i][dim][j] << endl << "  "
                 << outputimhihi[i][dim][j] << "  "
                 << outputimlohi[i][dim][j] << endl << "  "
                 << outputimhilo[i][dim][j] << "  "
                 << outputimlolo[i][dim][j] << endl;
      }
   }
}

void CPU_dbl4_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double **acchihi, double **acclohi, double **acchilo, double **acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhihi[i][j] = 0.0; funvallohi[i][j] = 0.0;
         funvalhilo[i][j] = 0.0; funvallolo[i][j] = 0.0;
      }
   funvalhihi[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalhihi[k][i][j] = 0.0; jacvallohi[k][i][j] = 0.0;
            jacvalhilo[k][i][j] = 0.0; jacvallolo[k][i][j] = 0.0;
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffhihi[i][j][k] << "  " << cfflohi[i][j][k] << endl
                    << "  "
                    << cffhilo[i][j][k] << "  " << cfflolo[i][j][k] << endl;
         }
      }
      cout << "dim : " << dim << endl;
      for(int i=0; i<nbrcol; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";
            cout << "  idx[" << i << "][" << j << "] :";
            for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
            cout << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl
                 << "  "
                 << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
      }
   }
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // evaluate j-th monomial in column i
         {
            for(int k=0; k<=dim; k++)
               for(int L=0; L<degp1; L++)
               {
                  acchihi[k][L] = 0.0; acclohi[k][L] = 0.0;
                  acchilo[k][L] = 0.0; acclolo[k][L] = 0.0;
               }

            CPU_dbl4_evaldiff
               (dim,nvr[i][j],deg,idx[i][j],
                cffhihi[i][j],cfflohi[i][j],cffhilo[i][j],cfflolo[i][j],
                inputhihi,inputlohi,inputhilo,inputlolo,
                acchihi,acclohi,acchilo,acclolo);

            for(int L=0; L<degp1; L++) // funval[j][L] += acc[dim][L];
            {
               qdf_inc(&funvalhihi[j][L],&funvallohi[j][L],
                       &funvalhilo[j][L],&funvallolo[j][L],
                       acchihi[dim][L],acclohi[dim][L],
                       acchilo[dim][L],acclolo[dim][L]);
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  // jacval[L][j][idxval] += acc[idxval][L];
                  qdf_inc(&jacvalhihi[L][j][idxval],
                          &jacvallohi[L][j][idxval],
                          &jacvalhilo[L][j][idxval],
                          &jacvallolo[L][j][idxval],
                          acchihi[idxval][L],acclohi[idxval][L],
                          acchilo[idxval][L],acclolo[idxval][L]);
               }
            }
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<degp1; j++)
            cout << funvalhihi[i][j] << "  " << funvallohi[i][j] << endl
                 << "  "
                 << funvalhilo[i][j] << "  " << funvallolo[i][j] << endl;
      }
   }
}

void CPU_cmplx4_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **accrehihi, double **accrelohi,
   double **accrehilo, double **accrelolo,
   double **accimhihi, double **accimlohi,
   double **accimhilo, double **accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehihi[i][j] = 0.0; funvalrelohi[i][j] = 0.0;
         funvalrehilo[i][j] = 0.0; funvalrelolo[i][j] = 0.0;
         funvalimhihi[i][j] = 0.0; funvalimlohi[i][j] = 0.0;
         funvalimhilo[i][j] = 0.0; funvalimlolo[i][j] = 0.0;
      }
   funvalrehihi[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalrehihi[k][i][j] = 0.0; jacvalrelohi[k][i][j] = 0.0;
            jacvalrehilo[k][i][j] = 0.0; jacvalrelolo[k][i][j] = 0.0;
            jacvalimhihi[k][i][j] = 0.0; jacvalimlohi[k][i][j] = 0.0;
            jacvalimhilo[k][i][j] = 0.0; jacvalimlolo[k][i][j] = 0.0;
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffrehihi[i][j][k] << "  "
                    << cffrelohi[i][j][k] << endl
                    << "  "
                    << cffrehilo[i][j][k] << "  "
                    << cffrelolo[i][j][k] << endl
                    << "  "
                    << cffimhihi[i][j][k] << "  "
                    << cffimlohi[i][j][k] << endl
                    << "  "
                    << cffimhilo[i][j][k] << "  "
                    << cffimlolo[i][j][k] << endl;
         }
      }
      cout << "dim : " << dim << endl;
      for(int i=0; i<nbrcol; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";
            cout << "  idx[" << i << "][" << j << "] :";
            for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
            cout << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehihi[i][j] << "  " << inputrelohi[i][j] << endl
                 << "  "
                 << inputrehilo[i][j] << "  " << inputrelolo[i][j] << endl
                 << "  "
                 << inputimhihi[i][j] << "  " << inputimlohi[i][j] << endl
                 << "  "
                 << inputimhilo[i][j] << "  " << inputimlolo[i][j] << endl;
      }
   }
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // evaluate j-th monomial in column i
         {
            CPU_cmplx4_evaldiff
               (dim,nvr[i][j],deg,idx[i][j],
                cffrehihi[i][j],cffrelohi[i][j],
                cffrehilo[i][j],cffrelolo[i][j],
                cffimhihi[i][j],cffimlohi[i][j],
                cffimhilo[i][j],cffimlolo[i][j],
                inputrehihi,inputrelohi,inputrehilo,inputrelolo,
                inputimhihi,inputimlohi,inputimhilo,inputimlolo,
                accrehihi,accrelohi,accrehilo,accrelolo,
                accimhihi,accimlohi,accimhilo,accimlolo);

            for(int L=0; L<degp1; L++)
            {
               // funvalre[j][L] += accre[dim][L];
               qdf_inc(&funvalrehihi[j][L],&funvalrelohi[j][L],
                       &funvalrehilo[j][L],&funvalrelolo[j][L],
                       accrehihi[dim][L],accrelohi[dim][L],
                       accrehilo[dim][L],accrelolo[dim][L]);
               // funvalim[j][L] += accim[dim][L];
               qdf_inc(&funvalimhihi[j][L],&funvalimlohi[j][L],
                       &funvalimhilo[j][L],&funvalimlolo[j][L],
                       accimhihi[dim][L],accimlohi[dim][L],
                       accimhilo[dim][L],accimlolo[dim][L]);
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  // jacvalre[L][j][idxval] += accre[idxval][L];
                  qdf_inc(&jacvalrehihi[L][j][idxval],
                          &jacvalrelohi[L][j][idxval],
                          &jacvalrehilo[L][j][idxval],
                          &jacvalrelolo[L][j][idxval],
                          accrehihi[idxval][L],accrelohi[idxval][L],
                          accrehilo[idxval][L],accrelolo[idxval][L]);
                  // jacvalim[L][j][idxval] += accim[idxval][L];
                  qdf_inc(&jacvalimhihi[L][j][idxval],
                          &jacvalimlohi[L][j][idxval],
                          &jacvalimhilo[L][j][idxval],
                          &jacvalimlolo[L][j][idxval],
                          accimhihi[idxval][L],accimlohi[idxval][L],
                          accimhilo[idxval][L],accimlolo[idxval][L]);
               }
            }
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << funvalrehihi[i][j] << "  "
                 << funvalrelohi[i][j] << endl << "  "
                 << funvalrehilo[i][j] << "  "
                 << funvalrelolo[i][j] << endl << "  "
                 << funvalimhihi[i][j] << "  "
                 << funvalimlohi[i][j] << endl << "  "
                 << funvalimhilo[i][j] << "  "
                 << funvalimlolo[i][j] << endl;
      }
   }
}

void dbl4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbhihi, double **mblohi, double **mbhilo, double **mblolo,
   double damper,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi, 
   double **funvalhilo, double **funvallolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl )
{
   double acchihi,acclohi,acchilo,acclolo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhihi[i][j] = outputhihi[i][dim][j];
         funvallohi[i][j] = outputlohi[i][dim][j];
         funvalhilo[i][j] = outputhilo[i][dim][j];
         funvallolo[i][j] = outputlolo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << funvalhihi[i][0] << "  "
                            << funvallohi[i][0] << endl;
         cout << "     " << funvalhilo[i][0] << "  "
                         << funvallolo[i][0] << endl;
      }
   }
   /*
    * Linearize the function values in the rhs and swap sign,
    * but keep in mind that the right hand side is mbre + I*mbim,
    * as the monomial system is x^U = rhs, or x^U - rhs = 0,
    * so from the funval we substract rhs and then flip sign
    * in the computation of the rhs of the linear system.
    */
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         // rhs[i][j] = -(funval[j][i] - mb[j][i]);
         qdf_sub(funvalhihi[j][i],funvallohi[j][i],
                 funvalhilo[j][i],funvallolo[j][i],
                 mbhihi[j][i],mblohi[j][i],mbhilo[j][i],mblolo[j][i],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhshihi[i][j] = -acchihi;
         rhslohi[i][j] = -acclohi;
         rhshilo[i][j] = -acchilo;
         rhslolo[i][j] = -acclolo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << rhshihi[i][j] << "  " << rhslohi[i][j] << endl;
            cout << rhshilo[i][j] << "  " << rhslolo[i][j] << endl;
         }
      }
   }
   // output[i][idx[i][k]] for k from 0 to nvr[i] has the coefficients
   // of the series of the derivative of the i-th monomial with respect
   // to the variable idx[i][k].
   for(int i=0; i<dim; i++)          // the i-th monomial
   {
      int *indexes = idx[i];         // indices of the variables
      for(int k=0; k<nvr[i]; k++)    // derivative w.r.t. idx[i][k]
      {                              // has j-th coefficient
         int idxval = indexes[k];

         for(int j=0; j<degp1; j++) 
         {
            jacvalhihi[j][i][idxval] = outputhihi[i][idxval][j];
            jacvallohi[j][i][idxval] = outputlohi[i][idxval][j];
            jacvalhilo[j][i][idxval] = outputhilo[i][idxval][j];
            jacvallolo[j][i][idxval] = outputlolo[i][idxval][j];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
         {
            cout << jacvalhihi[0][i][j]
                 << "  " << jacvallohi[0][i][j] << endl;
            cout << jacvalhilo[0][i][j]
                 << "  " << jacvallolo[0][i][j] << endl;
         }
      }
   }
}

void cmplx4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
   double damper,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl )
{
   double acchihi,acclohi,acchilo,acclolo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehihi[i][j] = outputrehihi[i][dim][j];
         funvalrelohi[i][j] = outputrelohi[i][dim][j];
         funvalrehilo[i][j] = outputrehilo[i][dim][j];
         funvalrelolo[i][j] = outputrelolo[i][dim][j];
         funvalimhihi[i][j] = outputimhihi[i][dim][j];
         funvalimlohi[i][j] = outputimlohi[i][dim][j];
         funvalimhilo[i][j] = outputimhilo[i][dim][j];
         funvalimlolo[i][j] = outputimlolo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalrehihi[i][0] << "  " << funvalrelohi[i][0] << endl
              << "  "
              << funvalrehilo[i][0] << "  " << funvalrelolo[i][0] << endl
              << "  "
              << funvalimhihi[i][0] << "  " << funvalimlohi[i][0] << endl
              << "  "
              << funvalimhilo[i][0] << "  " << funvalimlolo[i][0] << endl;
   }
   /*
    * Linearize the function values in the rhs and swap sign,
    * but keep in mind that the right hand side is mbre + I*mbim,
    * as the monomial system is x^U = rhs, or x^U - rhs = 0,
    * so from the funval we substract rhs and then flip sign
    * in the computation of the rhs of the linear system.
    */

   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         // rhsre[i][j] = -(funvalre[j][i] - mbre[j][i]);
         qdf_sub(funvalrehihi[j][i],funvalrelohi[j][i],
                 funvalrehilo[j][i],funvalrelolo[j][i],
                 mbrehihi[j][i],mbrelohi[j][i],
                 mbrehilo[j][i],mbrelolo[j][i],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhsrehihi[i][j] = -acchihi;
         rhsrelohi[i][j] = -acclohi;
         rhsrehilo[i][j] = -acchilo;
         rhsrelolo[i][j] = -acclolo;
         // rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
         qdf_sub(funvalimhihi[j][i],funvalimlohi[j][i],
                 funvalimhilo[j][i],funvalimlolo[j][i],
                 mbimhihi[j][i],mbimlohi[j][i],
                 mbimhilo[j][i],mbimlolo[j][i],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhsimhihi[i][j] = -acchihi;
         rhsimlohi[i][j] = -acclohi;
         rhsimhilo[i][j] = -acchilo;
         rhsimlolo[i][j] = -acclolo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsrehihi[i][j] << "  " << rhsrelohi[i][j] << endl << "  "
                 << rhsrehilo[i][j] << "  " << rhsrelolo[i][j] << endl << "  "
                 << rhsimhihi[i][j] << "  " << rhsimlohi[i][j] << endl << "  "
                 << rhsimhilo[i][j] << "  " << rhsimlolo[i][j] << endl;
      }
   }
   // output[i][idx[i][k]] for k from 0 to nvr[i] has the coefficients
   // of the series of the derivative of the i-th monomial with respect
   // to the variable idx[i][k].
   for(int i=0; i<dim; i++)          // the i-th monomial
   {
      int *indexes = idx[i];         // indices of the variables
      for(int k=0; k<nvr[i]; k++)    // derivative w.r.t. idx[i][k]
      {                              // has j-th coefficient
         int idxval = indexes[k];

         for(int j=0; j<degp1; j++) 
         {
            jacvalrehihi[j][i][idxval] = outputrehihi[i][idxval][j];
            jacvalrelohi[j][i][idxval] = outputrelohi[i][idxval][j];
            jacvalrehilo[j][i][idxval] = outputrehilo[i][idxval][j];
            jacvalrelolo[j][i][idxval] = outputrelolo[i][idxval][j];
            jacvalimhihi[j][i][idxval] = outputimhihi[i][idxval][j];
            jacvalimlohi[j][i][idxval] = outputimlohi[i][idxval][j];
            jacvalimhilo[j][i][idxval] = outputimhilo[i][idxval][j];
            jacvalimlolo[j][i][idxval] = outputimlolo[i][idxval][j];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalrehihi[0][i][j] << "  "
                 << jacvalrelohi[0][i][j] << endl << "  "
                 << jacvalrehilo[0][i][j] << "  "
                 << jacvalrelolo[0][i][j] << endl << "  "
                 << jacvalimhihi[0][i][j] << "  "
                 << jacvalimlohi[0][i][j] << endl << "  "
                 << jacvalimhilo[0][i][j] << "  "
                 << jacvalimlolo[0][i][j] << endl;
      }
   }
}

void dbl4_define_rhs
 ( int dim, int degp1,
   double **mbhihi, double **mblohi, double **mbhilo, double **mblolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double **rhshihi, double **rhslohi,
   double **rhshilo, double **rhslolo, int vrblvl )
{
   double acchihi,acclohi,acchilo,acclolo;

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalhihi[i][0] << "  " << funvallohi[i][0] << endl
              << "  "
              << funvalhilo[i][0] << "  " << funvallolo[i][0] << endl;
   }
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++) // rhs[i][j] = -(funval[j][i] - mb[j][i]);
      {
         qdf_sub(funvalhihi[j][i],funvallohi[j][i],
                 funvalhilo[j][i],funvallolo[j][i],
                 mbhihi[j][i],mblohi[j][i],mbhilo[j][i],mblolo[j][i],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhshihi[i][j] = -acchihi;
         rhslohi[i][j] = -acclohi;
         rhshilo[i][j] = -acchilo;
         rhslolo[i][j] = -acclolo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhshihi[i][j] << "  " << rhslohi[i][j] << endl << "  "
                 << rhshilo[i][j] << "  " << rhslolo[i][j] << endl;
      }
   }
}

void cmplx4_define_rhs
 ( int dim, int degp1,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo, int vrblvl )
{
   double acchihi,acclohi,acchilo,acclolo;

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " 
              << funvalrehihi[i][0] << "  " << funvalrelohi[i][0] << endl
              << "  " 
              << funvalrehilo[i][0] << "  " << funvalrelolo[i][0] << endl
              << "  " 
              << funvalimhihi[i][0] << "  " << funvalimlohi[i][0] << endl
              << "  " 
              << funvalimhilo[i][0] << "  " << funvalimlolo[i][0] << endl;
   }
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         // rhsre[i][j] = -(funvalre[j][i] - mbre[j][i]);
         qdf_sub(funvalrehihi[j][i],funvalrelohi[j][i],
                 funvalrehilo[j][i],funvalrelolo[j][i],
                 mbrehihi[j][i],mbrelohi[j][i],
                 mbrehilo[j][i],mbrelolo[j][i],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhsrehihi[i][j] = -acchihi;
         rhsrelohi[i][j] = -acclohi;
         rhsrehilo[i][j] = -acchilo;
         rhsrelolo[i][j] = -acclolo;
         // rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
         qdf_sub(funvalimhihi[j][i],funvalimlohi[j][i],
                 funvalimhilo[j][i],funvalimlolo[j][i],
                 mbimhihi[j][i],mbimlohi[j][i],
                 mbimhilo[j][i],mbimlolo[j][i],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhsimhihi[i][j] = -acchihi;
         rhsimlohi[i][j] = -acclohi;
         rhsimhilo[i][j] = -acchilo;
         rhsimlolo[i][j] = -acclolo;
      }
 
   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsrehihi[i][j] << "  " << rhsrelohi[i][j] << endl
                 << "  "
                 << rhsrehilo[i][j] << "  " << rhsrelolo[i][j] << endl
                 << "  "
                 << rhsimhihi[i][j] << "  " << rhsimlohi[i][j] << endl
                 << "  "
                 << rhsimhilo[i][j] << "  " << rhsimlolo[i][j] << endl;
      }
   }
}

void dbl4_map_evaldiff_output
 ( int dim, int deg,
    double ***outputhihi, double ***outputlohi,
    double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         funvalhihi[i][j] = outputhihi[i][dim][j];
         funvallohi[i][j] = outputlohi[i][dim][j];
         funvalhilo[i][j] = outputhilo[i][dim][j];
         funvallolo[i][j] = outputlolo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the evaluated series :" << endl;

      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalhihi[i][0] << "  "
                            << funvallohi[i][0] << endl
                            << funvalhilo[i][0] << "  "
                            << funvallolo[i][0] << endl;
   }
   // output[i][j][k] is the k-th coefficient in the series
   // the derivative of the i-th polynomial with respect to variable j.

   for(int i=0; i<dim; i++)          // the i-th polynomial
   {
      for(int j=0; j<dim; j++)       // derivative w.r.t. j-th variable
      {
         for(int k=0; k<=deg; k++) 
         {
            jacvalhihi[k][i][j] = outputhihi[i][j][k];
            jacvallohi[k][i][j] = outputlohi[i][j][k];
            jacvalhilo[k][i][j] = outputhilo[i][j][k];
            jacvallolo[k][i][j] = outputlolo[i][j][k];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalhihi[0][i][j] << "  "
                 << jacvallohi[0][i][j] << endl
                 << jacvalhilo[0][i][j] << "  "
                 << jacvallolo[0][i][j] << endl;
      }
   }
}

void cmplx4_map_evaldiff_output
 ( int dim, int deg,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         funvalrehihi[i][j] = outputrehihi[i][dim][j];
         funvalrelohi[i][j] = outputrelohi[i][dim][j];
         funvalrehilo[i][j] = outputrehilo[i][dim][j];
         funvalrelolo[i][j] = outputrelolo[i][dim][j];
         funvalimhihi[i][j] = outputimhihi[i][dim][j];
         funvalimlohi[i][j] = outputimlohi[i][dim][j];
         funvalimhilo[i][j] = outputimhilo[i][dim][j];
         funvalimlolo[i][j] = outputimlolo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the evaluated series :" << endl;

      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalrehihi[i][0] << "  "
                            << funvalrelohi[i][0] << endl
                            << funvalrehilo[i][0] << "  "
                            << funvalrelolo[i][0] << endl
                            << funvalimhihi[i][0] << "  "
                            << funvalimlohi[i][0] << endl
                            << funvalimhilo[i][0] << "  "
                            << funvalimlolo[i][0] << endl;
   }
   // output[i][j][k] is the k-th coefficient in the series
   // the derivative of the i-th polynomial with respect to variable j.

   for(int i=0; i<dim; i++)          // the i-th polynomial
   {
      for(int j=0; j<dim; j++)       // derivative w.r.t. j-th variable
      {
         for(int k=0; k<=deg; k++) 
         {
            jacvalrehihi[k][i][j] = outputrehihi[i][j][k];
            jacvalrelohi[k][i][j] = outputrelohi[i][j][k];
            jacvalrehilo[k][i][j] = outputrehilo[i][j][k];
            jacvalrelolo[k][i][j] = outputrelolo[i][j][k];
            jacvalimhihi[k][i][j] = outputimhihi[i][j][k];
            jacvalimlohi[k][i][j] = outputimlohi[i][j][k];
            jacvalimhilo[k][i][j] = outputimhilo[i][j][k];
            jacvalimlolo[k][i][j] = outputimlolo[i][j][k];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalrehihi[0][i][j] << "  "
                 << jacvalrelohi[0][i][j] << endl
                 << jacvalrehilo[0][i][j] << "  "
                 << jacvalrelolo[0][i][j] << endl
                 << jacvalimhihi[0][i][j] << "  "
                 << jacvalimlohi[0][i][j] << endl
                 << jacvalimhilo[0][i][j] << "  "
                 << jacvalimlolo[0][i][j] << endl;
      }
   }
}
