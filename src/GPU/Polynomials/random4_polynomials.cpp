// The file random4_polynomials.cpp defines functions specified
// in random4_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "random4_monomials.h"
#include "random_polynomials.h"

bool make_real4_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_quad_double(&csthihi[i],&cstlohi[i],&csthilo[i],&cstlolo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real4_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],
                                 cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i]);
      if(fail) return true;
   }
   return fail;
}

void random_quad_double_complex
 ( double *rehihi, double *relohi, double *rehilo, double *relolo,
   double *imhihi, double *imlohi, double *imhilo, double *imlolo )
{
   double rndhihi,rndlohi,rndhilo,rndlolo;
   double sinhihi,sinlohi,sinhilo,sinlolo;

   random_quad_double(rehihi,relohi,rehilo,relolo);    // random cos

   qdf_sqr(*rehihi,*relohi,*rehilo,*relolo,
           &sinhihi,&sinlohi,&sinhilo,&sinlolo);       // cos^2(angle)
   qdf_minus(&sinhihi,&sinlohi,&sinhilo,&sinlolo);     // -cos^2(angle)
   qdf_inc_d(&sinhihi,&sinlohi,&sinhilo,&sinlolo,1.0); // 1-cos^2(angle)
   qdf_sqrt(sinhihi,sinlohi,sinhilo,sinlolo,
            imhihi,imlohi,imhilo,imlolo);              // sin is sqrt
}

bool make_complex4_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehihi, double *cstrelohi, double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi, double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_quad_double_complex
         (&cstrehihi[i],&cstrelohi[i],&cstrehilo[i],&cstrelolo[i],
          &cstimhihi[i],&cstimlohi[i],&cstimhilo[i],&cstimlolo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_complex4_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
                 cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i]);

      if(fail) return true;
   }
   return fail;
}

void make_real4_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   for(int i=0; i<=deg; i++)
      random_quad_double(&csthihi[i],&cstlohi[i],&csthilo[i],&cstlolo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_quad_double(&cffhihi[k][i],&cfflohi[k][i],
                            &cffhilo[k][i],&cfflolo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_complex4_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehihi, double *cstrelohi, double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi, double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo )
{
   for(int i=0; i<=deg; i++)
      random_quad_double_complex
         (&cstrehihi[i],&cstrelohi[i],&cstrehilo[i],&cstrelolo[i],
          &cstimhihi[i],&cstimlohi[i],&cstimhilo[i],&cstimlolo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_quad_double_complex
            (&cffrehihi[k][i],&cffrelohi[k][i],
             &cffrehilo[k][i],&cffrelolo[k][i],
             &cffimhihi[k][i],&cffimlohi[k][i],
             &cffimhilo[k][i],&cffimlolo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}


void make_real4_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   for(int i=0; i<=deg; i++)
      random_quad_double(&csthihi[i],&cstlohi[i],&csthilo[i],&cstlolo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_quad_double(&cffhihi[k][i],&cfflohi[k][i],
                            &cffhilo[k][i],&cfflolo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}

void make_complex4_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehihi, double *cstrelohi, double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi, double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo )
{
   for(int i=0; i<=deg; i++)
      random_quad_double_complex
         (&cstrehihi[i],&cstrelohi[i],&cstrehilo[i],&cstrelolo[i],
          &cstimhihi[i],&cstimlohi[i],&cstimhilo[i],&cstimlolo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_quad_double_complex
            (&cffrehihi[k][i],&cffrelohi[k][i],
             &cffrehilo[k][i],&cffrelolo[k][i],
             &cffimhihi[k][i],&cffimlohi[k][i],
             &cffimhilo[k][i],&cffimlolo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
