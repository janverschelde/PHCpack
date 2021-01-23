// The file random2_polynomials.cpp defines functions specified
// in random2_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "double_double_functions.h"
#include "random2_vectors.h"
#include "random2_monomials.h"
#include "random_polynomials.h"

bool make_real2_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthi, double *cstlo, double **cffhi, double **cfflo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++) random_double_double(&csthi[i],&cstlo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real2_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],
                                 cffhi[i],cfflo[i]);
      if(fail) return true;
   }
   return fail;
}

void random_double_double_complex
 ( double *rehi, double *relo, double *imhi, double *imlo )
{
   double rndhi,rndlo,sinhi,sinlo;

   random_double_double(rehi,relo);           // random cos

   ddf_sqr(*rehi,*relo,&sinhi,&sinlo);        // cos^2(angle)
   ddf_minus(&sinhi,&sinlo);                  // -cos^2(angle)
   ddf_inc_d(&sinhi,&sinlo,1.0);              // 1-cos^2(angle)
   ddf_sqrt(sinhi,sinlo,imhi,imlo);           // sin is sqrt
}

bool make_complex2_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_double_double_complex
         (&cstrehi[i],&cstrelo[i],&cstimhi[i],&cstimlo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_complex2_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i]);

      if(fail) return true;
   }
   return fail;
}

void make_real2_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo )
{
   for(int i=0; i<=deg; i++) random_double_double(&csthi[i],&cstlo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_double_double(&cffhi[k][i],&cfflo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_complex2_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo )
{
   for(int i=0; i<=deg; i++)
      random_double_double_complex
         (&cstrehi[i],&cstrelo[i],&cstimhi[i],&cstimlo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_double_double_complex
            (&cffrehi[k][i],&cffrelo[k][i],&cffimhi[k][i],&cffimlo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real2_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo )
{
   for(int i=0; i<=deg; i++) random_double_double(&csthi[i],&cstlo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_double_double(&cffhi[k][i],&cfflo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}

void make_complex2_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo )
{
   for(int i=0; i<=deg; i++)
      random_double_double_complex
         (&cstrehi[i],&cstrelo[i],&cstimhi[i],&cstimlo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_double_double_complex
            (&cffrehi[k][i],&cffrelo[k][i],&cffimhi[k][i],&cffimlo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
