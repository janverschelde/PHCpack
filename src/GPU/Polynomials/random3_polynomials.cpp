// The file random3_polynomials.cpp defines functions specified
// in random3_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "triple_double_functions.h"
#include "random3_vectors.h"
#include "random3_monomials.h"
#include "random_polynomials.h"

bool make_real3_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_triple_double(&csthi[i],&cstmi[i],&cstlo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real3_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],
                                 cffhi[i],cffmi[i],cfflo[i]);
      if(fail) return true;
   }
   return fail;
}

void random_triple_double_complex
 ( double *rehi, double *remi, double *relo,
   double *imhi, double *immi, double *imlo )
{
   double rndhi,rndmi,rndlo,sinhi,sinmi,sinlo;

   random_triple_double(rehi,remi,relo);             // random cos

   tdf_sqr(*rehi,*remi,*relo,&sinhi,&sinmi,&sinlo);  // cos^2(angle)
   tdf_minus(&sinhi,&sinmi,&sinlo);                  // -cos^2(angle)
   tdf_inc_d(&sinhi,&sinmi,&sinlo,1.0);              // 1-cos^2(angle)
   tdf_sqrt(sinhi,sinmi,sinlo,imhi,immi,imlo);       // sin is sqrt
}

bool make_complex3_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_triple_double_complex
         (&cstrehi[i],&cstremi[i],&cstrelo[i],
          &cstimhi[i],&cstimmi[i],&cstimlo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_complex3_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffrehi[i],cffremi[i],cffrelo[i],
                 cffimhi[i],cffimmi[i],cffimlo[i]);

      if(fail) return true;
   }
   return fail;
}

void make_real3_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo )
{
   for(int i=0; i<=deg; i++)
      random_triple_double(&csthi[i],&cstmi[i],&cstlo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_triple_double(&cffhi[k][i],&cffmi[k][i],&cfflo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_complex3_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo )
{
   for(int i=0; i<=deg; i++)
      random_triple_double_complex
         (&cstrehi[i],&cstremi[i],&cstrelo[i],
          &cstimhi[i],&cstimmi[i],&cstimlo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_triple_double_complex
            (&cffrehi[k][i],&cffremi[k][i],&cffrelo[k][i],
             &cffimhi[k][i],&cffimmi[k][i],&cffimlo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real3_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo )
{
   for(int i=0; i<=deg; i++)
      random_triple_double(&csthi[i],&cstmi[i],&cstlo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_triple_double(&cffhi[k][i],&cffmi[k][i],&cfflo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}

void make_complex3_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo )
{
   for(int i=0; i<=deg; i++)
      random_triple_double_complex
         (&cstrehi[i],&cstremi[i],&cstrelo[i],
          &cstimhi[i],&cstimmi[i],&cstimlo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_triple_double_complex
            (&cffrehi[k][i],&cffremi[k][i],&cffrelo[k][i],
             &cffimhi[k][i],&cffimmi[k][i],&cffimlo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
