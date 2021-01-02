// The file random8_polynomials.cpp defines functions specified
// in random8_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "random8_vectors.h"
#include "random8_monomials.h"
#include "random_polynomials.h"

bool make_real8_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_octo_double
         (&csthihihi[i],&csthilohi[i],&csthihilo[i],&csthilolo[i],
          &cstlohihi[i],&cstlolohi[i],&cstlohilo[i],&cstlololo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real8_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffhihihi[i],cffhilohi[i],cffhihilo[i],cffhilolo[i],
                 cfflohihi[i],cfflolohi[i],cfflohilo[i],cfflololo[i]);

      if(fail) return true;
   }
   return fail;
}

void make_real8_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo )
{
   for(int i=0; i<=deg; i++)
      random_octo_double
         (&csthihihi[i],&csthilohi[i],&csthihilo[i],&csthilolo[i],
          &cstlohihi[i],&cstlolohi[i],&cstlohilo[i],&cstlololo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_octo_double(&cffhihihi[k][i],&cffhilohi[k][i],
                            &cffhihilo[k][i],&cffhilolo[k][i],
                            &cfflohihi[k][i],&cfflolohi[k][i],
                            &cfflohilo[k][i],&cfflololo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real8_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo )
{
   for(int i=0; i<=deg; i++)
      random_octo_double
         (&csthihihi[i],&csthilohi[i],&csthihilo[i],&csthilolo[i],
          &cstlohihi[i],&cstlolohi[i],&cstlohilo[i],&cstlololo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_octo_double(&cffhihihi[k][i],&cffhilohi[k][i],
                            &cffhihilo[k][i],&cffhilolo[k][i],
                            &cfflohihi[k][i],&cfflolohi[k][i],
                            &cfflohilo[k][i],&cfflololo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
