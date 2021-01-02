// The file random3_polynomials.cpp defines functions specified
// in random3_polynomials.h.

#include <cstdlib>
#include <iostream>
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
