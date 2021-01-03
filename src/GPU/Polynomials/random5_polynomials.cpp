// The file random5_polynomials.cpp defines functions specified
// in random5_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "random5_vectors.h"
#include "random5_monomials.h"
#include "random_polynomials.h"

bool make_real5_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csttb, double *cstix, double *cstmi, double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, double **cffrg,
   double **cffpk )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_penta_double(&csttb[i],&cstix[i],&cstmi[i],&cstrg[i],&cstpk[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real5_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],
                                 cfftb[i],cffix[i],cffmi[i],
                                 cffrg[i],cffpk[i]);
      if(fail) return true;
   }
   return fail;
}

void make_real5_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csttb, double *cstix, double *cstmi, double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, double **cffrg,
   double **cffpk )
{
   for(int i=0; i<=deg; i++)
      random_penta_double(&csttb[i],&cstix[i],&cstmi[i],&cstrg[i],&cstpk[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_penta_double(&cfftb[k][i],&cffix[k][i],&cffmi[k][i],
                             &cffrg[k][i],&cffpk[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real5_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csttb, double *cstix, double *cstmi, double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, double **cffrg,
   double **cffpk )
{
   for(int i=0; i<=deg; i++)
      random_penta_double(&csttb[i],&cstix[i],&cstmi[i],&cstrg[i],&cstpk[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_penta_double(&cfftb[k][i],&cffix[k][i],&cffmi[k][i],
                             &cffrg[k][i],&cffpk[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
