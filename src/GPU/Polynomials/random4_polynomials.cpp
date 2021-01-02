// The file random4_polynomials.cpp defines functions specified
// in random4_polynomials.h.

#include <cstdlib>
#include <iostream>
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
