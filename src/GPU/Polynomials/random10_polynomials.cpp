// The file random10_polynomials.cpp defines functions specified
// in random10_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "random10_vectors.h"
#include "random10_monomials.h"
#include "random_polynomials.h"

bool make_real10_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrtb, double *cstrix, double *cstrmi, double *cstrrg,
   double *cstrpk, double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, double **cffrrg,
   double **cffrpk, double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_deca_double
         (&cstrtb[i],&cstrix[i],&cstrmi[i],&cstrrg[i],&cstrpk[i],
          &cstltb[i],&cstlix[i],&cstlmi[i],&cstlrg[i],&cstlpk[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real10_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],
                cffrtb[i],cffrix[i],cffrmi[i],cffrrg[i],cffrpk[i],
                cffltb[i],cfflix[i],cfflmi[i],cfflrg[i],cfflpk[i]);

      if(fail) return true;
   }
   return fail;
}

void make_real10_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrtb, double *cstrix, double *cstrmi, double *cstrrg,
   double *cstrpk, double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, double **cffrrg,
   double **cffrpk, double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk )
{
   for(int i=0; i<=deg; i++)
      random_deca_double
         (&cstrtb[i],&cstrix[i],&cstrmi[i],&cstrrg[i],&cstrpk[i],
          &cstltb[i],&cstlix[i],&cstlmi[i],&cstlrg[i],&cstlpk[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_deca_double
            (&cffrtb[k][i],&cffrix[k][i],&cffrmi[k][i],&cffrrg[k][i],
             &cffrpk[k][i],&cffltb[k][i],&cfflix[k][i],&cfflmi[k][i],
             &cfflrg[k][i],&cfflpk[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real10_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrtb, double *cstrix, double *cstrmi, double *cstrrg,
   double *cstrpk, double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, double **cffrrg,
   double **cffrpk, double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk )
{
   for(int i=0; i<=deg; i++)
      random_deca_double
         (&cstrtb[i],&cstrix[i],&cstrmi[i],&cstrrg[i],&cstrpk[i],
          &cstltb[i],&cstlix[i],&cstlmi[i],&cstlrg[i],&cstlpk[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_deca_double
            (&cffrtb[k][i],&cffrix[k][i],&cffrmi[k][i],&cffrrg[k][i],
             &cffrpk[k][i],&cffltb[k][i],&cfflix[k][i],&cfflmi[k][i],
             &cfflrg[k][i],&cfflpk[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
