// The file random5_polynomials.cpp defines functions specified
// in random5_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "penta_double_functions.h"
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

void random_penta_double_complex
 ( double *retb, double *reix, double *remi, double *rerg, double *repk,
   double *imtb, double *imix, double *immi, double *imrg, double *impk )
{
   double rndtb,rndix,rndmi,rndrg,rndpk;
   double sintb,sinix,sinmi,sinrg,sinpk;

   random_penta_double(retb,reix,remi,rerg,repk);     // random cos

   pdf_sqr(*retb,*reix,*remi,*rerg,*repk,
           &sintb,&sinix,&sinmi,&sinrg,&sinpk);       // cos^2(angle)
   pdf_minus(&sintb,&sinix,&sinmi,&sinrg,&sinpk);     // -cos^2(angle)
   pdf_inc_d(&sintb,&sinix,&sinmi,&sinrg,&sinpk,1.0); // 1-cos^2(angle)
   pdf_sqrt(sintb,sinix,sinmi,sinrg,sinpk,
            imtb,imix,immi,imrg,impk);                // sin is sqrt
}

bool make_complex5_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_penta_double_complex
         (&cstretb[i],&cstreix[i],&cstremi[i],&cstrerg[i],&cstrepk[i],
          &cstimtb[i],&cstimix[i],&cstimmi[i],&cstimrg[i],&cstimpk[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_complex5_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffretb[i],cffreix[i],cffremi[i],cffrerg[i],cffrepk[i],
                 cffimtb[i],cffimix[i],cffimmi[i],cffimrg[i],cffimpk[i]);

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

void make_complex5_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk )
{
   for(int i=0; i<=deg; i++)
      random_penta_double_complex
         (&cstretb[i],&cstreix[i],&cstremi[i],&cstrerg[i],&cstrepk[i],
          &cstimtb[i],&cstimix[i],&cstimmi[i],&cstimrg[i],&cstimpk[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_penta_double_complex
            (&cffretb[k][i],&cffreix[k][i],&cffremi[k][i],
             &cffrerg[k][i],&cffrepk[k][i],
             &cffimtb[k][i],&cffimix[k][i],&cffimmi[k][i],
             &cffimrg[k][i],&cffimpk[k][i]);

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

void make_complex5_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi, 
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk )
{
   for(int i=0; i<=deg; i++)
      random_penta_double_complex
         (&cstretb[i],&cstreix[i],&cstremi[i],&cstrerg[i],&cstrepk[i],
          &cstimtb[i],&cstimix[i],&cstimmi[i],&cstimrg[i],&cstimpk[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_penta_double_complex
            (&cffretb[k][i],&cffreix[k][i],&cffremi[k][i],
             &cffrerg[k][i],&cffrepk[k][i],
             &cffimtb[k][i],&cffimix[k][i],&cffimmi[k][i],
             &cffimrg[k][i],&cffimpk[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
