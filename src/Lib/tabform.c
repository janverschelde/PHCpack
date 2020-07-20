/* This file "tabform.c" contains the definitions of the operations
 * declared in the file "tabform.h". */

#include <stdio.h>
#include "syscon.h"
#include "tabform.h"

int number_of_standard_terms
 ( int neq, int *nbterms, int *nbtsum, int verbose )
{
   int fail,idx,nbt;

   *nbtsum = 0;

   for(idx=1; idx <= neq; idx++)
   {
      fail = syscon_number_of_standard_terms(idx,&nbt);
      if(verbose > 0) printf("  #terms in polynomial %d : %d\n", idx, nbt);
      nbterms[idx-1] = nbt;
      *nbtsum = *nbtsum + nbt;
   }
   if(verbose > 0)
   {
      printf("Total number of terms : %d = %d", *nbtsum, nbterms[0]);
      for(idx=1; idx<neq; idx++) printf(" + %d", nbterms[idx]);
      printf("\n");
   }
   return fail;
}

int standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents,
   int verbose )
{
   int wrk[nvr];  // work space to hold the exponents of a term
   double cff[2]; // holds real and imaginary part of a coefficient
   int cffidx = 0;
   int expidx = 0;
   int fail,idx,idxtrm,k;

   if(verbose > 0) printf("%d\n", neq);
   for(idx=0; idx < neq; idx++)
   {
      if(verbose > 0) printf("%d\n", nbterms[idx]);
      for(idxtrm=1; idxtrm <= nbterms[idx]; idxtrm++)
      {
         fail = syscon_retrieve_standard_term(idx+1,idxtrm,neq,wrk,cff);
         coefficients[cffidx++] = cff[0];
         coefficients[cffidx++] = cff[1];
         if(verbose > 0)
         {
            printf("%22.15e  %22.15e", cff[0], cff[1]);
            for (k=0; k < nvr; k++) printf(" %d", wrk[k]);
            printf("\n");
         }
         for (k=0; k < nvr; k++) exponents[expidx++] = wrk[k];
      }
   }
   return fail;
}

void write_standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents )
{
   int cffidx = 0;
   int expidx = 0;
   int idx,idxtrm,k;

   printf("%d\n", neq);
   for(idx=0; idx < neq; idx++)
   {
      printf("%d\n", nbterms[idx]);
      for(idxtrm=1; idxtrm <= nbterms[idx]; idxtrm++)
      {
         printf("%22.15e  %22.15e",
                coefficients[cffidx], coefficients[cffidx+1]);
         cffidx = cffidx + 2;
         for (k=0; k < nvr; k++) printf(" %d", exponents[expidx++]);
         printf("\n");
      }
   }
}

int store_standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents,
   int verbose )
{
   int fail,idx,nbtsum;

   nbtsum = 0;
   for(idx=0; idx<neq; idx++) nbtsum = nbtsum + nbterms[idx];

   {
      int dim[neq+4];
   
      dim[0] = neq;
      dim[1] = nvr;
      dim[2] = nbtsum;

      for(idx=0; idx<neq; idx++) dim[3+idx] = nbterms[idx];

      dim[neq+3] = verbose;

      fail = _ada_use_c2phc4c(889,dim,exponents,coefficients,0);
   }
   return fail;
}

int load_standard_tableau_dimensions ( int *neq, int *nvr, int *nbt )
{
   int fail,*b;
   double *c;
   int dim[3];

   fail = _ada_use_c2phc4c(890,dim,b,c,0);

   *neq = dim[0];
   *nvr = dim[1];
   *nbt = dim[2];

   return fail;
}
