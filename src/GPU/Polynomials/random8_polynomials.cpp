// The file random8_polynomials.cpp defines functions specified
// in random8_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "random8_vectors.h"
#include "random8_monomials.h"
#include "random_polynomials.h"

bool make_real8_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_octo_double
         (&csthihihi[i],&cstlohihi[i],&csthilohi[i],&cstlolohi[i],
          &csthihilo[i],&cstlohilo[i],&csthilolo[i],&cstlololo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_real8_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
                 cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i]);

      if(fail) return true;
   }
   return fail;
}

bool make_cmplx8_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo )
{
   bool fail = false;

   for(int i=0; i<=deg; i++)
      random_octo_complex
         (&cstrehihihi[i],&cstrelohihi[i],&cstrehilohi[i],&cstrelolohi[i],
          &cstrehihilo[i],&cstrelohilo[i],&cstrehilolo[i],&cstrelololo[i],
          &cstimhihihi[i],&cstimlohihi[i],&cstimhilohi[i],&cstimlolohi[i],
          &cstimhihilo[i],&cstimlohilo[i],&cstimhilolo[i],&cstimlololo[i]);

   for(int i=0; i<nbr; i++)
   {
      fail = make_complex8_monomial
                (dim,nvr[i],pwr,deg,idx[i],exp[i],
                 cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
                 cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
                 cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
                 cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i]);

      if(fail) return true;
   }
   return fail;
}

void make_real8_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   for(int i=0; i<=deg; i++)
      random_octo_double
         (&csthihihi[i],&cstlohihi[i],&csthilohi[i],&cstlolohi[i],
          &csthihilo[i],&cstlohilo[i],&csthilolo[i],&cstlololo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_octo_double(&cffhihihi[k][i],&cfflohihi[k][i],
                            &cffhilohi[k][i],&cfflolohi[k][i],
                            &cffhihilo[k][i],&cfflohilo[k][i],
                            &cffhilolo[k][i],&cfflololo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_cmplx8_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo )
{
   for(int i=0; i<=deg; i++)
      random_octo_complex
         (&cstrehihihi[i],&cstrelohihi[i],&cstrehilohi[i],&cstrelolohi[i],
          &cstrehihilo[i],&cstrelohilo[i],&cstrehilolo[i],&cstrelololo[i],
          &cstimhihihi[i],&cstimlohihi[i],&cstimhilohi[i],&cstimlolohi[i],
          &cstimhihilo[i],&cstimlohilo[i],&cstimhilolo[i],&cstimlololo[i]);

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
         random_octo_complex(&cffrehihihi[k][i],&cffrelohihi[k][i],
                             &cffrehilohi[k][i],&cffrelolohi[k][i],
                             &cffrehihilo[k][i],&cffrelohilo[k][i],
                             &cffrehilolo[k][i],&cffrelololo[k][i],
                             &cffimhihihi[k][i],&cffimlohihi[k][i],
                             &cffimhilohi[k][i],&cffimlolohi[k][i],
                             &cffimhihilo[k][i],&cffimlohilo[k][i],
                             &cffimhilolo[k][i],&cffimlololo[k][i]);

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real8_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   for(int i=0; i<=deg; i++)
      random_octo_double
         (&csthihihi[i],&cstlohihi[i],&csthilohi[i],&cstlolohi[i],
          &csthihilo[i],&cstlohilo[i],&csthilolo[i],&cstlololo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_octo_double(&cffhihihi[k][i],&cfflohihi[k][i],
                            &cffhilolo[k][i],&cfflolohi[k][i],
                            &cffhihilo[k][i],&cfflohilo[k][i],
                            &cffhilolo[k][i],&cfflololo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}

void make_cmplx8_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo )
{
   for(int i=0; i<=deg; i++)
      random_octo_complex
         (&cstrehihihi[i],&cstrelohihi[i],&cstrehilohi[i],&cstrelolohi[i],
          &cstrehihilo[i],&cstrelohilo[i],&cstrehilolo[i],&cstrelololo[i],
          &cstimhihihi[i],&cstimlohihi[i],&cstimhilohi[i],&cstimlolohi[i],
          &cstimhihilo[i],&cstimlohilo[i],&cstimhilolo[i],&cstimlololo[i]);

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
         random_octo_complex(&cffrehihihi[k][i],&cffrelohihi[k][i],
                             &cffrehilohi[k][i],&cffrelolohi[k][i],
                             &cffrehihilo[k][i],&cffrelohilo[k][i],
                             &cffrehilolo[k][i],&cffrelololo[k][i],
                             &cffimhihihi[k][i],&cffimlohihi[k][i],
                             &cffimhilohi[k][i],&cffimlolohi[k][i],
                             &cffimhihilo[k][i],&cffimlohilo[k][i],
                             &cffimhilolo[k][i],&cffimlololo[k][i]);

   make_cyclic_exponents(dim,nva,idx);
}
