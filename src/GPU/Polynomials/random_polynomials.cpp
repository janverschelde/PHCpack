// The file random_polynomials.cpp defines functions specified
// in random_polynomials.h.

#include <cstdlib>
#include <iostream>
#include "random_numbers.h"
#include "random_monomials.h"

void make_supports ( int dim, int nbr, int *nvr )
{
   int rnd;

   for(int i=0; i<nbr; i++)
   {
      rnd = rand() % dim;  // in range 0..dim-1
      nvr[i] = 1 + rnd;    // in range 1..dim
   }
}

int duplicate_index ( int dim, int nbr, int *nvr, int **idx, int monidx )
{
   int result = -1;

   if(monidx < nbr)
   {
      for(int i=0; i<monidx; i++)
      {
         if(nvr[i] == nvr[monidx])
         {
            result = i;
            for(int j=0; j<nvr[i]; j++)
            {
               if(idx[i][j] != idx[monidx][j]) 
               {
                  result = -1; break;
               }
            }
            if(result != -1) break;
         }
      }
   }
   return result;
}

bool duplicate_supports
 ( int dim, int nbr, int *nvr, int **idx, bool verbose )
{
   bool result = false;
   int dupidx;

   using namespace std;

   for(int i=1; i<nbr; i++)
   {
      dupidx = duplicate_index(dim,nbr,nvr,idx,i);

      if(dupidx != -1)
      {
         if(verbose)
         {
            if(!result) cout << "duplicate indices :";
            cout << " " << dupidx << "==" << i;
         }
         result = true;
      }
      if(!verbose && result) break;
   }
   if(verbose && result) cout << endl;

   return result;
}

bool make_real_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cst, double **cff )
{
   bool fail = false;

   for(int i=0; i<=deg; i++) cst[i] = random_double();

   for(int i=0; i<nbr; i++)
   {
      fail = make_real_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],cff[i]);
      if(fail) return true;
   }
   return fail;
}

int minors_count ( int dim, int nbr )
{
   if(nbr > dim)
      return 0;
   else if(nbr == dim)
      return 1;
   else
   {
      int result = 1;

      for(int i=(dim-nbr+1); i<=dim; i++) result = result*i;

      for(int i=2; i<=nbr; i++) result = result/i;

      return result;
   }
}

void make_exponents
 ( int lvl, int dim, int nbv, int *accu, int *moncnt, int **idx )
{
   if(lvl == nbv)
   {
      for(int i=0; i<nbv; i++) idx[*moncnt][i] = accu[i];
      *moncnt = *moncnt + 1;
   }
   else
   {
      if(lvl == 0)
      {
         for(int i=0; i<dim; i++)
         {
            accu[lvl] = i;
            make_exponents(lvl+1,dim,nbv,accu,moncnt,idx);
         }
      }
      else
      {
         for(int i=accu[lvl-1]+1; i<dim; i++)
         {
            accu[lvl] = i;
            make_exponents(lvl+1,dim,nbv,accu,moncnt,idx);
         }
      }
   }
}

void make_real_minors
 ( int dim, int nbr, int nbv, int deg, int **idx, double *cst, double **cff )
{
   for(int i=0; i<=deg; i++) cst[i] = random_double();

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++) cff[k][i] = random_double();

   int moncnt = 0;
   int *accu = new int[nbv];

   make_exponents(0,dim,nbv,accu,&moncnt,idx);

   free(accu);
}
