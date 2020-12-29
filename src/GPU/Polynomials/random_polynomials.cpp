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
