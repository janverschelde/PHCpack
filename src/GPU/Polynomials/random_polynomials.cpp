// The file random_polynomials.cpp defines functions specified
// in random_polynomials.h.

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random_monomials.h"

using namespace std;

void make_supports ( int dim, int nbr, int *nvr )
{
   int rnd;

   for(int i=0; i<nbr; i++)
   {
      rnd = rand() % dim;  // in range 0..dim-1
      nvr[i] = 1 + rnd;    // in range 1..dim
   }
}

void read_supports ( int dim, int nbr, int *nvr )
{
   int inpnvr;

   for(int i=0; i<nbr; i++)
   {
      do
      {
         cout << "Give number of variables in monomial " << i << " : ";
         cin >> inpnvr;
         if(inpnvr < 1)
            cout << "-> entered " << inpnvr << " < 1, retry" << endl; 
         else if (inpnvr > dim)
            cout << "-> entered " << inpnvr << " > " << dim
                 << ", retry" << endl; 
      }
      while((inpnvr < 1) || (inpnvr > dim));
      nvr[i] = inpnvr;    // in range 1..dim
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

bool make_complex_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstre, double *cstim, double **cffre, double **cffim )
{
   bool fail = false;
   double rnd;

   for(int i=0; i<=deg; i++)
   {
      rnd = random_angle();        
      cstre[i] = cos(rnd);
      cstim[i] = sin(rnd);
   }
   for(int i=0; i<nbr; i++)
   {
      fail = make_complex_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],
                                   cffre[i],cffim[i]);
      if(fail) return true;
   }
   return fail;
}

int products_count ( int dim, int nbr )
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

void make_product_exponents
 ( int lvl, int dim, int nva, int *accu, int *moncnt, int **idx )
{
   if(lvl == nva)
   {
      for(int i=0; i<nva; i++) idx[*moncnt][i] = accu[i];
      *moncnt = *moncnt + 1;
   }
   else
   {
      if(lvl == 0)
      {
         for(int i=0; i<dim; i++)
         {
            accu[lvl] = i;
            make_product_exponents(lvl+1,dim,nva,accu,moncnt,idx);
         }
      }
      else
      {
         for(int i=accu[lvl-1]+1; i<dim; i++)
         {
            accu[lvl] = i;
            make_product_exponents(lvl+1,dim,nva,accu,moncnt,idx);
         }
      }
   }
}

void insert_sort ( int dim, int *data )
{
   int smallest,tmp;

   for(int i=0; i<dim; i++)
   {
      smallest = i;

      for(int j=i+1; j<dim; j++) 
      {
         if(data[j] < data[smallest]) smallest = j;
      }
      if(smallest != i)
      {
         tmp = data[i];
         data[i] = data[smallest];
         data[smallest] = tmp;
      }
   }
}

void make_cyclic_exponents ( int dim, int nva, int **idx )
{
   for(int i=0; i<nva; i++) idx[0][i] = i;

   for(int k=1; k<dim; k++)
   {
      if(idx[k-1][nva-1] < dim-1)
         idx[k][nva-1] = idx[k-1][nva-1] + 1;
      else
         idx[k][nva-1] = 0;
     
      for(int i=1; i<nva; i++) idx[k][i-1] = idx[k-1][i];
   }
   for(int k=1; k<dim; k++) insert_sort(nva,idx[k]);
}

void make_real_products
 ( int dim, int nbr, int nva, int deg, int **idx, double *cst, double **cff )
{
   for(int i=0; i<=deg; i++) cst[i] = random_double();

   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++) cff[k][i] = random_double();

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_complex_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim )
{
   double rnd;

   for(int i=0; i<=deg; i++)
   {
      rnd = random_angle();        
      cstre[i] = cos(rnd);
      cstim[i] = sin(rnd);
   }
   for(int k=0; k<nbr; k++)
      for(int i=0; i<=deg; i++)
      {
         rnd = random_angle();
         cffre[k][i] = cos(rnd);
         cffim[k][i] = sin(rnd);
      }

   int moncnt = 0;
   int *accu = new int[nva];

   make_product_exponents(0,dim,nva,accu,&moncnt,idx);

   free(accu);
}

void make_real_cyclic
 ( int dim, int nva, int deg, int **idx, double *cst, double **cff )
{
   for(int i=0; i<=deg; i++) cst[i] = random_double();

   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++) cff[k][i] = random_double();

   make_cyclic_exponents(dim,nva,idx);
}

void make_complex_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim )
{
   double rnd;

   for(int i=0; i<=deg; i++)
   {
      rnd = random_angle();        
      cstre[i] = cos(rnd);
      cstim[i] = sin(rnd);
   }
   for(int k=0; k<dim; k++)
      for(int i=0; i<=deg; i++)
      {
         rnd = random_angle();
         cffre[k][i] = cos(rnd);
         cffim[k][i] = sin(rnd);
      }

   make_cyclic_exponents(dim,nva,idx);
}
