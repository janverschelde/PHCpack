// Collects the addition jobs to evaluate and differentiate
// one polynomial in several variables.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "random_polynomials.h"
#include "addition_job.h"
#include "addition_jobs.h"

using namespace std;

void differential_index_count
 ( int dim, int nbr, int *nvr, int **idx, int *cnt, bool verbose );
/*
 * DESCRIPTION :
 *   Counts the number of monomials where the variables appear.
 *
 * ON ENTRY :
 *   dim      the total number of variables;
 *   nbr      number of monomials;
 *   nvr      nbr integers, nvr[k] counts the number of variables
 *            in the k-th monomial;
 *   idx      array of nbr support vectors,
 *            idx[k] has nvr[k] integers,
 *            idx[k][i] is the index of the i-th variable in monomial k;
 *   cnt      space for dim integers;
 *   verbose  indicates if output needs to be written during the count.
 *
 * ON RETURN :
 *   cnt      cnt[k] counts the number of monomials that contain k. */ 

void differential_indices
 ( int dim, int nbr, int *nvr, int **idx, int *cnt, int **difidx,
   bool verbose );
/*
 * DESCRIPTION :
 *   Defines the rearrangment of the monomial indices according
 *   to which monomials the variables appear.
 *
 * ON ENTRY :
 *   dim      the total number of variables;
 *   nbr      number of monomials;
 *   nvr      nbr integers, nvr[k] counts the number of variables
 *            in the k-th monomial;
 *   idx      array of nbr support vectors,
 *            idx[k] has nvr[k] integers,
 *            idx[k][i] is the index of the i-th variable in monomial k;
 *   cnt      dim counts, cnt[k] counts the number of monomials
 *            that contain k;
 *   difidx   array of dim integer vectors,
 *            difidx[k] has space for cnt[k]+1 numbers;
 *   verbose  indicates if output needs to be written;
 *
 * ON RETURN :
 *   cnt      cnt[k] is reduced by one if there is a monomial where
 *            only the variable k occurs;
 *   difidx   differential indices,
 *            difidx[k] lists all monomial indices that contain variable k,
 *            difidx[k][0] == -1 if monomial coefficients play no role,
 *                         == i, if cff[i] contributes to derivative k,
 *            difidx[k][i] is the index of the i-th monomial that has k. */

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;


   cout << "Give the number of terms : ";
   int nbr; cin >> nbr;

   int seedused;

   if(seed != 0)
   {
      srand(seed);
      seedused = seed;
   }
   else
   {
      const int timevalue = time(NULL); // for a random seed
      srand(timevalue);
      seedused = timevalue;
   }
   const int deg = 0;
   const int pwr = 1;

   double *cst = new double[deg+1]; // constant coefficient series
   double **cff = new double*[nbr]; // coefficient series of terms
   for(int i=0; i<nbr; i++) cff[i] = new double[deg+1];
   int *nvr = new int[nbr]; // number of variables in each monomial

   make_supports(dim,nbr,nvr); // define supports of polynomial

   int **idx = new int*[nbr];  // indices of variables in monomials
   for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   int **exp = new int*[nbr];  // exponents of the variables
   for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

   bool fail = make_real_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,cst,cff);

   for(int i=0; i<nbr; i++)
   {
      cout << "Indices of monomial " << i << " :";
      for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j]; cout << endl;
   }
   int *idxcnt = new int[dim];
   differential_index_count(dim,nbr,nvr,idx,idxcnt,true);
   cout << "Index count :";
   for(int i=0; i<dim; i++) cout << " " << idxcnt[i];
   cout << endl;
   int **difidx = new int*[dim];
   for(int i=0; i<dim; i++) difidx[i] = new int[idxcnt[i]+1];
   differential_indices(dim,nbr,nvr,idx,idxcnt,difidx,true);
   cout << "The differential indices :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "variable " << i << " :";
      for(int j=0; j<=idxcnt[i]; j++) cout << " " << difidx[i][j];
      cout << endl;
   }
/*
   AdditionJobs jobs(nbr); // initialize with the number of monomials

   jobs.make(nbr,nvr,true);

   cout << "number of addition jobs : " << jobs.get_count() << endl;
   cout << "number of layers : " << jobs.get_depth() << endl;
   cout << "frequency of layer counts :" << endl;
   int checksum = 0;
   for(int i=0; i<jobs.get_depth(); i++)
   {
      cout << i << " : " << jobs.get_layer_count(i) << endl;
      checksum = checksum + jobs.get_layer_count(i); 
   }
   cout << "layer count sum : " << checksum << endl;

   for(int k=0; k<jobs.get_depth(); k++)
   {
      cout << "jobs at layer " << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
         cout << jobs.get_job(k,i) << endl;
   }
 */
   cout << "seed used : " << seedused << endl;

   return 0;
}

void differential_index_count
 ( int dim, int nbr, int *nvr, int **idx, int *cnt, bool verbose )
{
   for(int i=0; i<dim; i++)
   {
      if(verbose) cout << "Variable " <<  i << " occurs in monomials"; 

      cnt[i] = 0;
      for(int j=0; j<nbr; j++)
      {
         for(int k=0; k<nvr[j]; k++)
            if(idx[j][k] == i)
            {
               if(verbose)
               {
                  cout << " " << j;
                  if(nvr[j] == 1) cout << "(cff!)";
               }
               cnt[i] = cnt[i] + 1; break;
            }
      }
      if(verbose) cout << endl;
   }
}

void differential_indices
 ( int dim, int nbr, int *nvr, int **idx, int *cnt, int **difidx,
   bool verbose )
{
   int pos;

   for(int i=0; i<dim; i++)
   {
      if(verbose) cout << "Variable " <<  i << " occurs in monomials"; 
      
      difidx[i][0] = -1;
      pos = 1;
      for(int j=0; j<nbr; j++)
      {
         for(int k=0; k<nvr[j]; k++)
            if(idx[j][k] == i)
            {
               if(verbose)
               {
                  cout << " " << j;
                  if(nvr[j] == 1) cout << "(cff!)";
               }
               if(nvr[j] == 1)
               {
                  difidx[i][0] = j;
                  cnt[i] = cnt[i] - 1;
               }
               else
               {
                  difidx[i][pos++] = j;
               }
               break;
            }
      }
      if(verbose) cout << endl;
   }
}
