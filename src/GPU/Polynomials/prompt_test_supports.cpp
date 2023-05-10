// The file prompt_test_supports.cpp defines the functions specified in
// the file prompt_test_supports.h.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "random_polynomials.h"
#include "prompt_test_supports.h"

using namespace std;

void prompt_testpoly_dimensions ( int *seed, int *dim, int *nva, int *nbr )
{
   cout << "Give the seed (0 for time) : ";
   int s; cin >> s;

   cout << "Give the dimension : "; cin >> *dim;

   cout << "Give the variables per monomial (0 for random polynomial) : ";
   cin >> *nva;

   if(*nva <= 0)
   {
      cout << "Give the number of terms : "; cin >> *nbr;
   }
   else
   {
      cout << "Enter 0 for products, other number of cyclic : ";
      cin >> *nbr;

      if(*nbr == 0)
         *nbr = products_count(*dim,*nva);
      else
         *nbr = *dim;

      cout << "-> number of monomials : " << *nbr << endl;
   }
   if(s != 0)
   {
      srand(s);
      *seed = s;
   }
   else
   {
      const int timevalue = time(NULL); // for a random seed
      srand(timevalue);
      *seed = timevalue;
   }
   cout << "  seed used : " << *seed << endl;
}

void make_test_supports ( int dim, int nva, int nbr, int *nvr, int **idx )
{
   const int deg = 0;
   const int pwr = 1;

   double *cst = new double[deg+1]; // constant coefficient series
   double **cff = new double*[nbr]; // coefficient series of terms
   for(int i=0; i<nbr; i++) cff[i] = new double[deg+1];

   if(nva == 0) make_supports(dim,nbr,nvr); // define supports of polynomial

   if(nva == 0)
   {
      for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   }
   else 
   {
      for(int i=0; i<nbr; i++)
      {
         idx[i] = new int[nva];
         nvr[i] = nva;
      }
   }
   if(nva > 0)
   {
      if(nbr == dim)
         make_real_cyclic(dim,nva,deg,idx,cst,cff);
      else
         make_real_products(dim,nbr,nva,deg,idx,cst,cff);
   }
   else
   {
      int **exp = new int*[nbr];  // exponents of the variables
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,cst,cff);
   }
   for(int i=0; i<nbr; i++)
   {
      cout << "Indices of monomial " << i << " :";
      for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j]; cout << endl;
   }
   bool dup = duplicate_supports(dim,nbr,nvr,idx,false);
   if(dup)
      cout << "Duplicate supports found." << endl;
   else
      cout << "No duplicate supports found." << endl;
}
