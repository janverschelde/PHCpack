/* The file prompt_for_setup.cpp contains the definition of a function. */

#include <iostream>
#include "random_polynomials.h"

using namespace std;

void prompt_for_setup
 ( int *seed, int *dim, int *nbr, int *nva, int *pwr, int *deg, int *vrb,
   int *mode )
{
   cout << "Give the seed (0 for time) : "; cin >> *seed;

   cout << "Give the dimension : "; cin >> *dim;

   cout << "Give the variables per monomial (0 for random polynomial) : ";
   cin >> *nva;

   if(*nva > 0)
   {
      cout << "Enter 0 for products, other number of cyclic : ";
      cin >> *nbr;

      if(*nbr == 0)
         *nbr = products_count(*dim,*nva);
      else
         *nbr = *dim;

      cout << "-> number of monomials : " << *nbr << endl;
   }
   else
   {
      cout << "Give the number of terms : "; cin >> *nbr;
   }
   // cout << "Give the largest power of each variable : "; cin >> pwr;
   *pwr = 1;

   cout << "Give the degree of the series : "; cin >> *deg;

   cout << "Give the verbose level : "; cin >> *vrb;

   cout << "Enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> *mode;
}
