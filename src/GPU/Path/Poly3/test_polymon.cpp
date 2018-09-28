// testing the operations on monomials

#include <iostream>
#include <string>
#include "complexH.h"
#include "polymon.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( int dim );
/*
 * Runs an interactive test on monomials with
 * as many variables as the value of the dimension dim. */

int main ( void )
{
   cout << endl << "Testing the monomials ..." << endl;

   char choice;

   cout << endl << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   int dimension;

   cout << endl << "Give the dimension : "; cin >> dimension;

   if(choice == '0')
      test<complexH<double>,double>(dimension);
   else if(choice == '1')
      test<complexH<dd_real>,dd_real>(dimension);
   else if(choice == '2')
      test<complexH<qd_real>,qd_real>(dimension);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return 0;
}

template <class ComplexType, class RealType>
int test ( int dim )
{
   RealType* coefficient = new RealType[2];
   int exponents[dim];

   for(int idx=0; idx<dim; idx++)
   {
      cout << "-> give the exponent for variable " << idx+1 << " : ";
      cin >> exponents[idx];
   }
   coefficient[0] = 1.0;
   coefficient[1] = 0.0;

   PolyMon<ComplexType,RealType> monomial
      = PolyMon<ComplexType,RealType>(dim, exponents, coefficient);

   cout << endl << "dim : " << monomial.dim;
   cout << endl << "n_var : " << monomial.n_var;
   cout << endl << "n_base : " << monomial.n_base << endl;

   cout << endl << "The position array of the monomial :";
   for(int idx=0; idx<monomial.n_var; idx++)
      cout << " " << monomial.pos[idx];
   cout << endl;

   cout << endl << "The exponent array of the monomial :";
   for(int idx=0; idx<monomial.n_var; idx++)
      cout << " " << monomial.exp[idx];
   cout << endl;

   cout << endl << "The base position array of the monomial :";
   for(int idx=0; idx<monomial.n_base; idx++)
      cout << " " << monomial.pos_base[idx];
   cout << endl;

   cout << endl << "The base exponent array of the monomial :";
   for(int idx=0; idx<monomial.n_base; idx++)
      cout << " " << monomial.exp_base[idx];
   cout << endl;

   return 0;
}
