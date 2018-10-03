// Tests the operations on monomials as defined by polyeq.

#include <iostream>
#include <cstdlib>
#include <string>
#include "complexH.h"
#include "polymon.h"
#include "polyeq.h"
#include "test_utils.h"

using namespace std;

template <class ComplexType, class RealType>
void write_polynomial ( PolyEq<ComplexType,RealType>& p );
/*
 * Writes the terms in p in tableau style format. */

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
void write_polynomial ( PolyEq<ComplexType,RealType>& p )
{
   cout << p.constant.real << "  " << p.constant.imag << " ";
   for(int idx=0; idx<p.dim; idx++) cout << " 0";
   cout << endl;

   for(int idx=0; idx<p.n_mon; idx++)
   {
      cout << p.mon[idx]->coef.real << "  "
           << p.mon[idx]->coef.imag << endl;
   }
}

template <class ComplexType, class RealType>
int test ( int dim )
{
   int nbterms;
   RealType* coefficient = new RealType[2];
   int exponents[dim];

   for(int idx=0; idx<dim; idx++) exponents[idx] = 0;

   cout << "Give the number of terms : "; cin >> nbterms;

   srand(time(NULL));

   PolyEq<ComplexType,RealType> polynomial(dim);
   polynomial.constant = random_complex<ComplexType,RealType>();

   polynomial.n_mon = nbterms-1;

   for(int idx=0; idx<nbterms-1; idx++)
   {
      ComplexType ran = random_complex<ComplexType,RealType>();
      coefficient[0] = ran.real;
      coefficient[1] = ran.imag;

      PolyMon<ComplexType,RealType> *monomial;

      monomial = new PolyMon<ComplexType,RealType>(dim,exponents,coefficient);
   
      polynomial.mon.push_back(monomial);
   }

   cout << endl << "The terms in a random polynomial :" << endl;

   write_polynomial(polynomial);

   return 0;
}
