// Tests the operations on monomials as defined by polyeq.

#include <iostream>
#include <cstdlib>
#include <string>
#include "complexH.h"
#include "polymon.h"
#include "polyeq.h"
#include "test_utils.h"

using namespace std;

int maximum ( int dim, int* values );
/*
 * Returns the maximum value of dim values. */

template <class ComplexType, class RealType>
void test_evaluation ( PolyEq<ComplexType,RealType> &polynomial );
/*
 * Tests the evaluation and differentation of the polynomial
 * at a random point. */

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

int maximum ( int dim, int* values )
{
   int result = values[0];

   for(int idx=0; idx<dim; idx++)
      if(values[idx] > result) result = values[idx];

   return result;
}

template <class ComplexType, class RealType>
void test_evaluation ( PolyEq<ComplexType,RealType> &polynomial )
{
   const int dim = polynomial.dim;

   cout << endl << "Testing the evaluation ..." << endl;
   ComplexType* point = new ComplexType[dim];
   random_point<ComplexType,RealType>(dim, point);

   ComplexType* derivatives = new ComplexType[dim];
   ComplexType* plainderivs = new ComplexType[dim];

   ComplexType val1 = polynomial.eval(point);
   cout << "the value at a random point : " << val1;

   ComplexType val2 = plain_eval<ComplexType,RealType>(polynomial,point);
   cout << "the value at a random point : " << val2;

   int maxdegrees[dim];
   for(int idx=0; idx<dim; idx++) maxdegrees[idx] = 0;

   polynomial.update_max_deg(maxdegrees);
   cout << "The maximum degrees :";
   for(int idx=0; idx<dim; idx++) cout << " " << maxdegrees[idx];
   int maxdeg = maximum(dim,maxdegrees);
   cout << " max deg : " << maxdeg << endl;

   plain_diff<ComplexType,RealType>(polynomial,point,plainderivs);

   ComplexType val3;

   if(maxdeg == 1)
      val3 = polynomial.eval(point,derivatives);
   else
   {
      ComplexType** powers = new ComplexType*[dim];
      for(int idx=0; idx<dim; idx++)
         powers[idx] = new ComplexType[maxdegrees[idx]+1];
      powertable<ComplexType, RealType>(dim,maxdegrees,point,powers);
      val3 = polynomial.eval(point,derivatives,powers);
   }
   cout << "the value at a random point : " << val3;

   cout << "The derivatives computed plainly and with Speelpenning :"
        << endl;
   for(int idx=0; idx<polynomial.dim; idx++)
   {
      cout << " at " << idx << " : " << plainderivs[idx];
      cout << " at " << idx << " : " << derivatives[idx];
   }
}

template <class ComplexType, class RealType>
int test ( int dim )
{
   int nbterms,expmax;

   cout << "Give the number of terms : "; cin >> nbterms;
   cout << "Give the largest exponent : "; cin >> expmax;

   srand(time(NULL));

   PolyEq<ComplexType,RealType> polynomial
     = random_polynomial<ComplexType,RealType>(dim,nbterms,expmax);

   cout << endl << "The terms in a random polynomial :" << endl;
   write_polynomial<ComplexType,RealType>(polynomial);

   test_evaluation<ComplexType,RealType>(polynomial);

   return 0;
}
