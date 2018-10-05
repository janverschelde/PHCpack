// Tests the operations on polynomial systems as defined by polysys.

#include <iostream>
#include <cstdlib>
#include <string>
#include "complexH.h"
#include "polymon.h"
#include "polyeq.h"
#include "polysys.h"
#include "test_utils.h"

using namespace std;

template <class ComplexType, class RealType>
PolySys<ComplexType,RealType> random_polynomials 
 ( int dim, int nbterms, int expmax );
/*
 * Returns a random polynomial system in a space with dim variables.
 * Every polynomial in the system on return has as many terms as nbterms.
 * The largest exponent is defined by the value of expmax. */

template <class ComplexType, class RealType>
void write_polynomials ( PolySys<ComplexType,RealType>& p );
/*
 * Writes the polynomials in p in tableau style format. */

template <class ComplexType, class RealType>
ComplexType* plain_eval
 ( PolySys<ComplexType,RealType> &polsys, ComplexType *x );
/*
 * Applies the straightforward algorithm to evaluate polsys at x. */

template <class ComplexType, class RealType>
void test_evaluation ( PolySys<ComplexType,RealType> &polsys );
/*
 * Tests the evaluation and differentation of the polynomial system
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

template <class ComplexType, class RealType>
PolySys<ComplexType,RealType> random_polynomials 
 ( int dim, int nbterms, int expmax )
{
   PolySys<ComplexType,RealType> result;

   result.n_eq = dim;
   result.dim = dim;

   for(int idx=0; idx<dim; idx++)
   {
      PolyEq<ComplexType,RealType>* pol
         = new PolyEq<ComplexType,RealType>
                  (random_polynomial<ComplexType,RealType>(dim,nbterms,expmax));
      result.eq.push_back(pol);
   }

   return result;
}

template <class ComplexType, class RealType>
void write_polynomials ( PolySys<ComplexType,RealType>& p )
{
   cout << "The dimension : " << p.dim << endl;

   for(int idx=0; idx<p.n_eq; idx++)
   {
      write_polynomial<ComplexType,RealType>(*(p.eq[idx]));
   }
}

template <class ComplexType, class RealType>
ComplexType* plain_eval
 ( PolySys<ComplexType,RealType> &polsys, ComplexType *x )
{
   return x;
}

template <class ComplexType, class RealType>
void test_evaluation ( PolySys<ComplexType,RealType> &polsys )
{
   const int neq = polsys.n_eq;
   const int dim = polsys.dim;

   cout << endl << "Testing the evaluation ..." << endl;
   ComplexType* point = new ComplexType[dim];
   random_point<ComplexType,RealType>(dim, point);

   ComplexType* derivatives = new ComplexType[dim];
   ComplexType* plainderivs = new ComplexType[dim];
/*
   ComplexType val1 = polsys.eval(point);
   cout << "the value at a random point : " << val1;

   ComplexType val2 = plain_eval<ComplexType,RealType>(polsys,point);
   cout << "the value at a random point : " << val2;
 */
}

template <class ComplexType, class RealType>
int test ( int dim )
{
   int nbterms,expmax;

   cout << "Give the number of terms : "; cin >> nbterms;
   cout << "Give the largest exponent : "; cin >> expmax;

   srand(time(NULL));

   PolySys<ComplexType,RealType> polynomials
     = random_polynomials<ComplexType,RealType>(dim,nbterms,expmax);

   cout << endl << "The terms in a random polynomial system :" << endl;
   write_polynomials<ComplexType,RealType>(polynomials);

   test_evaluation<ComplexType,RealType>(polynomials);

   return 0;
}
