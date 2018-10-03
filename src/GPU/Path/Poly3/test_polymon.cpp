// Tests the operations on monomials as defined by polymon.

#include <iostream>
#include <cstdlib>
#include <string>
#include "complexH.h"
#include "polymon.h"
#include "test_utils.h"

using namespace std;

template <class ComplexType, class RealType>
void test_evaluation ( PolyMon<ComplexType,RealType>& m );
/*
 * Test the application of the straightforward and Speelpenning algorithm 
 * to evaluate m at x. */

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
void test_evaluation ( PolyMon<ComplexType,RealType>& m )
{
   print_data<ComplexType,RealType>(m);

   ComplexType* point = new ComplexType[m.dim];
   random_point<ComplexType,RealType>(m.dim, point);

   ComplexType val1 = m.eval(point);
   cout << "The value at a random point : " << val1;

   ComplexType val2 = plain_eval<ComplexType,RealType>(m, point);
   cout << "The value at a random point : " << val2;

   ComplexType* derivatives = new ComplexType[m.n_var];
   ComplexType* plainderivs = new ComplexType[m.n_var];

   plain_diff<ComplexType,RealType>(m, point, plainderivs);

   if(m.n_base == 0)
   {
      ComplexType val3 = m.speel(point, derivatives);
      cout << "The value at a random point : " << val3 << endl;
   }
   else
   {
      ComplexType** powers = new ComplexType*[m.dim];
      for(int idx=0; idx<m.n_var; idx++)
      {
         int varidx = m.pos[idx]; // index of the variable
         int size = m.exp[idx]+1;
         powers[varidx] = new ComplexType[size];
         powers[varidx][0] = point[varidx];
         for(int powidx=1; powidx<size; powidx++)
            powers[varidx][powidx] = powers[varidx][powidx-1]*point[varidx];
      }
      ComplexType valbase = m.eval_base(point, powers);
      ComplexType val4 = m.speel_with_base(point, derivatives, valbase);    
      cout << "The value at a random point : " << val4 << endl;
   }
   cout << "The derivatives computed plainly and with Speelpenning :" << endl;
   for(int idx=0; idx<m.n_var; idx++)
   {
      cout << " at " << m.pos[idx] << " : " << plainderivs[idx];
      cout << " at " << m.pos[idx] << " : " << derivatives[idx];
   }
}

template <class ComplexType, class RealType>
int test ( int dim )
{
   cout << "Test a random monomial ? (y/n) ";
   char answer; cin >> answer;

   if(answer != 'y')
   {
      PolyMon<ComplexType,RealType> monomial
         = prompted_monomial<ComplexType,RealType>(dim);

      test_evaluation<ComplexType,RealType>(monomial);
   }
   else
   {
      srand(time(NULL));

      cout << "Give the largest exponent in a monomial : ";
      int expmax; cin >> expmax;

      PolyMon<ComplexType,RealType> monomial
         = random_monomial<ComplexType,RealType>(dim, expmax);

      test_evaluation<ComplexType,RealType>(monomial);
   }
   return 0;
}
