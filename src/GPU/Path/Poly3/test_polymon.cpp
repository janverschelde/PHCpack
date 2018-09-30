// Tests the operations on monomials as defined by polymon.

#include <iostream>
#include <cstdlib>
#include <string>
#include "complexH.h"
#include "polymon.h"

using namespace std;

template <class ComplexType, class RealType>
ComplexType random_complex ();
/*
 * Returns a random complex number. */

template <class ComplexType, class RealType>
void random_point ( int dim, ComplexType* coordinates );
/*
 * Generates the coordinates of a random complex point.
 * Space must be allocated for the dim many coordinates. */

template <class ComplexType, class RealType>
ComplexType plain_eval ( PolyMon<ComplexType,RealType>& m, ComplexType *x );
/*
 * Applies the straightforward algorithm to evaluate m at x. */

template <class ComplexType, class RealType>
void test_plain_eval ( PolyMon<ComplexType,RealType>& m );
/*
 * Test the application of the 
 * straightforward algorithm to evaluate m at x. */

template <class ComplexType, class RealType>
void print_data ( PolyMon<ComplexType,RealType>& monomial );
/*
 * Prints the content of the data stored in monomial. */

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
ComplexType random_complex ()
{
   double angle = 3.14*((double) rand())/RAND_MAX;
   double re_part = cos(angle);
   double im_part = sin(angle);

   RealType real_part = RealType(re_part);
   RealType imag_part = RealType(im_part);

   ComplexType result = ComplexType(real_part, imag_part);

   return result;
}

template <class ComplexType, class RealType>
void random_point ( int dim, ComplexType* coordinates )
{
   for(int idx=0; idx<dim; idx++)
      coordinates[idx] = random_complex<ComplexType,RealType>();
}

template <class ComplexType, class RealType>
ComplexType plain_eval ( PolyMon<ComplexType,RealType>& m, ComplexType* x )
{
   ComplexType result = m.coef;

   for(int varidx=0; varidx<m.n_var; varidx++)
   {
      int idx = m.pos[varidx];

      for(int expidx=0; expidx<m.exp[varidx]; expidx++)
         result = result*x[idx];
   }
   return result;
}

template <class ComplexType, class RealType>
void print_data ( PolyMon<ComplexType,RealType>& monomial )
{
   cout << endl << "The coefficient : " << monomial.coef << endl;

   cout << "dim : " << monomial.dim;
   cout << "  n_var : " << monomial.n_var;
   cout << "  n_base : " << monomial.n_base << endl;

   cout << endl << "The position array of the monomial :";
   for(int idx=0; idx<monomial.n_var; idx++)
      cout << " " << monomial.pos[idx];
   cout << endl;

   cout << "The exponent array of the monomial :";
   for(int idx=0; idx<monomial.n_var; idx++)
      cout << " " << monomial.exp[idx];
   cout << endl;

   cout << endl << "The base position array of the monomial :";
   for(int idx=0; idx<monomial.n_base; idx++)
      cout << " " << monomial.pos_base[idx];
   cout << endl;

   cout << "The base exponent array of the monomial :";
   for(int idx=0; idx<monomial.n_base; idx++)
      cout << " " << monomial.exp_base[idx];
   cout << endl << endl;
}

template <class ComplexType, class RealType>
void test_plain_eval ( PolyMon<ComplexType,RealType>& m )
{
   ComplexType* point = new ComplexType[m.dim];
   random_point<ComplexType,RealType>(m.dim, point);

   ComplexType val1 = m.eval(point);
   cout << "The value at a random point : " << val1;

   ComplexType val2 = plain_eval<ComplexType,RealType>(m, point);
   cout << "The value at a random point : " << val2 << endl;
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
   srand(time(NULL));
   ComplexType ran = random_complex<ComplexType,RealType>();

   coefficient[0] = ran.real;
   coefficient[1] = ran.imag;

   PolyMon<ComplexType,RealType> monomial
      = PolyMon<ComplexType,RealType>(dim, exponents, coefficient);

   print_data<ComplexType,RealType>(monomial);
   test_plain_eval<ComplexType,RealType>(monomial);

   return 0;
}
