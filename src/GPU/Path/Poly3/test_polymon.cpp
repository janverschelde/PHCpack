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

void random_exponents ( int dim, int expmax, int* exponents );
/*
 * Generates dim random positive integers as exponents in [0, expmax],
 * uniformly distributed.
 * Space must be allocated in exponents for dim integers. */

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> random_monomial ( int dim, int expmax );
/*
 * Returns a random monomial, with a random complex coefficient,
 * for dim variables, with exponents in [0, expmax]. */

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> prompted_monomial ( int dim );
/*
 * Prompts the user for a coefficient and dim exponents
 * to define a monomial which is returned. */

template <class ComplexType, class RealType>
ComplexType plain_eval ( PolyMon<ComplexType,RealType>& m, ComplexType *x );
/*
 * Applies the straightforward algorithm to evaluate m at x. */

template <class ComplexType, class RealType>
void test_evaluation ( PolyMon<ComplexType,RealType>& m );
/*
 * Test the application of the straightforward and Speelpenning algorithm 
 * to evaluate m at x. */

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

void random_exponents ( int dim, int expmax, int* exponents )
{
   for(int idx=0; idx<dim; idx++)
      exponents[idx] = rand() % (expmax+1);
}

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> random_monomial ( int dim, int expmax )
{
   ComplexType ran = random_complex<ComplexType,RealType>();
   RealType coefficient[2];
   coefficient[0] = ran.real;
   coefficient[1] = ran.imag;

   int exponents[dim];
   random_exponents(dim,expmax,exponents);

   PolyMon<ComplexType,RealType> result
      = PolyMon<ComplexType,RealType>(dim, exponents, coefficient);

   return result;
}

template <class ComplexType, class RealType>
PolyMon<ComplexType,RealType> prompted_monomial ( int dim )
{
   RealType coefficient[2];

   cout << "-> give the real part of the coefficient : ";
   cin >> coefficient[0];
   cout << "-> give the imaginary part of the coefficient : ";
   cin >> coefficient[1];

   int exponents[dim];

   for(int idx=0; idx<dim; idx++)
   {
      cout << "-> give the exponent for variable " << idx+1 << " : ";
      cin >> exponents[idx];
   }

   PolyMon<ComplexType,RealType> result
      = PolyMon<ComplexType,RealType>(dim, exponents, coefficient);

   return result;
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
   ComplexType val3 = m.speel(point, derivatives);
   cout << "The value at a random point : " << val3 << endl;
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

      PolyMon<ComplexType,RealType> monomial
         = random_monomial<ComplexType,RealType>(dim, 1);

      test_evaluation<ComplexType,RealType>(monomial);
   }
   return 0;
}
