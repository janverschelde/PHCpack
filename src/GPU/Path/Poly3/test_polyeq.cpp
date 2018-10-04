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
PolyEq<ComplexType,RealType> random_polynomial
 ( int dim, int nbterms, int expmax );
/*
 * Returns a random polynomial in a space with dim variables
 * and with as many terms as the value of nbterms.
 * The largest exponent is defined by the value of expmax.  */

template <class ComplexType, class RealType>
void write_polynomial ( PolyEq<ComplexType,RealType>& p );
/*
 * Writes the terms in p in tableau style format. */

template <class ComplexType, class RealType>
ComplexType plain_eval
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x );
/*
 * Applies the straightforward algorithm to evaluate p at x. */

template <class ComplexType, class RealType>
void plain_diff
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x, ComplexType *deri );
/*
 * Applies the straightforward algorithm to differentiate p at x.
 * The function returns in deri the values of all derivatives. */

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
      // print_data<ComplexType,RealType>(*p.mon[idx]);

      cout << p.mon[idx]->coef.real << "  "
           << p.mon[idx]->coef.imag;

      int varidx = 0;
      for(int posidx=0; posidx<p.mon[idx]->n_var; posidx++)
      {
         while(varidx < p.mon[idx]->pos[posidx])
         {
            cout << " 0";
            varidx = varidx+1;
         }
         cout << " " << p.mon[idx]->exp[posidx];
         varidx = varidx+1;
      }
      while(varidx < p.mon[idx]->dim)
      {
         cout << " 0";
         varidx = varidx + 1;
      }
      cout << endl;
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

   cout << endl << "Testing the evaluation ..." << endl;
   ComplexType* point = new ComplexType[dim];
   random_point<ComplexType,RealType>(dim, point);

   ComplexType* derivatives = new ComplexType[dim];
   ComplexType* plainderivs = new ComplexType[dim];

   ComplexType val1 = polynomial.eval(point);
   cout << "the value at a random point : " << val1;

   ComplexType val2 = plain_eval<ComplexType,RealType>(polynomial,point);
   cout << "the value at a random point : " << val2;

   ComplexType val3 = polynomial.eval(point,derivatives);
   cout << "the value at a random point : " << val3;

   plain_diff<ComplexType,RealType>(polynomial,point,plainderivs);

   cout << "The derivatives computed plainly and with Speelpenning :" << endl;
   for(int idx=0; idx<polynomial.dim; idx++)
   {
      cout << " at " << idx << " : " << plainderivs[idx];
      cout << " at " << idx << " : " << derivatives[idx];
   }

   return 0;
}

template <class ComplexType, class RealType>
PolyEq<ComplexType,RealType> random_polynomial
 ( int dim, int nbterms, int expmax )
{
   PolyEq<ComplexType,RealType> polynomial(dim);
   polynomial.constant = random_complex<ComplexType,RealType>();

   polynomial.n_mon = nbterms-1;

   for(int idx=0; idx<nbterms-1; idx++)
   {
      PolyMon<ComplexType,RealType>* term
         = new PolyMon<ComplexType,RealType>
                  (random_monomial<ComplexType,RealType>(dim,expmax));
      polynomial.mon.push_back(term);
   }
   return polynomial;
}

template <class ComplexType, class RealType>
ComplexType plain_eval
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x )
{
   ComplexType result = p.constant;

   for(int idx=0; idx<p.n_mon; idx++)
      result += plain_eval(*p.mon[idx],x);

   return result;
}

template <class ComplexType, class RealType>
void plain_diff
 ( PolyEq<ComplexType,RealType>& p, ComplexType* x, ComplexType *deri )
{
   ComplexType monderi[p.dim];

   for(int idx=0; idx<p.dim; idx++) deri[idx].init(0.0,0.0);

   for(int idx=0; idx<p.n_mon; idx++)
   {
      PolyMon<ComplexType,RealType> *m = p.mon[idx];

      plain_diff<ComplexType,RealType>(*m, x, monderi);

      for(int jdx=0; jdx<m->n_var; jdx++)
         deri[m->pos[jdx]] += monderi[jdx];
   }
}
