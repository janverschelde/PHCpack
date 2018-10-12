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
void plain_diff
 ( PolySys<ComplexType,RealType> &polsys, ComplexType *x,
   ComplexType** derivatives );
/*
 * Applies the straightforward algorithm to differentiate polsys at x.
 * All evaluated partial derivatives of the polynomial with index idx
 * are returned in derivatives[idx].
 */

template <class ComplexType, class RealType>
void test_evaluation ( PolySys<ComplexType,RealType> &polsys );
/*
 * Tests the evaluation and differentation of the polynomial system
 * at a random point. */

template <class ComplexType, class RealType>
int random_test ( int dim );
/*
 * Runs an interactive test on a randomly generated system with
 * as many variables as the value of the dimension dim. */

template <class ComplexType, class RealType>
int polsys_test ( void );
/*
 * Prompts the user for a file name for a system and then
 * evaluates and differentiates the system at a random point. */

int main ( void )
{
   srand(time(NULL));

   cout << endl << "Testing the polynomial systems ..." << endl;

   char choice;

   cout << endl << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   cout << endl << "Generate and test random system ? (y/n) ";
   char rantest; cin >> rantest;

   if(rantest == 'y')
   {
      int dimension;

      cout << endl << "Give the dimension : "; cin >> dimension;

      if(choice == '0')
         random_test<complexH<double>,double>(dimension);
      else if(choice == '1')
         random_test<complexH<dd_real>,dd_real>(dimension);
      else if(choice == '2')
         random_test<complexH<qd_real>,qd_real>(dimension);
      else
         cout << "Invalid choice " << choice << " for the precision." << endl; 
   }
   else
   {
      if(choice == '0')
         polsys_test<complexH<double>,double>();
      else if(choice == '1')
         polsys_test<complexH<dd_real>,dd_real>();
      else if(choice == '2')
         polsys_test<complexH<qd_real>,qd_real>();
      else
         cout << "Invalid choice " << choice << " for the precision." << endl; 
   }

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

   if(expmax > 1)              // not all monomials are products of variables
      result.eval_base = true; // there are factors common to all derivatives

   return result;
}

template <class ComplexType, class RealType>
void write_polynomials ( PolySys<ComplexType,RealType>& p )
{
   cout << "The dimension : " << p.dim << endl;

   for(int idx=0; idx<p.n_eq; idx++)
   {
      write_polynomial<ComplexType,RealType>(*(p.eq[idx]));
      cout << endl;
   }
}

template <class ComplexType, class RealType>
ComplexType* plain_eval
 ( PolySys<ComplexType,RealType>& polsys, ComplexType *x )
{
   const int neq = polsys.n_eq;

   ComplexType* result = new ComplexType[neq];

   for(int idx=0; idx<neq; idx++)
      result[idx] = plain_eval<ComplexType,RealType>(*polsys.eq[idx],x);

   return result;
}

template <class ComplexType, class RealType>
void plain_diff
 ( PolySys<ComplexType,RealType> &polsys, ComplexType *x,
   ComplexType** derivatives )
{
   for(int idx=0; idx<polsys.n_eq; idx++)
      plain_diff<ComplexType, RealType>(*polsys.eq[idx],x,derivatives[idx]);
}

template <class ComplexType, class RealType>
void test_evaluation ( PolySys<ComplexType,RealType>& polsys )
{
   const int neq = polsys.n_eq;
   const int dim = polsys.dim;

   cout << "Testing the evaluation ..." << endl;
   ComplexType* point = new ComplexType[dim];
   random_point<ComplexType,RealType>(dim, point);

   ComplexType** derivatives = new ComplexType*[dim];
   ComplexType** plainderivs = new ComplexType*[dim];

   for(int idx=0; idx<neq; idx++)
   {
      derivatives[idx] = new ComplexType[dim];
      plainderivs[idx] = new ComplexType[dim];
   }
   ComplexType* val1 = polsys.eval(point);
   ComplexType* val2 = plain_eval<ComplexType,RealType>(polsys,point);
   plain_diff<ComplexType,RealType>(polsys,point,plainderivs);

   ComplexType* val3;
   if(polsys.eval_base == false)
      val3 = polsys.eval(point,derivatives);
   else
   {  
      val3 = new ComplexType[dim];
      polsys.update_max_deg_base();

      ComplexType** deg_table = polsys.allocate_deg_table();
      polsys.compute_deg_table(point,deg_table);
      polsys.eval(point,val3,derivatives,deg_table);
   }

   for(int idx=0; idx<neq; idx++)
   {
      cout << "at " << idx << " : " << val1[idx];
      cout << "at " << idx << " : " << val2[idx];
      cout << "at " << idx << " : " << val3[idx];
   }

   cout << "derivatives computed plainly and by Speelpenning :" << endl;
   for(int polidx=0; polidx<neq; polidx++)
      for(int varidx=0; varidx<dim; varidx++)
      {
         cout << "at " << polidx << "," << varidx << " : "
              << plainderivs[polidx][varidx];
         cout << "at " << polidx << "," << varidx << " : "
              << derivatives[polidx][varidx];
      }
}

template <class ComplexType, class RealType>
int random_test ( int dim )
{
   int nbterms,expmax;

   cout << "Give the number of terms : "; cin >> nbterms;
   cout << "Give the largest exponent : "; cin >> expmax;

   PolySys<ComplexType,RealType> polynomials
     = random_polynomials<ComplexType,RealType>(dim,nbterms,expmax);

   cout << endl << "The terms in a random polynomial system :" << endl;
   write_polynomials<ComplexType,RealType>(polynomials);

   test_evaluation<ComplexType,RealType>(polynomials);

   return 0;
}

template <class ComplexType, class RealType>
int polsys_test ( void )
{
   cout << "-> give a file name : ";
   string name; cin >> name;

   PolySys<ComplexType,RealType> polynomials;
   polynomials.read_file(name,1);
   cout << "The polynomials on the file " << name << " :" << endl;
   polynomials.print();

   if(polynomials.eval_base)
      cout << "There is a common factor in the polynomials." << endl;
   else
      cout << "All monomials are products of variables." << endl;

   cout << endl << "The terms in a random polynomial system :" << endl;
   write_polynomials<ComplexType,RealType>(polynomials);

   test_evaluation<ComplexType,RealType>(polynomials);

   return 0;
}
