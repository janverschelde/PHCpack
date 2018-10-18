// Tests the start system by evaluation of all its solutions.

#include <iostream>
#include <string>
#include "complexH.h"
#include "polysys.h"
#include "polysolset.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( string startfile, int precision );
/*
 * Reads start system from file, calls start_test.
 * The precision should be 16, 32, or 64 for printing solutions.
 */

template <class ComplexType, class RealType>
int start_test
 ( PolySys<ComplexType,RealType>& startpols,
   PolySolSet<ComplexType,RealType>& sols,
   int precision );
/*
 * Evaluates all start solutions in sols in startpols.
 * The precision should be 16 (double), 32 (double double),
 * or 64 (quad double), for use in the printing of the solutions.
 */

int main ( void )
{
   cout << "Testing the start system ..." << endl;

   string name4start;
   cout << "\n-> give a file name for the start system : ";
   cin >> name4start;

   char choice;

   cout << endl << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   int fail = 0;

   if(choice == '0')
      fail = test<complexH<double>,double>(name4start,16);
   else if(choice == '1')
      fail = test<complexH<dd_real>,dd_real>(name4start,32);
   else if(choice == '2')
      fail = test<complexH<qd_real>,qd_real>(name4start,64);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return fail;
}

template <class ComplexType, class RealType>
int test ( string startfile, int precision )
{
   PolySys<ComplexType,RealType> startsys;

   cout << endl << "Reading from file " << startfile << " ..." << endl;
   startsys.read_file(startfile);
   cout << "The polynomials on the file " << startfile << " :" << endl;
   startsys.print();

   ifstream solfile(startfile.c_str());
   PolySolSet<ComplexType,RealType> sols(solfile);

   cout << endl << "-> read " << sols.n_sol
        << " solutions from " << startfile << endl;

   // cout << "The solutions on the file " << startfile << " :" << endl;
   // sols.print();

   const int neq = startsys.n_eq;
   const int dim = startsys.dim;

   ComplexType* funval = new ComplexType[neq];
   ComplexType** jacmat = new ComplexType*[neq];
   for(int idx=0; idx<neq; idx++) jacmat[idx] = new ComplexType[dim];

   ComplexType* point=sols.get_sol(0);

   startsys.eval(point, funval, jacmat);

   ComplexType sum(0.0,0.0);

   for(int cnt=0; cnt<sols.n_sol; cnt++)
      for(int idx=0; idx<dim; idx++) sum = sum + funval[idx];

   cout << "The sum of all function values : " << sum << endl;

   return 0;
}
