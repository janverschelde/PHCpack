// Tests Newton's method with algorithmic differentiation.

#include <iostream>
#include <string>
#include "complexH.h"
#include "eval_host.h"
#include "newton_host.h"

using namespace std;

template <class ComplexType, class RealType>
int newton_test ( PolySys<ComplexType,RealType>& polynomials );
// calls Newton's method on the polynomials

template <class ComplexType, class RealType>
int test ( string filename );
// reads a double precision system from file, calls newton_test

int main ( void )
{
   cout << "Testing the polynomials ..." << endl;

   string name; // = "/Users/jan/PHCv2/Demo/cyclic5";
   cout << "-> give a file name : "; cin >> name;

   char choice;

   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;

   cout << "\nReading from file " << name << " ..." << endl;

   int fail = 0;

   if(choice == '0')
      fail = test<complexH<double>,double>(name);
   else if(choice == '1')
      fail = test<complexH<dd_real>,dd_real>(name);
   else if(choice == '2')
      fail = test<complexH<qd_real>,qd_real>(name);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return fail;
}

template <class ComplexType, class RealType>
int test ( string filename )
{
   PolySys<ComplexType,RealType> polynomials;

   polynomials.read_file(filename);
   cout << "The polynomials on the file " << filename << " :" << endl;
   polynomials.print();

   return newton_test<ComplexType,RealType>(polynomials);
}

template <class ComplexType, class RealType>
int newton_test ( PolySys<ComplexType,RealType>& polynomials )
{
   cout << "The dimension : " << polynomials.dim << endl;

   ComplexType* arg = new ComplexType[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(i,0.0);

   CPUInstHom<ComplexType,RealType> ped; // data for eval_host
   Workspace<ComplexType> wrk;
   ComplexType alpha,t;
   alpha.init(0.0,0.0); // initialize the data for eval_host
   t.init(0.0,0.0);
   ped.init(polynomials,polynomials.dim,polynomials.n_eq,0,alpha);
   ped.init_workspace(wrk);
   Parameter pars;

   double tSec_Eval = 0;
   double tSec_MGS = 0;

   bool success_cpu
      = CPU_Newton<ComplexType,RealType>(wrk,ped,pars,tSec_Eval,tSec_MGS);

   return 0;
}
