// Tests the operations on evaluation and differentiation.
// When prompted for a system, cyclic5 is a good test case.

#include <iostream>
#include <string>
#include "complexH.h"
#include "eval_host.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( string filename );
// reads a double precision system from file, calls eval_test

template <class ComplexType, class RealType>
int eval_test ( PolySys<ComplexType,RealType>& polynomials, ComplexType* arg );
// evaluates and differentiates the polynomials at 0, 1, 2, ...
// once with poly methods and once with eval_host methods

template <class ComplexType, class RealType>
int basic_eval ( PolySys<ComplexType,RealType>& polynomials, ComplexType* arg );
// evaluates and differentiates the polynomials at arg with basic methods

template <class ComplexType, class RealType>
int ade_eval ( PolySys<ComplexType,RealType>& polynomials, ComplexType* arg );
// evaluates and differentiates the polynomials at arg
// with methods of algorithmic differentiation

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

   return eval_test<ComplexType,RealType>(polynomials);
}

template <class ComplexType, class RealType>
int eval_test ( PolySys<ComplexType,RealType>& polynomials )
{
   cout << "The dimension : " << polynomials.dim << endl;

   ComplexType* arg = new ComplexType[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i] = ComplexType(i,0.0);
   // arg[i].init(i,0.0);

   int fail = basic_eval<ComplexType,RealType>(polynomials,arg);

   if(fail != 0)
      fail = ade_eval<ComplexType,RealType>(polynomials,arg);

   return fail;
}

template <class ComplexType, class RealType>
int basic_eval ( PolySys<ComplexType,RealType>& polynomials, ComplexType* arg )
{
   ComplexType* val = new ComplexType[polynomials.dim]; // function value
   ComplexType** jac = new ComplexType*[polynomials.dim]; // Jacobian
   for(int i=0; i<polynomials.dim; i++)
      jac[i] = new ComplexType[polynomials.dim];

   polynomials.eval(arg,val,jac); // evaluate and differentiate

   cout << "The value of the system at 0, 1, 2, ... :" << endl;
   for(int i=0; i<polynomials.dim; i++) cout << val[i];
   cout << "The derivatives at 0, 1, 2, ... :" << endl;
   for(int i=0; i<polynomials.dim; i++)
   {
      cout << "All derivatives of polynomial " << i << " :" << endl;
      for(int j=0; j<polynomials.dim; j++)
      {
         int idx = i*polynomials.dim + j;
         cout << jac[i][j];
      }
   }
   return 0;
}

template <class ComplexType, class RealType>
int ade_eval ( PolySys<ComplexType,RealType>& polynomials, ComplexType* arg )
{
   CPUInstHom<ComplexType,RealType> ped; // data for eval_host
   Workspace<ComplexType> wrk;
   ComplexType alpha,t;

   alpha.init(0.0,0.0); // initialize the data for eval_host
   t.init(0.0,0.0);
   ped.init(polynomials,polynomials.dim,polynomials.n_eq,0,alpha);
   ped.init_workspace(wrk);

   ped.eval(wrk,arg,t); // evaluate and differentiate

   cout << "The value of the system at 0, 1, 2, ... :" << endl;
   for(int j=0; j<polynomials.dim; j++)
      cout << wrk.matrix[polynomials.dim*polynomials.dim + j];
   for(int j=0; j<polynomials.dim; j++)
   {
      cout << "All derivatives of polynomial " << j << " :" << endl;
      for(int i=0; i<polynomials.dim; i++)
      {
         int idx = i*polynomials.dim + j;
         cout << wrk.matrix[idx];
      }
   }
   return 0;
}
