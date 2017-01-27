// Tests the operations on evaluation and differentiation.
// When prompted for a system, cyclic5 is a good test case.

#include <iostream>
#include <string>
#include "complexH.h"
#include "eval_host.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( PolySys<ComplexType,RealType> polynomials );
// evaluates and differentiates the polynomials at 0, 1, 2, ...
// once with poly methods and once with eval_host methods

int double_test ( string filename );
// reads a double precision system from file, calls test

int double_double_test ( string filename );
// reads a double double precision system from file, calls test

int quad_double_test ( string filename );
// reads a quad double precision system from file, calls test

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

   if(choice == '0')
      double_test(name);
   else if(choice == '1')
      double_double_test(name);
   else if(choice == '2')
      quad_double_test(name);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return 0;
}

template <class ComplexType, class RealType>
int test ( PolySys<ComplexType,RealType> polynomials )
{
   cout << "The dimension : " << polynomials.dim << endl;

   ComplexType* arg = new ComplexType[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(i,0.0);
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

int double_test ( string filename )
{
   PolySys<complexH<double>,double> polynomials;

   polynomials.read_file(filename);
   cout << "The polynomials on the file " << filename << " :" << endl;
   polynomials.print();

   return test<complexH<double>,double>(polynomials);
}

int double_double_test ( string filename )
{
   PolySys< complexH<dd_real>, dd_real > polynomials;

   polynomials.read_file(filename);
   cout << "The polynomials on the file " << filename << " :" << endl;
   polynomials.print();

   return test<complexH<dd_real>,dd_real>(polynomials);
}

int quad_double_test ( string filename )
{
   PolySys< complexH<qd_real>, qd_real > polynomials;

   polynomials.read_file(filename);
   cout << "The polynomials on the file " << filename << " :" << endl;
   polynomials.print();

   return test<complexH<qd_real>,qd_real>(polynomials);
}
