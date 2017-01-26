// testing the operations on evaluation and differentiation

#include <iostream>
#include <string>
#include "complexH.h"
#include "poly.h"

using namespace std;

int double_test ( string filename );
// test in double precision, reads a system from file

int double_double_test ( string filename );
// test in double double precision, reads a system from file

int quad_double_test ( string filename );
// test in quad double precision, reads a system from file

int main ( void )
{
   cout << "Testing the polynomials ..." << endl;

   string name;

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

int double_test ( string filename )
{
   PolySys< complexH<double>, double > polynomials;

   polynomials.read_file(filename);

   cout << "The polynomials on the file " << filename << " :" << endl;

   polynomials.print();

   cout << "The dimension : " << polynomials.dim << endl;
 
   complexH<double> *arg = new complexH<double>[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(1.0,0.0);

   complexH<double> *val = new complexH<double>[polynomials.dim];
   val = polynomials.eval(arg);
  
   cout << "The value of the system at a vector of ones :" << endl;
   for(int i=0; i<polynomials.dim; i++) cout << val[i];

   return 0;
}

int double_double_test ( string filename )
{
   PolySys< complexH<dd_real>, dd_real > polynomials;

   polynomials.read_file(filename);

   cout << "The polynomials on the file " << filename << " :" << endl;

   polynomials.print();

   cout << "The dimension : " << polynomials.dim << endl;
 
   complexH<dd_real> *arg = new complexH<dd_real>[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(1.0,0.0);

   complexH<dd_real> *val = new complexH<dd_real>[polynomials.dim];
   val = polynomials.eval(arg);
  
   cout << "The value of the system at a vector of ones :" << endl;
   for(int i=0; i<polynomials.dim; i++) cout << val[i];

   return 0;
}

int quad_double_test ( string filename )
{
   PolySys< complexH<qd_real>, qd_real > polynomials;

   polynomials.read_file(filename);

   cout << "The polynomials on the file " << filename << " :" << endl;

   polynomials.print();

   cout << "The dimension : " << polynomials.dim << endl;
 
   complexH<qd_real> *arg = new complexH<qd_real>[polynomials.dim];
   for(int i=0; i<polynomials.dim; i++) arg[i].init(1.0,0.0);

   complexH<qd_real> *val = new complexH<qd_real>[polynomials.dim];
   val = polynomials.eval(arg);
  
   cout << "The value of the system at a vector of ones :" << endl;
   for(int i=0; i<polynomials.dim; i++) cout << val[i];

   return 0;
}
