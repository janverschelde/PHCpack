// testing the operations on polynomials

#include <iostream>
#include <string>
#include "poly.h"

using namespace std;

int main ( void )
{
   cout << "Testing the polynomials ..." << endl;

   string filename = "../../../Demo/cyclic5";
   PolySys polynomials;

   cout << "Reading from file " << filename << endl;

   polynomials.read_file(filename);

   cout << "The polynomials on the file :" << endl;

   polynomials.print();

   return 0;
}
