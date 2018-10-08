// testing the operations on polynomials and solutions

#include <iostream>
#include <string>
#include "complexH.h"
#include "polysys.h"
#include "polysolset.h"

using namespace std;

template <class ComplexType, class RealType>
int test ( string filename );
// reads a system from file and writes it to screen

int main ( void )
{
   cout << "Testing the solutions ..." << endl;

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
      test<complexH<double>,double>(name);
   else if(choice == '1')
      test<complexH<dd_real>,dd_real>(name);
   else if(choice == '2')
      test<complexH<qd_real>,qd_real>(name);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return 0;
}

template <class ComplexType, class RealType>
int test ( string filename )
{
   ifstream solfile(filename.c_str());

   PolySolSet<ComplexType,RealType> sols(solfile);

   cout << "The solutions on the file " << filename << " :" << endl;

   sols.print();

   return 0;
}
