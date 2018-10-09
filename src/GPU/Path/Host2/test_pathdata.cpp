// Test on the methods on the class Path, as defined in path_data.h. 

#include <iostream>
#include <string>
#include <cmath>
#include "complexH.h"
#include "path_data.h"

using namespace std;

template <class ComplexType, class RealType>
void test ( string filename );
/*
 * Opens the file with name filename.
 * The file should be the output file of phc -p.
 * The solutions on file are read into an instance
 * of the path class. */

int main ( void )
{
   cout << "Testing the Path class ..." << endl;

   string name;
   cout << "-> give a file name : "; cin >> name;

   char choice;

   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;
   cout << endl;

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
void test ( string filename )
{
   ifstream file(filename.c_str());

   int dim; file >> dim;

   cout << "The dimension read from file : " << dim << endl;

   Path<ComplexType,RealType> p(dim,file);

   p.print_phc();
}
