// Test on the methods on the class Workspace, 
// as defined in workspace_host.h. 

#include <iostream>
#include <string>
#include <cmath>
#include "complexH.h"
#include "workspace_host.h"

using namespace std;

template <class ComplexType>
void test ( void );
/*
 * Runs a test. */

int main ( void )
{
   cout << "Testing the class Workspace ..." << endl;

   char choice;

   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;
   cout << endl;

   if(choice == '0')
      test< complexH<double> >();
   else if(choice == '1')
      test< complexH<dd_real> >();
   else if(choice == '2')
      test< complexH<qd_real> >();
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return 0;
}

template <class ComplexType>
void test ( void )
{
   Workspace<ComplexType> wrk;
}
