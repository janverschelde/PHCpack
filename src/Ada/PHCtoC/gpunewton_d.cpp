// C++ function called by Ada routine to accelerate Newton's method

#include <iostream>
#include "syscon.h"
#include "solcon.h"

using namespace std;

extern "C" int gpunewton_d ( void )
/*
 * DESCRIPTION :
 *   A C++ function to accelerate Newton's method,
 *   encapsulated as a C function for to be called from Ada. */
{
   int fail,dim,len;

   cout << endl;
   cout << "Acceleration of Newton's method ..." << endl;

   fail = syscon_number_of_polynomials(&dim);
   cout << "number of polynomials : " << dim << endl;
   fail = solcon_number_of_solutions(&len);
   cout << "number of solutions : " << len << endl;

   return 0;
}
