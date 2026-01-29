/* Tests the collection of functions to split doubles. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "splitting_doubles.h"

int test_half_split ( void );
/*
 * Generates a random number, splits in two equal sized halves, and then
 * checks if adding the parts gives the original number. */

int test_quarter_sums
 ( double x, double x0, double x1, double x2, double x3 );
/*
 * Tests whether the sums of the quarters x0, x1, x2, and x3
 * add up to be equal to x.  Returns 0 if okay, 1 if not. */

int test_quarter_split ( void );
/*
 * Generates a random number, splits in four equal sized halves, and then
 * checks if adding the parts gives the original number. */

int test_quarter_balancer ( void );
/*
 * Generates a random number, splits in four equal sized halves, and then
 * checks if the quarters are balanced. */

using namespace std;

int main ( void )
{
   int fail;

   fail = test_half_split();
   if(fail == 1)
      cout << "\nTest on half split failed?!!!\n\n";
   else
      cout << "\nTest on half split succeeded.\n\n";

   fail = test_quarter_split();
   if(fail == 1)
      cout << "\nTest on quarter split failed?!!!\n\n";
   else
      cout << "\nTest on quarter split succeeded.\n\n";

   fail = test_quarter_balancer();
   if(fail == 1)
      cout << "\nTest on quarter balancer failed?!!!\n\n";
   else
      cout << "\nTest on quarter balancer succeeded.\n\n";

   return fail;
}

int test_half_split ( void )
{
   double x,x0,x1,s,e;

   srand(time(NULL));

   cout << scientific << setprecision(16);

   x = ((double) rand())/RAND_MAX;
   cout << " x : " << x << endl;

   half_split(x, &x0, &x1, 1);

   cout << "x0 : " << x0 << endl;
   cout << "x1 : " << x1 << endl;
   s = x0 + x1;
   cout << "x0 + x1 : " << s << endl;
   cout << "      x : " << x << endl;
   e = abs(x - s);

   cout << scientific << setprecision(3);

   cout << "  error : " << e << endl;

   return not(e == 0);
}

int test_quarter_sums
 ( double x, double x0, double x1, double x2, double x3 )
{
   double s,e;

   cout << scientific << setprecision(16);

   cout << "x0 : " << x0 << endl;
   cout << " b : "; write_52double(x0);
   cout << "x1 : " << x1 << endl;
   cout << " b : "; write_52double(x1);
   cout << "x2 : " << x2 << endl;
   cout << " b : "; write_52double(x2);
   cout << "x3 : " << x3 << endl;
   cout << " b : "; write_52double(x3);

   cout << "                x : " << x << endl;
   e = fabs(x - x0);
   cout << "               x0 : " << x0;
   cout << ", error : " << scientific << setprecision(3) << e << endl;
   s = x0 + x1; e = fabs(x - s);
   cout << scientific << setprecision(16);
   cout << "          x0 + x1 : " << s;
   cout << ", error : " << scientific << setprecision(3) << e << endl;
   s = s + x2; e = fabs(x - s);
   cout << scientific << setprecision(16);
   cout << "     x0 + x1 + x2 : " << s;
   cout << ", error : " << scientific << setprecision(3) << e << endl;
   s = s + x3; e = fabs(x - s);
   cout << scientific << setprecision(16);
   cout << "x0 + x1 + x2 + x3 : " << s;
   cout << ", error : " << scientific << setprecision(3) << e << endl;

   return not(e == 0.0);
}

int test_quarter_split ( void )
{
   double x,x0,x1,x2,x3;

   cout << "Give seed (0 for default) : ";
   int seed; cin >> seed;

   if(seed == 0) seed = time(NULL);
   srand(seed);

   x = ((double) rand())/RAND_MAX;
   cout << scientific << setprecision(16);
   cout << " x : " << x << endl;

   quarter_split(x, &x0, &x1, &x2, &x3, 1);

   int fail = test_quarter_sums(x, x0, x1, x2, x3);

   cout << endl << "seed used : " << seed << endl;

   return fail;
}

int test_quarter_balancer ( void )
{
   double x,x0,x1,x2,x3;

   cout << "Give seed (0 for default) : ";
   int seed; cin >> seed;

   if(seed == 0) seed = time(NULL);
   srand(seed);

   x = ((double) rand())/RAND_MAX;
   cout << scientific << setprecision(16);
   cout << " x : " << x << endl;

   double factor;
   make_exponent_zero(&x, &factor, 1);

   quarter_split(x, &x0, &x1, &x2, &x3, 0);

   bool b01 = is_quarter_balanced(x0, x1, 1);
   bool b12 = is_quarter_balanced(x1, x2, 1);
   bool b23 = is_quarter_balanced(x2, x3, 1);

   if(b01 and b12 and b23)
      cout << "The quarters are balanced." << endl;
   else
      cout << "The quarters are not balanced." << endl;

   if(not b01) 
   {
      quarter_balance(&x0, &x1, 1); 
      b12 = is_quarter_balanced(x1, x2, 1); // must recompute b12
   }
   if(not b12)
   {
      quarter_balance(&x1, &x2, 1); 
      b23 = is_quarter_balanced(x2, x3, 1); // must recompute b23
   }
   if(not b23) quarter_balance(&x2, &x3, 1); 

   int fail = test_quarter_sums(x, x0, x1, x2, x3);

   cout << endl << "seed used : " << seed << endl;

   return fail;
}
