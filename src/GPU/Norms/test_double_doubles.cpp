// Test on double double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double.h"
#include "double_double_functions.h"

using namespace std;

int my_sqrt ( void );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root. */

int main ( void )
{
   int fail = my_sqrt();

   return 0;
}

int my_sqrt ( void )
{
   double n[2],x[2],y[2],z[2],e[2],a[2];
   double x_hi,x_lo,z_hi,z_lo,y_hi,y_lo,e_hi,e_lo;
   const int max_steps = 7;
   const char sqrt2[] = "1.4142135623730950488016887242097\0";
   int i;

   x_hi = 2.0; x_lo = 0.0; n[0] = 2.0; n[1] = 0.0; // dd_copy(n,x);

   cout << "\nrunning Newton's method for sqrt(2) ...\n";
   dd_read(sqrt2,y);

   cout << scientific << setprecision(16);
   cout << "y hi : " << y[0] << endl;
   cout << "y lo : " << y[1] << endl;
   ddf_sqrt(x_hi,x_lo,&y_hi,&y_lo);
   cout << "y_hi : " << y_hi << endl;
   cout << "y_lo : " << y_lo << endl;
   a[0] = y_hi; a[1] = y_lo;
   cout << "    sqrt2 : " << sqrt2 << endl;
   cout << "ddf_sqrt2 : "; dd_write(a,32); cout << endl;

   x[0] = x_hi; x[1] = x_lo;
   cout << "step 0 : "; dd_write(x,32); cout << endl;
   for(i=1; i <= max_steps; i++)
   {
      ddf_sqr(x_hi,x_lo,&z_hi,&z_lo);           // z = x*x
      ddf_inc(&z_hi,&z_lo,n[0],n[1]);           // z += n
      ddf_div(z_hi,z_lo,x_hi,x_lo,&z_hi,&z_lo); // z = z/x
      z[0] = z_hi; z[1] = z_lo;
      ddf_mlt_d(&z_hi,&z_lo,0.5);               // z *= 0.5
      cout << "step " << i << " : "; dd_write(z,32);
      x_hi = z_hi; x_lo = z_lo;                 // dd_copy(z,x);
      ddf_sub(x_hi,x_lo,y[0],y[1],&e_hi,&e_lo); // dd_sub(x,y,e);
      ddf_abs(e_hi,e_lo,&a[0],&a[1]);
      cout << "  error : "; dd_write(a,3); cout << endl;
   }

   return 0;
}
