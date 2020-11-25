// Test on double double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double.h"
#include "double_double_functions.h"

using namespace std;

int my_sqrt ( double *hi, double *lo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root
 *   of the double double with high part in hi and low part in lo,
 *   in as many iterations as the value of max_steps.
 *   Returns in hi and lo the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 7;
   double twohi = 2.0;
   double twolo = 0.0;

   int fail = my_sqrt(&twohi,&twolo,max);

   if(fail == 0)
      cout << "Test passed." << endl;
   else
      cout << "Test failed!" << endl;

   return 0;
}

int my_sqrt ( double *hi, double *lo, int max_steps )
{
   const double tol = 1.0e-28;

   double n[2],x[2],y[2],z[2],e[2],a[2];
   double x_hi,x_lo,z_hi,z_lo,y_hi,y_lo,e_hi,e_lo;

   const char sqrt2[] = "1.4142135623730950488016887242097\0";
   int i;

   x_hi = *hi; x_lo = *lo; n[0] = x_hi; n[1] = x_lo; // dd_copy(n,x);

   cout << "\nRunning Newton's method for sqrt ...\n";
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
   *hi = x_hi; *lo = x_lo;

   return int(a[0] + a[1] > tol);
}
