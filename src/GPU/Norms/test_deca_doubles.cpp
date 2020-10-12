// Test on deca double functions/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "deca_double_functions.h"

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
   double x_rtb,x_rix,x_rmi,x_rrg,x_rpk,z_rtb,z_rix,z_rmi,z_rrg,z_rpk;
   double x_ltb,x_lix,x_lmi,x_lrg,x_lpk,z_ltb,z_lix,z_lmi,z_lrg,z_lpk;
   double y_rtb,y_rix,y_rmi,y_rrg,y_rpk,e_rtb,e_rix,e_rmi,e_rrg,e_rpk;
   double y_ltb,y_lix,y_lmi,y_lrg,y_lpk,e_ltb,e_lix,e_lmi,e_lrg,e_lpk;
   double a_rtb,a_rix,a_rmi,a_rrg,a_rpk,a_ltb,a_lix,a_lmi,a_lrg,a_lpk;
   const int max_steps = 9;
   int i;

   x_rtb = 2.0; x_rix = 0.0; x_rmi = 0.0; x_rrg = 0.0; x_rpk = 0.0;
   x_ltb = 0.0; x_lix = 0.0; x_lmi = 0.0; x_lrg = 0.0; x_lpk = 0.0;

   cout << "\nrunning Newton's method for sqrt(2) ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      daf_copy(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,x_ltb,x_lix,x_lmi,x_lrg,x_lpk,
               &y_rtb,&y_rix,&y_rmi,&y_rrg,&y_rpk,
               &y_ltb,&y_lix,&y_lmi,&y_lrg,&y_lpk);
      daf_sqr(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,x_ltb,x_lix,x_lmi,x_lrg,x_lpk,
              &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
              &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);
      cout << "x*x : " << endl;
      daf_write_doubles(z_rtb,z_rix,z_rmi,z_rrg,z_rpk,
                        z_ltb,z_lix,z_lmi,z_lrg,z_lpk);
      daf_inc(&z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
              &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk,
              2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
      daf_div(z_rtb,z_rix,z_rmi,z_rrg,z_rpk,z_ltb,z_lix,z_lmi,z_lrg,z_lpk,
              x_rtb,x_rix,x_rmi,x_rrg,x_rpk,x_ltb,x_lix,x_lmi,x_lrg,x_lpk,
              &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
              &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);
      daf_mul_da_d(z_rtb,z_rix,z_rmi,z_rrg,z_rpk,
                   z_ltb,z_lix,z_lmi,z_lrg,z_lpk,0.5,
                   &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
                   &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);
      daf_copy(z_rtb,z_rix,z_rmi,z_rrg,z_rpk,z_ltb,z_lix,z_lmi,z_lrg,z_lpk,
               &x_rtb,&x_rix,&x_rmi,&x_rrg,&x_rpk,
               &x_ltb,&x_lix,&x_lmi,&x_lrg,&x_lpk);
      cout << "after step " << i << " : " << endl;
      daf_write_doubles(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
                        x_ltb,x_lix,x_lmi,x_lrg,x_lpk);
      daf_sub(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,x_ltb,x_lix,x_lmi,x_lrg,x_lpk,
              y_rtb,y_rix,y_rmi,y_rrg,y_rpk,y_ltb,y_lix,y_lmi,y_lrg,y_lpk,
              &e_rtb,&e_rix,&e_rmi,&e_rrg,&e_rpk,
              &e_ltb,&e_lix,&e_lmi,&e_lrg,&e_lpk); 
      daf_abs(e_rtb,e_rix,e_rmi,e_rrg,e_rpk,e_ltb,e_lix,e_lmi,e_lrg,e_lpk,
              &a_rtb,&a_rix,&a_rmi,&a_rrg,&a_rpk,
              &a_ltb,&a_lix,&a_lmi,&a_lrg,&a_lpk);
      cout << "  error : "<< a_rtb << endl;
   }
   cout << "sqrt(2) : " << endl;
   daf_sqrt(2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
            &y_rtb,&y_rix,&y_rmi,&y_rrg,&y_rpk,
            &y_ltb,&y_lix,&y_lmi,&y_lrg,&y_lpk);
   daf_write_doubles(y_rtb,y_rix,y_rmi,y_rrg,y_rpk,
                     y_ltb,y_lix,y_lmi,y_lrg,y_lpk);
   daf_sub(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,x_ltb,x_lix,x_lmi,x_lrg,x_lpk,
           y_rtb,y_rix,y_rmi,y_rrg,y_rpk,y_ltb,y_lix,y_lmi,y_lrg,y_lpk,
           &e_rtb,&e_rix,&e_rmi,&e_rrg,&e_rpk,
           &e_ltb,&e_lix,&e_lmi,&e_lrg,&e_lpk); 
   daf_abs(e_rtb,e_rix,e_rmi,e_rrg,e_rpk,e_ltb,e_lix,e_lmi,e_lrg,e_lpk,
           &a_rtb,&a_rix,&a_rmi,&a_rrg,&a_rpk,
           &a_ltb,&a_lix,&a_lmi,&a_lrg,&a_lpk);
   cout << "  error : "<< a_rtb << endl;

   return 0;
}
