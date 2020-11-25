// Test on penta double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "penta_double_functions.h"
#include "dbl5_sqrt_kernels.h"

using namespace std;

int my_sqrt
 ( double *tb, double *ix, double *mi, double *rg, double *pk,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the penta double
 *   given by five parts respectively in tb, ix, mi, rg, and pk,
 *   in as many interations as the value of max_steps.
 *   Returns in tb, ix, mi, rg, and pk the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 8;

   double twotb_h = 2.0;
   double twoix_h = 0.0;
   double twomi_h = 0.0;
   double tworg_h = 0.0;
   double twopk_h = 0.0;
   double twotb_d = 2.0;
   double twoix_d = 0.0;
   double twomi_d = 0.0;
   double tworg_d = 0.0;
   double twopk_d = 0.0;

   int fail = my_sqrt(&twotb_h,&twoix_h,&twomi_h,&tworg_h,&twopk_h,max);

   if(fail != 0)
      cout << "Test failed!" << endl;
   else
   {
      cout << "Test passed." << endl;

      GPU_dbl5_sqrt(&twotb_d,&twoix_d,&twomi_d,&tworg_d,&twopk_d,max);

      cout << "GPU computed sqrt :" << endl;
      pdf_write_doubles(twotb_d,twoix_d,twomi_d,tworg_d,twopk_d);

      double err = abs(twotb_h - twotb_d)
                 + abs(twoix_h - twoix_d)
                 + abs(twomi_h - twomi_d)
                 + abs(tworg_h - tworg_d)
                 + abs(twopk_h - twopk_d);
      cout << "  error : " << err << endl; 

      if(err < 1.0e-76)
         cout << "GPU test on penta doubles passed." << endl;
      else
         cout << "GPU test on penta doubles failed!" << endl;
   }
   return 0;
}

int my_sqrt
 ( double *tb, double *ix, double *mi, double *rg, double *pk,
   int max_steps )
{
   const double tol = 1.0e-76;

   double x_tb,x_ix,x_mi,x_rg,x_pk,z_tb,z_ix,z_mi,z_rg,z_pk;
   double y_tb,y_ix,y_mi,y_rg,y_pk,e_tb,e_ix,e_mi,e_rg,e_pk;
   double a_tb,a_ix,a_mi,a_rg,a_pk;
   int i;

   x_tb = 2.0; x_ix = 0.0; x_mi = 0.0; x_rg = 0.0; x_pk = 0.0;

   cout << "\nRunning Newton's method for ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      pdf_copy(x_tb,x_ix,x_mi,x_rg,x_pk,&y_tb,&y_ix,&y_mi,&y_rg,&y_pk);
      pdf_sqr(x_tb,x_ix,x_mi,x_rg,x_pk,&z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
      cout << "x*x : " << endl;
      pdf_write_doubles(z_tb,z_ix,z_mi,z_rg,z_pk); cout << endl;
      pdf_inc(&z_tb,&z_ix,&z_mi,&z_rg,&z_pk,2.0,0.0,0.0,0.0,0.0);
      pdf_div(z_tb,z_ix,z_mi,z_rg,z_pk,x_tb,x_ix,x_mi,x_rg,x_pk,
              &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
      pdf_mul_pd_d(z_tb,z_ix,z_mi,z_rg,z_pk,0.5,
                   &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
      pdf_copy(z_tb,z_ix,z_mi,z_rg,z_pk,&x_tb,&x_ix,&x_mi,&x_rg,&x_pk);
      cout << "after step " << i << " : " << endl;
      pdf_write_doubles(x_tb,x_ix,x_mi,x_rg,x_pk);
      pdf_sub(x_tb,x_ix,x_mi,x_rg,x_pk,y_tb,y_ix,y_mi,y_rg,y_pk,
              &e_tb,&e_ix,&e_mi,&e_rg,&e_pk); 
      pdf_abs(e_tb,e_ix,e_mi,e_rg,e_pk,&a_tb,&a_ix,&a_mi,&a_rg,&a_pk);
      cout << "  error : "<< a_tb << endl;
   }

   *tb = x_tb; *ix = x_ix; *mi = x_mi, *rg = x_rg; *pk = x_pk;

   return int(a_tb + a_ix + a_mi + a_rg + a_pk > tol);
}
