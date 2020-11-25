// Test on penta double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "penta_double_functions.h"

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

   double twotb = 2.0;
   double twoix = 0.0;
   double twomi = 0.0;
   double tworg = 0.0;
   double twopk = 0.0;

   int fail = my_sqrt(&twotb,&twoix,&twomi,&tworg,&twopk,max);

   if(fail == 0)
      cout << "Test passed." << endl;
   else
      cout << "Test failed!" << endl;

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
   pdf_sqrt(2.0,0.0,0.0,0.0,0.0,&y_tb,&y_ix,&y_mi,&y_rg,&y_pk);
   cout << "sqrt(2) :" << endl;
   pdf_write_doubles(y_tb,y_ix,y_mi,y_rg,y_pk);
   pdf_sub(x_tb,x_ix,x_mi,x_rg,x_pk,y_tb,y_ix,y_mi,y_rg,y_pk,
           &e_tb,&e_ix,&e_mi,&e_rg,&e_pk);
   pdf_abs(e_tb,e_ix,e_mi,e_rg,e_pk,&a_tb,&a_ix,&a_mi,&a_rg,&a_pk);
   cout << "  error : "<< a_tb << endl;

   *tb = x_tb; *ix = x_ix; *mi = x_mi, *rg = x_rg; *pk = x_pk;

   return int(a_tb + a_ix + a_mi + a_rg + a_pk > tol);
}
