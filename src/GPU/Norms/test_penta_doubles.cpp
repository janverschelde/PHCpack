// Test on penta double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "penta_double_functions.h"

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
   double x_tb,x_ix,x_mi,x_rg,x_pk,z_tb,z_ix,z_mi,z_rg,z_pk;
   double y_tb,y_ix,y_mi,y_rg,y_pk,e_tb,e_ix,e_mi,e_rg,e_pk;
   double a_tb,a_ix,a_mi,a_rg,a_pk;
   const int max_steps = 8;
   int i;

   x_tb = 2.0; x_ix = 0.0; x_mi = 0.0; x_rg = 0.0; x_pk = 0.0;

   cout << "\nrunning Newton's method for sqrt(2) ...\n";

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

   return 0;
}
