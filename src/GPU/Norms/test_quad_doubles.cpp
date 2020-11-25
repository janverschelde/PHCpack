// Test on quad double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"

using namespace std;

int my_sqrt
 ( double *hihi, double *lohi, double *hilo, double *lolo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the quad double
 *   given by four parts respectively in hihi, lohi, hilo, and lolo,
 *   in as many interations as the value of max_steps.
 *   Returns in hihi, lohi, hilo, and lolo the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 8;
   double twohihi = 2.0;
   double twolohi = 0.0;
   double twohilo = 0.0;
   double twololo = 0.0;

   int fail = my_sqrt(&twohihi,&twolohi,&twohilo,&twololo,max);

   if(fail == 0)
      cout << "Test passed." << endl;
   else
      cout << "Test failed!" << endl;

   return 0;
}

int my_sqrt
 ( double *hihi, double *lohi, double *hilo, double *lolo, int max_steps )
{
   const double tol = 1.0e-60;

   double x_hihi,x_lohi,x_hilo,x_lolo,y_hihi,y_lohi,y_hilo,y_lolo;
   double z_hihi,z_lohi,z_hilo,z_lolo;
   double e_hihi,e_lohi,e_hilo,e_lolo,a_hihi,a_lohi,a_hilo,a_lolo;
   int i;

   x_hihi = 2.0; x_lohi = 0.0; x_hilo = 0.0; x_lolo = 0.0;

   cout << "\nrunning Newton's method for sqrt ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      qdf_copy(x_hihi,x_lohi,x_hilo,x_lolo,&y_hihi,&y_lohi,&y_hilo,&y_lolo);
      qdf_sqr(x_hihi,x_lohi,x_hilo,x_lolo,&z_hihi,&z_lohi,&z_hilo,&z_lolo);
      cout << "x*x : " << endl;
      qdf_write_doubles(z_hihi,z_lohi,z_hilo,z_lolo);
      qdf_inc(&z_hihi,&z_lohi,&z_hilo,&z_lolo,2.0,0.0,0.0,0.0);
      qdf_div(z_hihi,z_lohi,z_hilo,z_lolo,x_hihi,x_lohi,x_hilo,x_lolo,
              &z_hihi,&z_lohi,&z_hilo,&z_lolo);
      qdf_mul_qd_d(z_hihi,z_lohi,z_hilo,z_lolo,0.5,
                   &z_hihi,&z_lohi,&z_hilo,&z_lolo);
      qdf_copy(z_hihi,z_lohi,z_hilo,z_lolo,
               &x_hihi,&x_lohi,&x_hilo,&x_lolo);
      cout << "after step " << i << " : " << endl;
      qdf_write_doubles(x_hihi,x_lohi,x_hilo,x_lolo);
      qdf_sub(x_hihi,x_lohi,x_hilo,x_lolo,y_hihi,y_lohi,y_hilo,y_lolo,
              &e_hihi,&e_lohi,&e_hilo,&e_lolo);
      qdf_abs(e_hihi,e_lohi,e_hilo,e_lolo,&a_hihi,&a_lohi,&a_hilo,&a_lolo);
      cout << "  error : "<< a_hihi << endl;
   }
   qdf_sqrt(2.0,0.0,0,0.0,&y_hihi,&y_lohi,&y_hilo,&y_lolo);
   cout << "sqrt(2) :" << endl;
   qdf_write_doubles(y_hihi,y_lohi,y_hilo,y_lolo);
   qdf_sub(x_hihi,x_lohi,x_hilo,x_lolo,y_hihi,y_lohi,y_hilo,y_lolo,
           &e_hihi,&e_lohi,&e_hilo,&e_lolo);
   qdf_abs(e_hihi,e_lohi,e_hilo,e_lolo,&a_hihi,&a_lohi,&a_hilo,&a_lolo);
   cout << "  error : "<< a_hihi << endl;

   *hihi = x_hihi; *lohi = x_lohi; *hilo = x_hilo; *lolo = x_lolo;

   return int(a_hihi + a_lohi + a_hilo + a_lolo > tol);
}
