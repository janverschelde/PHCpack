// Test on quad double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"

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
   double x_hihi,x_lohi,x_hilo,x_lolo,y_hihi,y_lohi,y_hilo,y_lolo;
   double z_hihi,z_lohi,z_hilo,z_lolo;
   double e_hihi,e_lohi,e_hilo,e_lolo,a_hihi,a_lohi,a_hilo,a_lolo;
   const int max_steps = 8;
   int i;

   x_hihi = 2.0; x_lohi = 0.0; x_hilo = 0.0; x_lolo = 0.0;

   cout << "\nrunning Newton's method for sqrt(2) ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      qdf_copy(x_hihi,x_lohi,x_hilo,x_lolo,&y_hihi,&y_lohi,&y_hilo,&y_lolo);
      qdf_sqr(x_hihi,x_lohi,x_hilo,x_lolo,&z_hihi,&z_lohi,&z_hilo,&z_lolo);
      cout << "x*x : " << endl;
      cout << "  hihi : " << z_hihi;
      cout << "  lohi : " << z_lohi << endl;
      cout << "  hilo : " << z_hilo;
      cout << "  lolo : " << z_lolo << endl;
      qdf_inc(&z_hihi,&z_lohi,&z_hilo,&z_lolo,2.0,0.0,0.0,0.0);
      qdf_div(z_hihi,z_lohi,z_hilo,z_lolo,x_hihi,x_lohi,x_hilo,x_lolo,
              &z_hihi,&z_lohi,&z_hilo,&z_lolo);
      qdf_mul_qd_d(z_hihi,z_lohi,z_hilo,z_lolo,0.5,
                   &z_hihi,&z_lohi,&z_hilo,&z_lolo);
      qdf_copy(z_hihi,z_lohi,z_hilo,z_lolo,
               &x_hihi,&x_lohi,&x_hilo,&x_lolo);
      cout << "after step " << i << " : " << endl;
      cout << "  hihi : " << x_hihi;
      cout << "  lohi : " << x_lohi << endl;
      cout << "  hilo : " << x_hilo;
      cout << "  lolo : " << x_lolo << endl;
      qdf_sub(x_hihi,x_lohi,x_hilo,x_lolo,y_hihi,y_lohi,y_hilo,y_lolo,
              &e_hihi,&e_lohi,&e_hilo,&e_lolo);
      qdf_abs(e_hihi,e_lohi,e_hilo,e_lolo,&a_hihi,&a_lohi,&a_hilo,&a_lolo);
      cout << "  error : "<< a_hihi << endl;
   }
   return 0;
}
