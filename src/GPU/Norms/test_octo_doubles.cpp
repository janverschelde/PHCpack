// Test on octo double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"

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
   double x_hihihi,x_lohihi,x_hilohi,x_lolohi;
   double x_hihilo,x_lohilo,x_hilolo,x_lololo;
   double y_hihihi,y_lohihi,y_hilohi,y_lolohi;
   double y_hihilo,y_lohilo,y_hilolo,y_lololo;
   double z_hihihi,z_lohihi,z_hilohi,z_lolohi;
   double z_hihilo,z_lohilo,z_hilolo,z_lololo;
   double e_hihihi,e_lohihi,e_hilohi,e_lolohi;
   double e_hihilo,e_lohilo,e_hilolo,e_lololo;
   double a_hihihi,a_lohihi,a_hilohi,a_lolohi;
   double a_hihilo,a_lohilo,a_hilolo,a_lololo;
   const int max_steps = 9;
   int i;

   x_hihihi = 2.0; x_lohihi = 0.0; x_hilohi = 0.0; x_lolohi = 0.0;
   x_hihilo = 0.0; x_lohilo = 0.0; x_hilolo = 0.0; x_lololo = 0.0;

   cout << "\nrunning Newton's method for sqrt(2) ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      odf_copy(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
               x_hihilo,x_lohilo,x_hilolo,x_lololo,
               &y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
               &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo);
      odf_sqr(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
              x_hihilo,x_lohilo,x_hilolo,x_lololo,
              &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
              &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);
      cout << "x*x : " << endl;
      cout << "  hihihi : " << z_hihihi;
      cout << "  lohihi : " << z_lohihi << endl;
      cout << "  hilohi : " << z_hilohi;
      cout << "  lolohi : " << z_lolohi << endl;
      cout << "  hihilo : " << z_hihilo;
      cout << "  lohilo : " << z_lohilo << endl;
      cout << "  hilolo : " << z_hilolo;
      cout << "  lololo : " << z_lololo << endl;
      odf_inc(&z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
              &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo,
              2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
      odf_div(z_hihihi,z_lohihi,z_hilohi,z_lolohi,
              z_hihilo,z_lohilo,z_hilolo,z_lololo,
              x_hihihi,x_lohihi,x_hilohi,x_lolohi,
              x_hihilo,x_lohilo,x_hilolo,x_lololo,
              &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
              &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);
      odf_mul_od_d(z_hihihi,z_lohihi,z_hilohi,z_lolohi,
                   z_hihilo,z_lohilo,z_hilolo,z_lololo,0.5,
                   &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
                   &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);
      odf_copy(z_hihihi,z_lohihi,z_hilohi,z_lolohi,
               z_hihilo,z_lohilo,z_hilolo,z_lololo,
               &x_hihihi,&x_lohihi,&x_hilohi,&x_lolohi,
               &x_hihilo,&x_lohilo,&x_hilolo,&x_lololo);
      cout << "after step " << i << " : " << endl;
      cout << "  hihihi : " << x_hihihi;
      cout << "  lohihi : " << x_lohihi << endl;
      cout << "  hilohi : " << x_hilohi;
      cout << "  lolohi : " << x_lolohi << endl;
      cout << "  hihilo : " << x_hihilo;
      cout << "  lohilo : " << x_lohilo << endl;
      cout << "  hilolo : " << x_hilolo;
      cout << "  lololo : " << x_lololo << endl;
      odf_sub(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
              x_hihilo,x_lohilo,x_hilolo,x_lololo,
              y_hihihi,y_lohihi,y_hilohi,y_lolohi,
              y_hihilo,y_lohilo,y_hilolo,y_lololo,
              &e_hihihi,&e_lohihi,&e_hilohi,&e_lolohi,
              &e_hihilo,&e_lohilo,&e_hilolo,&e_lololo);
      odf_abs(e_hihihi,e_lohihi,e_hilohi,e_lolohi,
              e_hihilo,e_lohilo,e_hilolo,e_lololo,
              &a_hihihi,&a_lohihi,&a_hilohi,&a_lolohi,
              &a_hihilo,&a_lohilo,&a_hilolo,&a_lololo);
      cout << "  error : "<< a_hihihi << endl;
   }
   return 0;
}
