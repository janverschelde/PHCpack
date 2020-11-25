// Test on octo double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"

using namespace std;

int my_sqrt
 ( double *hihihi, double *lohihi, double *hilohi, double *lolohi,
   double *hihilo, double *lohilo, double *hilolo, double *lololo,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the octo double
 *   given by eight parts respectively in hihihi, lohihi, hilohi, lolohi,
 *   hihilo, lohilo, hilolo, and lololo, in as many interations
 *   as the value of max_steps.  Returns in hihihi, lohihi, hilohi, lolohi,
 *   hihilo, lohilo, hilolo, and lololo the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 9;
   double twohihihi = 2.0;
   double twolohihi = 0.0;
   double twohilohi = 0.0;
   double twololohi = 0.0;
   double twohihilo = 0.0;
   double twolohilo = 0.0;
   double twohilolo = 0.0;
   double twolololo = 0.0;

   int fail = my_sqrt(&twohihihi,&twolohihi,&twohilohi,&twololohi,
                      &twohihilo,&twolohilo,&twohilolo,&twolololo,max);

   if(fail == 0)
      cout << "Test passed." << endl;
   else
      cout << "Test failed!" << endl;

   return 0;
}

int my_sqrt
 ( double *hihihi, double *lohihi, double *hilohi, double *lolohi,
   double *hihilo, double *lohilo, double *hilolo, double *lololo,
   int max_steps )
{
   const double tol = 1.0e-120;

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
   int i;

   x_hihihi = 2.0; x_lohihi = 0.0; x_hilohi = 0.0; x_lolohi = 0.0;
   x_hihilo = 0.0; x_lohilo = 0.0; x_hilolo = 0.0; x_lololo = 0.0;

   cout << "\nRunning Newton's method for ...\n";

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
      odf_write_doubles(z_hihihi,z_lohihi,z_hilohi,z_lolohi,
                        z_hihilo,z_lohilo,z_hilolo,z_lololo);
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
      odf_write_doubles(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
                        x_hihilo,x_lohilo,x_hilolo,x_lololo);
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
   odf_sqrt(2.0,0.0,0,0.0,0.0,0.0,0.0,0.0,
            &y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
            &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo);
   cout << "sqrt(2) :" << endl;
   odf_write_doubles(y_hihihi,y_lohihi,y_hilohi,y_lolohi,
                     y_hihilo,y_lohilo,y_hilolo,y_lololo);
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

   *hihihi = x_hihihi; *lohihi = x_lohihi;
   *hilohi = x_hilohi; *lolohi = x_lolohi;
   *hihilo = x_hihilo; *lohilo = x_lohilo;
   *hilolo = x_hilolo; *lololo = x_lololo;

   return int(a_hihihi + a_lohihi + a_hilohi + a_lolohi
            + a_hihilo + a_lohilo + a_hilolo + a_lololo > tol);
}
