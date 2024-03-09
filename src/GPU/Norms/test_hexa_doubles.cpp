// Test on octo double functions.

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "hexa_double_functions.h"
#include "dbl16_sqrt_kernels.h"

using namespace std;

int my_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the hexa double
 *   given by sixteen doubles respectively.
 *   Because this is a test on hexa double precision, all computations are
 *   done in the full hexa double precision, which is of course wasteful.
 *   Returns in sixteen doubles the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 10;
   double twohihihihi_h = 2.0;
   double twolohihihi_h = 0.0;
   double twohilohihi_h = 0.0;
   double twololohihi_h = 0.0;
   double twohihilohi_h = 0.0;
   double twolohilohi_h = 0.0;
   double twohilolohi_h = 0.0;
   double twolololohi_h = 0.0;
   double twohihihilo_h = 0.0;
   double twolohihilo_h = 0.0;
   double twohilohilo_h = 0.0;
   double twololohilo_h = 0.0;
   double twohihilolo_h = 0.0;
   double twolohilolo_h = 0.0;
   double twohilololo_h = 0.0;
   double twololololo_h = 0.0;
   double twohihihihi_d = 2.0;
   double twolohihihi_d = 0.0;
   double twohilohihi_d = 0.0;
   double twololohihi_d = 0.0;
   double twohihilohi_d = 0.0;
   double twolohilohi_d = 0.0;
   double twohilolohi_d = 0.0;
   double twolololohi_d = 0.0;
   double twohihihilo_d = 0.0;
   double twolohihilo_d = 0.0;
   double twohilohilo_d = 0.0;
   double twololohilo_d = 0.0;
   double twohihilolo_d = 0.0;
   double twolohilolo_d = 0.0;
   double twohilololo_d = 0.0;
   double twololololo_d = 0.0;

   int fail;

   fail = my_sqrt(&twohihihihi_h,&twolohihihi_h,&twohilohihi_h,&twololohihi_h,
                  &twohihilohi_h,&twolohilohi_h,&twohilolohi_h,&twolololohi_h,
                  &twohihihilo_h,&twolohihilo_h,&twohilohilo_h,&twololohilo_h,
                  &twohihilolo_h,&twolohilolo_h,&twohilololo_h,&twololololo_h,max);

   if(fail != 0)
      cout << "Test failed!" << endl;
   else
   {
      cout << "Test passed." << endl;

      GPU_dbl16_sqrt
         (&twohihihihi_d,&twolohihihi_d,&twohilohihi_d,&twololohihi_d,
          &twohihilohi_d,&twolohilohi_d,&twohilolohi_d,&twolololohi_d,
          &twohihihilo_d,&twolohihilo_d,&twohilohilo_d,&twololohilo_d,
          &twohihilolo_d,&twolohilolo_d,&twohilololo_d,&twololololo_d,max);

      cout << "GPU computed sqrt :" << endl;
      hdf_write_doubles
         (twohihihihi_d,twolohihihi_d,twohilohihi_d,twololohihi_d,
          twohihilohi_d,twolohilohi_d,twohilolohi_d,twolololohi_d,
          twohihihilo_d,twolohihilo_d,twohilohilo_d,twololohilo_d,
          twohihilolo_d,twolohilolo_d,twohilololo_d,twololololo_d);

      double err = abs(twohihihihi_h - twohihihihi_d)
                 + abs(twolohihihi_h - twolohihihi_d)
                 + abs(twohilohihi_h - twohilohihi_d)
                 + abs(twololohihi_h - twololohihi_d)
                 + abs(twohihilohi_h - twohihilohi_d)
                 + abs(twolohilohi_h - twolohilohi_d)
                 + abs(twohilolohi_h - twohilolohi_d)
                 + abs(twolololohi_h - twolololohi_d)
                 + abs(twohihihilo_h - twohihihilo_d)
                 + abs(twolohihilo_h - twolohihilo_d)
                 + abs(twohilohilo_h - twohilohilo_d)
                 + abs(twololohilo_h - twololohilo_d)
                 + abs(twohihilolo_h - twohihilolo_d)
                 + abs(twolohilolo_h - twolohilolo_d)
                 + abs(twohilololo_h - twohilololo_d)
                 + abs(twololololo_h - twololololo_d);
      cout << "  error : " << err << endl; 

      if(err < 1.0e-240)
         cout << "GPU test on hexa doubles passed." << endl;
      else
         cout << "GPU test on hexa doubles failed!" << endl;

   }
   return 0;
}

int my_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int max_steps )
{
   const double tol = 1.0e-240;

   double x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi;
   double x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi;
   double x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo;
   double x_hihilolo,x_lohilolo,x_hilololo,x_lolololo;
   double y_hihihihi,y_lohihihi,y_hilohihi,y_lolohihi;
   double y_hihilohi,y_lohilohi,y_hilolohi,y_lololohi;
   double y_hihihilo,y_lohihilo,y_hilohilo,y_lolohilo;
   double y_hihilolo,y_lohilolo,y_hilololo,y_lolololo;
   double z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi;
   double z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi;
   double z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo;
   double z_hihilolo,z_lohilolo,z_hilololo,z_lolololo;
   double e_hihihihi,e_lohihihi,e_hilohihi,e_lolohihi;
   double e_hihilohi,e_lohilohi,e_hilolohi,e_lololohi;
   double e_hihihilo,e_lohihilo,e_hilohilo,e_lolohilo;
   double e_hihilolo,e_lohilolo,e_hilololo,e_lolololo;
   double a_hihihihi,a_lohihihi,a_hilohihi,a_lolohihi;
   double a_hihilohi,a_lohilohi,a_hilolohi,a_lololohi;
   double a_hihihilo,a_lohihilo,a_hilohilo,a_lolohilo;
   double a_hihilolo,a_lohilolo,a_hilololo,a_lolololo;
   int i;

   x_hihihihi = *hihihihi; x_lohihihi = *lohihihi;
   x_hilohihi = *hilohihi; x_lolohihi = *lolohihi;
   x_hihilohi = *hihilohi; x_lohilohi = *lohilohi;
   x_hilolohi = *hilolohi; x_lololohi = *lololohi;
   x_hihihilo = *hihihilo; x_lohihilo = *lohihilo;
   x_hilohilo = *hilohilo; x_lolohilo = *lolohilo;
   x_hihilolo = *hihilolo; x_lohilolo = *lohilolo;
   x_hilololo = *hilololo; x_lolololo = *lolololo;

   cout << "\nRunning Newton's method for sqrt ...\n";

   cout << scientific << setprecision(16);

   for(i=1; i <= max_steps; i++)
   {
      cout << "step " << i << " : " << endl;
      hdf_copy(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
               x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
               x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
               x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,
               &y_hihihihi,&y_lohihihi,&y_hilohihi,&y_lolohihi,
               &y_hihilohi,&y_lohilohi,&y_hilolohi,&y_lololohi,
               &y_hihihilo,&y_lohihilo,&y_hilohilo,&y_lolohilo,
               &y_hihilolo,&y_lohilolo,&y_hilololo,&y_lolololo);
      hdf_sqr(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
              x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
              x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
              x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,
              &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
              &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
              &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
              &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);
      cout << "x*x : " << endl;
      hdf_write_doubles(z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi,
                        z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi,
                        z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo,
                        z_hihilolo,z_lohilolo,z_hilololo,z_lolololo);
      hdf_inc(&z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
              &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
              &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
              &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo,
              *hihihihi,*lohihihi,*hilohihi,*lolohihi,
              *hihilohi,*lohilohi,*hilolohi,*lololohi,
              *hihihilo,*lohihilo,*hilohilo,*lolohilo,
              *hihilolo,*lohilolo,*hilololo,*lolololo);
      hdf_div(z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi,
              z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi,
              z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo,
              z_hihilolo,z_lohilolo,z_hilololo,z_lolololo,
              x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
              x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
              x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
              x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,
              &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
              &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
              &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
              &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);
      hdf_mul_hd_d(z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi,
                   z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi,
                   z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo,
                   z_hihilolo,z_lohilolo,z_hilololo,z_lolololo,0.5,
                   &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
                   &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
                   &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
                   &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);
      hdf_copy(z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi,
               z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi,
               z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo,
               z_hihilolo,z_lohilolo,z_hilololo,z_lolololo,
               &x_hihihihi,&x_lohihihi,&x_hilohihi,&x_lolohihi,
               &x_hihilohi,&x_lohilohi,&x_hilolohi,&x_lololohi,
               &x_hihihilo,&x_lohihilo,&x_hilohilo,&x_lolohilo,
               &x_hihilolo,&x_lohilolo,&x_hilololo,&x_lolololo);
      cout << "after step " << i << " : " << endl;
      hdf_write_doubles(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
                        x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
                        x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
                        x_hihilolo,x_lohilolo,x_hilololo,x_lolololo);
      hdf_sub(x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi,
              x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi,
              x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo,
              x_hihilolo,x_lohilolo,x_hilololo,x_lolololo,
              y_hihihihi,y_lohihihi,y_hilohihi,y_lolohihi,
              y_hihilohi,y_lohilohi,y_hilolohi,y_lololohi,
              y_hihihilo,y_lohihilo,y_hilohilo,y_lolohilo,
              y_hihilolo,y_lohilolo,y_hilololo,y_lolololo,
              &e_hihihihi,&e_lohihihi,&e_hilohihi,&e_lolohihi,
              &e_hihilohi,&e_lohilohi,&e_hilolohi,&e_lololohi,
              &e_hihihilo,&e_lohihilo,&e_hilohilo,&e_lolohilo,
              &e_hihilolo,&e_lohilolo,&e_hilololo,&e_lolololo);
      hdf_abs(e_hihihihi,e_lohihihi,e_hilohihi,e_lolohihi,
              e_hihilohi,e_lohilohi,e_hilolohi,e_lololohi,
              e_hihihilo,e_lohihilo,e_hilohilo,e_lolohilo,
              e_hihilolo,e_lohilolo,e_hilololo,e_lolololo,
              &a_hihihihi,&a_lohihihi,&a_hilohihi,&a_lolohihi,
              &a_hihilohi,&a_lohilohi,&a_hilolohi,&a_lololohi,
              &a_hihihilo,&a_lohihilo,&a_hilohilo,&a_lolohilo,
              &a_hihilolo,&a_lohilolo,&a_hilololo,&a_lolololo);
      cout << "  error : "<< a_hihihihi << endl;
   }
   *hihihihi = x_hihihihi; *lohihihi = x_lohihihi;
   *hilohihi = x_hilohihi; *lolohihi = x_lolohihi;
   *hihilohi = x_hihilohi; *lohilohi = x_lohilohi;
   *hilolohi = x_hilolohi; *lololohi = x_lololohi;
   *hihihilo = x_hihihilo; *lohihilo = x_lohihilo;
   *hilohilo = x_hilohilo; *lolohilo = x_lolohilo;
   *hihilolo = x_hihilolo; *lohilolo = x_lohilolo;
   *hilololo = x_hilololo; *lolololo = x_lolololo;

   return int(a_hihihihi + a_lohihihi + a_hilohihi + a_lolohihi
            + a_hihilohi + a_lohilohi + a_hilolohi + a_lololohi
            + a_hihihilo + a_lohihilo + a_hilohilo + a_lolohilo
            + a_hihilolo + a_lohilolo + a_hilololo + a_lolololo > tol);
}
