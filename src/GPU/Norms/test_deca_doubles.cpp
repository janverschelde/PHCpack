// Test on deca double functions/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "deca_double_functions.h"
#include "dbl10_sqrt_kernels.h"

using namespace std;

int my_sqrt
 ( double *rtb, double *rix, double *rmi, double *rrg, double *rpk,
   double *ltb, double *lix, double *lmi, double *lrg, double *lpk,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the deca double
 *   given by ten parts respectively in rtb, rix, rmi, rrg, rpk, ltb, lix,
 *   lmi, lrg, and lpk, in as many interations as the value of max_steps.
 *   Returns in rtb, rix, rmi, rrg, rpk, ltb, lix, lmi, lrg, and lpk 
 *   the value of the square root.
 *   The integer on return is 1 if the error is too large,
 *   or 0 if Newton's method converged properly. */ 

int main ( void )
{
   const int max = 9;
   
   double twortb_h = 2.0;
   double tworix_h = 0.0;
   double twormi_h = 0.0;
   double tworrg_h = 0.0;
   double tworpk_h = 0.0;
   double twoltb_h = 0.0;
   double twolix_h = 0.0;
   double twolmi_h = 0.0;
   double twolrg_h = 0.0;
   double twolpk_h = 0.0;
   double twortb_d = 2.0;
   double tworix_d = 0.0;
   double twormi_d = 0.0;
   double tworrg_d = 0.0;
   double tworpk_d = 0.0;
   double twoltb_d = 0.0;
   double twolix_d = 0.0;
   double twolmi_d = 0.0;
   double twolrg_d = 0.0;
   double twolpk_d = 0.0;

   int fail = my_sqrt(&twortb_h,&tworix_h,&twormi_h,&tworrg_h,&tworpk_h,
                      &twoltb_h,&twolix_h,&twolmi_h,&twolrg_h,&twolpk_h,max);

   if(fail != 0)
      cout << "Test failed!" << endl;
   else
   {
      cout << "Test passed." << endl;

      GPU_dbl10_sqrt
         (&twortb_d,&tworix_d,&twormi_d,&tworrg_d,&tworpk_d,
          &twoltb_d,&twolix_d,&twolmi_d,&twolrg_d,&twolpk_d,max);

      cout << "GPU computed sqrt :" << endl;
      daf_write_doubles
        (twortb_d,tworix_d,twormi_d,tworrg_d,tworpk_d,
         twoltb_d,twolix_d,twolmi_d,twolrg_d,twolpk_d);

      double err = abs(twortb_h - twortb_d) + abs(tworix_h - tworix_d)
                 + abs(twormi_h - twormi_d) + abs(tworrg_h - tworrg_d)
                 + abs(tworpk_h - tworpk_d)
                 + abs(twoltb_h - twoltb_d) + abs(twolix_h - twolix_d)
                 + abs(twolmi_h - twolmi_d) + abs(twolrg_h - twolrg_d)
                 + abs(twolpk_h - twolpk_d);
      cout << "  error : " << err << endl; 

      if(err < 1.0e-154)
         cout << "GPU test on deca doubles passed." << endl;
      else
         cout << "GPU test on deca doubles failed!" << endl;
   }
   return 0;
}

int my_sqrt
 ( double *rtb, double *rix, double *rmi, double *rrg, double *rpk,
   double *ltb, double *lix, double *lmi, double *lrg, double *lpk,
   int max_steps )
{
   const double tol = 1.0e-154;

   double x_rtb,x_rix,x_rmi,x_rrg,x_rpk,z_rtb,z_rix,z_rmi,z_rrg,z_rpk;
   double x_ltb,x_lix,x_lmi,x_lrg,x_lpk,z_ltb,z_lix,z_lmi,z_lrg,z_lpk;
   double y_rtb,y_rix,y_rmi,y_rrg,y_rpk,e_rtb,e_rix,e_rmi,e_rrg,e_rpk;
   double y_ltb,y_lix,y_lmi,y_lrg,y_lpk,e_ltb,e_lix,e_lmi,e_lrg,e_lpk;
   double a_rtb,a_rix,a_rmi,a_rrg,a_rpk,a_ltb,a_lix,a_lmi,a_lrg,a_lpk;
   int i;

   x_rtb = 2.0; x_rix = 0.0; x_rmi = 0.0; x_rrg = 0.0; x_rpk = 0.0;
   x_ltb = 0.0; x_lix = 0.0; x_lmi = 0.0; x_lrg = 0.0; x_lpk = 0.0;

   cout << "\nRunning Newton's method for sqrt ...\n";

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
   *rtb = x_rtb; *rix = x_rix; *rmi = x_rmi, *rrg = x_rrg; *rpk = x_rpk;
   *ltb = x_ltb; *lix = x_lix; *lmi = x_lmi, *lrg = x_lrg; *lpk = x_lpk;

   return int(a_rtb + a_rix + a_rmi + a_rrg + a_rpk
            + a_ltb + a_lix + a_lmi + a_lrg + a_lpk > tol);

   return 0;
}
