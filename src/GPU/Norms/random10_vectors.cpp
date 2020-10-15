// The file random10_vectors.cpp defines the code for the functions
// specified in random10_vectors.h.

#include "random_numbers.h"
#include "random10_vectors.h"
#include "deca_double_functions.h"

void random_deca_double
 ( double *x_rtb, double *x_rix, double *x_rmi, double *x_rrg, double *x_rpk,
   double *x_ltb, double *x_lix, double *x_lmi, double *x_lrg, double *x_lpk )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double eps3 = eps*eps2;
   const double eps4 = eps*eps3;
   const double eps5 = eps*eps4;
   const double eps6 = eps*eps5;
   const double eps7 = eps*eps6;
   const double eps8 = eps*eps7;
   const double eps9 = eps*eps8;
   const double r0 = random_double();
   const double r1 = random_double();
   const double r2 = random_double();
   const double r3 = random_double();
   const double r4 = random_double();
   const double r5 = random_double();
   const double r6 = random_double();
   const double r7 = random_double();
   const double r8 = random_double();
   const double r9 = random_double();

   *x_rtb = r0; *x_rix = 0.0; *x_rmi = 0.0; *x_rrg = 0.0; *x_rpk = 0.0;
   *x_ltb = 0.0; *x_lix = 0.0; *x_lmi = 0.0; *x_lrg = 0.0; *x_lpk = 0.0;

   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r1*eps);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r2*eps2);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r3*eps3);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r4*eps4);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r5*eps5);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r6*eps6);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r7*eps7);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r8*eps8);
   daf_inc_d(x_rtb,x_rix,x_rmi,x_rrg,x_rpk,
             x_ltb,x_lix,x_lmi,x_lrg,x_lpk,r9*eps9);
}

void random_double10_vectors
 ( int dim, double *vrtb_host, double *vrix_host, double *vrmi_host,
   double *vrrg_host, double *vrpk_host, double *vltb_host,
   double *vlix_host, double *vlmi_host, double *vlrg_host,
   double *vlpk_host, double *vrtb_device, double *vrix_device,
   double *vrmi_device, double *vrrg_device, double *vrpk_device,
   double *vltb_device, double *vlix_device, double *vlmi_device,
   double *vlrg_device, double *vlpk_device )
{
   double r_rtb,r_rix,r_rmi,r_rrg,r_rpk,r_ltb,r_lix,r_lmi,r_lrg,r_lpk;

   for(int k=0; k<dim; k++)
   {
      random_deca_double(&r_rtb,&r_rix,&r_rmi,&r_rrg,&r_rpk,
                         &r_ltb,&r_lix,&r_lmi,&r_lrg,&r_lpk);
      vrtb_host[k] = r_rtb; vrtb_device[k] = r_rtb;
      vrix_host[k] = r_rix; vrix_device[k] = r_rix;
      vrmi_host[k] = r_rmi; vrmi_device[k] = r_rmi;
      vrrg_host[k] = r_rrg; vrrg_device[k] = r_rrg;
      vrpk_host[k] = r_rpk; vrpk_device[k] = r_rpk;
      vltb_host[k] = r_ltb; vltb_device[k] = r_ltb;
      vlix_host[k] = r_lix; vlix_device[k] = r_lix;
      vlmi_host[k] = r_lmi; vlmi_device[k] = r_lmi;
      vlrg_host[k] = r_lrg; vlrg_device[k] = r_lrg;
      vlpk_host[k] = r_lpk; vlpk_device[k] = r_lpk;
   }
}

void random_complex10_vectors
 ( int dim,
   double *vrertb_host, double *vrerix_host, double *vrermi_host,
   double *vrerrg_host, double *vrerpk_host,
   double *vreltb_host, double *vrelix_host, double *vrelmi_host,
   double *vrelrg_host, double *vrelpk_host,
   double *vimrtb_host, double *vimrix_host, double *vimrmi_host,
   double *vimrrg_host, double *vimrpk_host,
   double *vimltb_host, double *vimlix_host, double *vimlmi_host,
   double *vimlrg_host, double *vimlpk_host,
   double *vrertb_device, double *vrerix_device, double *vrermi_device,
   double *vrerrg_device, double *vrerpk_device,
   double *vreltb_device, double *vrelix_device, double *vrelmi_device,
   double *vrelrg_device, double *vrelpk_device,
   double *vimrtb_device, double *vimrix_device, double *vimrmi_device,
   double *vimrrg_device, double *vimrpk_device,
   double *vimltb_device, double *vimlix_device, double *vimlmi_device,
   double *vimlrg_device, double *vimlpk_device )
{
   double rnd_rtb,rnd_rix,rnd_rmi,rnd_rrg,rnd_rpk;
   double rnd_ltb,rnd_lix,rnd_lmi,rnd_lrg,rnd_lpk;
   double cosrnd_rtb,cosrnd_rix,cosrnd_rmi,cosrnd_rrg,cosrnd_rpk;
   double cosrnd_ltb,cosrnd_lix,cosrnd_lmi,cosrnd_lrg,cosrnd_lpk;
   double sinrnd_rtb,sinrnd_rix,sinrnd_rmi,sinrnd_rrg,sinrnd_rpk;
   double sinrnd_ltb,sinrnd_lix,sinrnd_lmi,sinrnd_lrg,sinrnd_lpk;

   for(int k=0; k<dim; k++)
   {
      random_deca_double(&rnd_rtb,&rnd_rix,&rnd_rmi,&rnd_rix,&rnd_rpk,
                         &rnd_ltb,&rnd_lix,&rnd_lmi,&rnd_lix,&rnd_lpk);
      sinrnd_rtb = rnd_rtb; sinrnd_rix = rnd_rix; sinrnd_rmi = rnd_rix;
      sinrnd_rrg = rnd_rrg; sinrnd_rpk = rnd_rpk; 
      sinrnd_ltb = rnd_ltb; sinrnd_lix = rnd_lix; sinrnd_lmi = rnd_lix;
      sinrnd_lrg = rnd_lrg; sinrnd_lpk = rnd_lpk; 

      double y_rtb,y_rix,y_rmi,y_rrg,y_rpk;  // work around to compute cos
      double y_ltb,y_lix,y_lmi,y_lrg,y_lpk;  

      daf_sqr(sinrnd_rtb,sinrnd_rix,sinrnd_rmi,sinrnd_rrg,sinrnd_rpk,
              sinrnd_ltb,sinrnd_lix,sinrnd_lmi,sinrnd_lrg,sinrnd_lpk,
              &y_rtb,&y_rix,&y_rmi,&y_rrg,&y_rpk,
              &y_ltb,&y_lix,&y_lmi,&y_lrg,&y_lpk);
      daf_minus(&y_rtb,&y_rix,&y_rmi,&y_rrg,&y_rpk,
                &y_ltb,&y_lix,&y_lmi,&y_lrg,&y_lpk);       // y = -sin^2
      daf_inc_d(&y_rtb,&y_rix,&y_rmi,&y_rrg,&y_rpk,
                &y_ltb,&y_lix,&y_lmi,&y_lrg,&y_lpk,1.0);   // y = 1 - sin^2
      daf_sqrt(y_rtb,y_rix,y_rmi,y_rrg,y_rpk,y_ltb,y_lix,y_lmi,y_lrg,y_lpk,
               &cosrnd_rtb,&cosrnd_rix,&cosrnd_rmi,&cosrnd_rrg,&cosrnd_rpk,
               &cosrnd_ltb,&cosrnd_lix,&cosrnd_lmi,&cosrnd_lrg,&cosrnd_lpk);
      // cos is sqrt(1-sin^2)

      vrertb_host[k] = cosrnd_rtb; vrertb_device[k] = cosrnd_rtb;
      vrerix_host[k] = cosrnd_rix; vrerix_device[k] = cosrnd_rix;
      vrermi_host[k] = cosrnd_rmi; vrermi_device[k] = cosrnd_rmi;
      vrerrg_host[k] = cosrnd_rrg; vrerrg_device[k] = cosrnd_rrg;
      vrerpk_host[k] = cosrnd_rpk; vrerpk_device[k] = cosrnd_rpk;
      vreltb_host[k] = cosrnd_ltb; vreltb_device[k] = cosrnd_ltb;
      vrelix_host[k] = cosrnd_lix; vrelix_device[k] = cosrnd_lix;
      vrelmi_host[k] = cosrnd_lmi; vrelmi_device[k] = cosrnd_lmi;
      vrelrg_host[k] = cosrnd_lrg; vrelrg_device[k] = cosrnd_lrg;
      vrelpk_host[k] = cosrnd_lpk; vrelpk_device[k] = cosrnd_lpk;

      vimrtb_host[k] = sinrnd_rtb; vimrtb_device[k] = sinrnd_rtb;
      vimrix_host[k] = sinrnd_rix; vimrix_device[k] = sinrnd_rix;
      vimrmi_host[k] = sinrnd_rmi; vimrmi_device[k] = sinrnd_rmi;
      vimrrg_host[k] = sinrnd_rrg; vimrrg_device[k] = sinrnd_rrg;
      vimrpk_host[k] = sinrnd_rpk; vimrpk_device[k] = sinrnd_rpk;
      vimltb_host[k] = sinrnd_ltb; vimltb_device[k] = sinrnd_ltb;
      vimlix_host[k] = sinrnd_lix; vimlix_device[k] = sinrnd_lix;
      vimlmi_host[k] = sinrnd_lmi; vimlmi_device[k] = sinrnd_lmi;
      vimlrg_host[k] = sinrnd_lrg; vimlrg_device[k] = sinrnd_lrg;
      vimlpk_host[k] = sinrnd_lpk; vimlpk_device[k] = sinrnd_lpk;
   }
}
