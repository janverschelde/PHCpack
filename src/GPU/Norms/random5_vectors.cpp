// The file random5_vectors.cpp defines the code for the functions
// specified in random5_vectors.h.

#include "random_numbers.h"
#include "random5_vectors.h"
#include "penta_double_functions.h"

void random_penta_double
 ( double *x_tb, double *x_ix, double *x_mi, double *x_rg, double *x_pk )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double eps3 = eps*eps2;
   const double eps4 = eps*eps3;
   const double r0 = random_double();
   const double r1 = random_double();
   const double r2 = random_double();
   const double r3 = random_double();
   const double r4 = random_double();

   *x_tb = r0; *x_ix = 0.0; *x_mi = 0.0; *x_rg = 0.0; *x_pk = 0.0;

   pdf_inc_d(x_tb,x_ix,x_mi,x_rg,x_pk,r1*eps);
   pdf_inc_d(x_tb,x_ix,x_mi,x_rg,x_pk,r2*eps2);
   pdf_inc_d(x_tb,x_ix,x_mi,x_rg,x_pk,r3*eps3);
   pdf_inc_d(x_tb,x_ix,x_mi,x_rg,x_pk,r4*eps4);
}

void random_double5_vectors
 ( int dim, double *vtb_host, double *vix_host, double *vmi_host,
   double *vrg_host, double *vpk_host,
   double *vtb_device, double *vix_device, double *vmi_device,
   double *vrg_device, double *vpk_device )
{
   double r_tb,r_ix,r_mi,r_rg,r_pk;

   for(int k=0; k<dim; k++)
   {
      random_penta_double(&r_tb,&r_ix,&r_mi,&r_rg,&r_pk);
      vtb_host[k] = r_tb; vtb_device[k] = r_tb;
      vix_host[k] = r_ix; vix_device[k] = r_ix;
      vmi_host[k] = r_mi; vmi_device[k] = r_mi;
      vrg_host[k] = r_rg; vrg_device[k] = r_rg;
      vpk_host[k] = r_pk; vpk_device[k] = r_pk;
   }
}

void random_complex5_vectors
 ( int dim, double *vretb_host, double *vreix_host, double *vremi_host,
            double *vrerg_host, double *vrepk_host,
            double *vimtb_host, double *vimix_host, double *vimmi_host,
            double *vimrg_host, double *vimpk_host,
   double *vretb_device, double *vreix_device, double *vremi_device,
   double *vrerg_device, double *vrepk_device,
   double *vimtb_device, double *vimix_device, double *vimmi_device,
   double *vimrg_device, double *vimpk_device )
{
   double rnd_tb,rnd_ix,rnd_mi,rnd_rg,rnd_pk;
   double cosrnd_tb,cosrnd_ix,cosrnd_mi,cosrnd_rg,cosrnd_pk;
   double sinrnd_tb,sinrnd_ix,sinrnd_mi,sinrnd_rg,sinrnd_pk;

   for(int k=0; k<dim; k++)
   {
      random_penta_double(&rnd_tb,&rnd_ix,&rnd_mi,&rnd_rg,&rnd_pk);
      sinrnd_tb = rnd_tb; sinrnd_ix = rnd_ix; sinrnd_mi = rnd_mi;
      sinrnd_rg = rnd_rg; sinrnd_pk = rnd_pk; 

      double y_tb,y_ix,y_mi,y_rg,y_pk;      // work around to compute cos

      pdf_sqr(sinrnd_tb,sinrnd_ix,sinrnd_mi,sinrnd_rg,sinrnd_pk,
              &y_tb,&y_ix,&y_mi,&y_rg,&y_pk);
      pdf_minus(&y_tb,&y_ix,&y_mi,&y_rg,&y_pk);       // y = -sin^2
      pdf_inc_d(&y_tb,&y_ix,&y_mi,&y_rg,&y_pk,1.0);   // y = 1 - sin^2
      pdf_sqrt(y_tb,y_ix,y_mi,y_rg,y_pk,
               &cosrnd_tb,&cosrnd_ix,&cosrnd_mi,&cosrnd_rg,&cosrnd_pk);
      // cos is sqrt(1-sin^2)

      vretb_host[k] = cosrnd_tb; vretb_device[k] = cosrnd_tb;
      vreix_host[k] = cosrnd_ix; vreix_device[k] = cosrnd_ix;
      vremi_host[k] = cosrnd_mi; vremi_device[k] = cosrnd_mi;
      vrerg_host[k] = cosrnd_rg; vrerg_device[k] = cosrnd_rg;
      vrepk_host[k] = cosrnd_pk; vrepk_device[k] = cosrnd_pk;
      vimtb_host[k] = sinrnd_tb; vimtb_device[k] = sinrnd_tb;
      vimix_host[k] = sinrnd_ix; vimix_device[k] = sinrnd_ix;
      vimmi_host[k] = sinrnd_mi; vimmi_device[k] = sinrnd_mi;
      vimrg_host[k] = sinrnd_rg; vimrg_device[k] = sinrnd_rg;
      vimpk_host[k] = sinrnd_pk; vimpk_device[k] = sinrnd_pk;
   }
}
