// The file random5_vectors.cpp defines the code for the functions
// specified in random5_vectors.h.

#include "random_numbers.h"
#include "random5_vectors.h"
#include "penta_double_functions.h"

void random_double5_vectors
 ( int dim, double *vtb_host, double *vix_host, double *vmi_host,
   double *vrg_host, double *vpk_host,
   double *vtb_device, double *vix_device, double *vmi_device,
   double *vrg_device, double *vpk_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      vtb_host[k] = r;   vtb_device[k] = r;
      vix_host[k] = 0.0; vix_device[k] = 0.0;
      vmi_host[k] = 0.0; vmi_device[k] = 0.0;
      vrg_host[k] = 0.0; vrg_device[k] = 0.0;
      vpk_host[k] = 0.0; vpk_device[k] = 0.0;
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
   double rnd;
   double cosrnd_tb,cosrnd_ix,cosrnd_mi,cosrnd_rg,cosrnd_pk;
   double sinrnd_tb,sinrnd_ix,sinrnd_mi,sinrnd_rg,sinrnd_pk;

   for(int k=0; k<dim; k++)
   {
      rnd = random_double(); // rnd is in [-1, +1]
      sinrnd_tb = rnd; sinrnd_ix = 0.0; sinrnd_mi = 0.0;
      sinrnd_rg = 0.0; sinrnd_pk = 0.0;

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
