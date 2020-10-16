// The file cmplx5_norm_host.cpp defines the code for the functions
// specified in cmplx5_norm_host.h.

#include "penta_double_functions.h"
#include "cmplx5_norm_host.h"

void make_copy
 ( int dim,
   double *orgretb, double *orgreix, double *orgremi, double *orgrerg,
   double *orgrepk,
   double *orgimtb, double *orgimix, double *orgimmi, double *orgimrg,
   double *orgimpk,
   double *dupretb, double *dupreix, double *dupremi, double *duprerg,
   double *duprepk,
   double *dupimtb, double *dupimix, double *dupimmi, double *dupimrg,
   double *dupimpk )
{
   for(int i=0; i<dim; i++)
   {
      dupretb[i] = orgretb[i];
      dupreix[i] = orgreix[i];
      dupremi[i] = orgremi[i];
      duprerg[i] = orgrerg[i];
      duprepk[i] = orgrepk[i];
      dupimtb[i] = orgimtb[i];
      dupimix[i] = orgimix[i];
      dupimmi[i] = orgimmi[i];
      dupimrg[i] = orgimrg[i];
      dupimpk[i] = orgimpk[i];
   }
}

void CPU_norm
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk )
{
   double sumtb = 0.0;
   double sumix = 0.0;
   double summi = 0.0;
   double sumrg = 0.0;
   double sumpk = 0.0;
   double prodtb,prodix,prodmi,prodrg,prodpk;

   for(int k=0; k<dim; k++) // sum = sum + vre[k]*vre[k] + vim[k]*vim[k]
   {
      pdf_sqr(  vretb[k],vreix[k],vremi[k],vrerg[k],vrepk[k],
              &prodtb, &prodix, &prodmi, &prodrg, &prodpk);
      pdf_inc(&sumtb,&sumix,&summi,&sumrg,&sumpk,
              prodtb,prodix,prodmi,prodrg,prodpk);
      pdf_sqr(  vimtb[k],vimix[k],vimmi[k],vimrg[k],vimpk[k],
              &prodtb, &prodix, &prodmi, &prodrg, &prodpk);
      pdf_inc(&sumtb,&sumix,&summi,&sumrg,&sumpk,
              prodtb,prodix,prodmi,prodrg,prodpk);
   }
   pdf_sqrt( sumtb, sumix, summi, sumrg, sumpk,
            normtb,normix,normmi,normrg,normpk);
}

void CPU_normalize
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double normtb, double normix, double normmi, double normrg, double normpk )
{
   for(int i=0; i<dim; i++)
   {
      pdf_div( vretb[i], vreix[i], vremi[i], vrerg[i], vrepk[i],
              normtb,   normix,   normmi,   normrg,   normpk,
              &vretb[i],&vreix[i],&vremi[i],&vrerg[i],&vrepk[i]);
      pdf_div( vimtb[i], vimix[i], vimmi[i], vimrg[i], vimpk[i],
              normtb,   normix,   normmi,   normrg,   normpk,
              &vimtb[i],&vimix[i],&vimmi[i],&vimrg[i],&vimpk[i]);
   }
}
