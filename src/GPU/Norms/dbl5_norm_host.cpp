// The file dbl5_norm_host.cpp defines the code for the functions
// specified in dbl5_norm_host.h.

#include "dbl5_norm_host.h"
#include "penta_double_functions.h"

void make_copy
 ( int dim,
   double *orgtb, double *orgix, double *orgmi, double *orgrg, double *orgpk,
   double *duptb, double *dupix, double *dupmi, double *duprg, double *duppk )
{
   for(int i=0; i<dim; i++)
   {
      duptb[i] = orgtb[i];
      dupix[i] = orgix[i];
      dupmi[i] = orgmi[i];
      duprg[i] = orgrg[i];
      duppk[i] = orgpk[i];
   }
}

void CPU_norm
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk,
   int dim, double *normtb, double *normix, double *normmi,
   double *normrg, double *normpk )
{
   double sumtb = 0.0;
   double sumix = 0.0;
   double summi = 0.0;
   double sumrg = 0.0;
   double sumpk = 0.0;
   double prodtb,prodix,prodmi,prodrg,prodpk;

   for(int i=0; i<dim; i++) // sum = sum + v[i]*v[i]
   {
      pdf_sqr(vtb[i],vix[i],vmi[i],vrg[i],vpk[i],
              &prodtb,&prodix,&prodmi,&prodrg,&prodpk);
      pdf_inc(&sumtb,&sumix,&summi,&sumrg,&sumpk,
              prodtb,prodix,prodmi,prodrg,prodpk);
   }
   pdf_sqrt(sumtb,sumix,summi,sumrg,sumpk,
            normtb,normix,normmi,normrg,normpk);
}

void CPU_normalize
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   double normtb, double normix, double normmi, double normrg, double normpk )
{
   for(int i=0; i<dim; i++)
      pdf_div(vtb[i],vix[i],vmi[i],vrg[i],vpk[i],
              normtb,normix,normmi,normrg,normpk,
              &vtb[i],&vix[i],&vmi[i],&vrg[i],&vpk[i]);
}
