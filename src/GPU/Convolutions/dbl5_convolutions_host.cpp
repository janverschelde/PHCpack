/* The file dbl5_convolutions_host.cpp defines functions
 * specified in dbl5_convolutions_host.h. */

#include "dbl5_convolutions_host.h"
#include "penta_double_functions.h"

void CPU_dbl5_product
 ( int deg, double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
            double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
            double *ztb, double *zix, double *zmi, double *zrg, double *zpk )
{
   int idx;
   double ptb,pix,pmi,prg,ppk;

   pdf_mul(xtb[0],xix[0],xmi[0],xrg[0],xpk[0],
           ytb[0],yix[0],ymi[0],yrg[0],ypk[0],
           &ztb[0],&zix[0],&zmi[0],&zrg[0],&zpk[0]);    // z[0] = x[0]*y[0]

   for(int k=1; k<=deg; k++)
   {
      pdf_mul(xtb[0],xix[0],xmi[0],xrg[0],xpk[0],
              ytb[k],yix[k],ymi[k],yrg[k],ypk[k],
              &ztb[k],&zix[k],&zmi[k],&zrg[k],&zpk[k]); // z[k] = x[0]*y[k]
      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         pdf_mul(xtb[i],xix[i],xmi[i],xrg[i],xpk[i],
                 ytb[idx],yix[idx],ymi[idx],yrg[idx],ypk[idx],
                 &ptb,&pix,&pmi,&prg,&ppk);             // x[i]*y[k-i]
         // z[k] = z[k] + x[i]*y[k-i]
         pdf_inc(&ztb[k],&zix[k],&zmi[k],&zrg[k],&zpk[k],
                 ptb,pix,pmi,prg,ppk); 
      }
   }
}

void CPU_cmplx5_product
 ( int deg,
   double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk )
{
   double rpatb,rpaix,rpami,rparg,rpapk;    // accumulates real parts
   double ipatb,ipaix,ipami,iparg,ipapk;    // accumulates imaginary parts
   double tmptb,tmpix,tmpmi,tmprg,tmppk;    // temporary penta double
   double xr0tb,xr0ix,xr0mi,xr0rg,xr0pk;    // real values in xr 
   double xi0tb,xi0ix,xi0mi,xi0rg,xi0pk;    // imaginary values in xi
   double yr0tb,yr0ix,yr0mi,yr0rg,yr0pk;    // real values in yr
   double yi0tb,yi0ix,yi0mi,yi0rg,yi0pk;    // imaginary values in yi
   int idx;

   xr0tb = xretb[0]; xr0ix = xreix[0]; xr0mi = xremi[0];
   xr0rg = xrerg[0]; xr0pk = xrepk[0];
   xi0tb = ximtb[0]; xi0ix = ximix[0]; xi0mi = ximmi[0];
   xi0rg = ximrg[0]; xi0pk = ximpk[0];
   yr0tb = yretb[0]; yr0ix = yreix[0]; yr0mi = yremi[0];
   yr0rg = yrerg[0]; yr0pk = yrepk[0];
   yi0tb = yimtb[0]; yi0ix = yimix[0]; yi0mi = yimmi[0];
   yi0rg = yimrg[0]; yi0pk = yimpk[0];
   // zre[0] = xr0*yr0 - xi0*yi0;
   pdf_mul(xr0tb,xr0ix,xr0mi,xr0rg,xr0pk,yr0tb,yr0ix,yr0mi,yr0rg,yr0pk,
           &zretb[0],&zreix[0],&zremi[0],&zrerg[0],&zrepk[0]);
   pdf_mul(xi0tb,xi0ix,xi0mi,xi0rg,xi0pk,yi0tb,yi0ix,yi0mi,yi0rg,yi0pk,
           &rpatb,&rpaix,&rpami,&rparg,&rpapk);
   pdf_minus(&rpatb,&rpaix,&rpami,&rparg,&rpapk);
   pdf_inc(&zretb[0],&zreix[0],&zremi[0],&zrerg[0],&zrepk[0],
           rpatb,rpaix,rpami,rparg,rpapk);
   // zim[0] = xi0*yr0 + xr0*yi0;
   pdf_mul(xi0tb,xi0ix,xi0mi,xi0rg,xi0pk,yr0tb,yr0ix,yr0mi,yr0rg,yr0pk,
           &zimtb[0],&zimix[0],&zimmi[0],&zimrg[0],&zimpk[0]);
   pdf_mul(xr0tb,xr0ix,xr0mi,xr0rg,xr0pk,yi0tb,yi0ix,yi0mi,yi0rg,yi0pk,
           &ipatb,&ipaix,&ipami,&iparg,&ipapk);
   pdf_inc(&zimtb[0],&zimix[0],&zimmi[0],&zimrg[0],&zimpk[0],
           ipatb,ipaix,ipami,iparg,ipapk);

   for(int k=1; k<=deg; k++)
   {
      xr0tb = xretb[0]; xr0ix = xreix[0]; xr0mi = xremi[0];
      xr0rg = xrerg[0]; xr0pk = xrepk[0];
      xi0tb = ximtb[0]; xi0ix = ximix[0]; xi0mi = ximmi[0];
      xi0rg = ximrg[0]; xi0pk = ximpk[0];
      yr0tb = yretb[k]; yr0ix = yreix[k]; yr0mi = yremi[k];
      yr0rg = yrerg[k]; yr0pk = yrepk[k];
      yi0tb = yimtb[k]; yi0ix = yimix[k]; yi0mi = yimmi[k];
      yi0rg = yimrg[k]; yi0pk = yimpk[k];
      // rpa = xr0*yr0 - xi0*yi0;
      pdf_mul(xr0tb,xr0ix,xr0mi,xr0rg,xr0pk,yr0tb,yr0ix,yr0mi,yr0rg,yr0pk,
              &rpatb,&rpaix,&rpami,&rparg,&rpapk);
      pdf_mul(xi0tb,xi0ix,xi0mi,xi0rg,xi0pk,yi0tb,yi0ix,yi0mi,yi0rg,yi0pk,
              &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_minus(&tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_inc(&rpatb,&rpaix,&rpami,&rparg,&rpapk,
              tmptb,tmpix,tmpmi,tmprg,tmppk);
      // ipa = xi0*yr0 + xr0*yi0;
      pdf_mul(xi0tb,xi0ix,xi0mi,xi0rg,xi0pk,yr0tb,yr0ix,yr0mi,yr0rg,yr0pk,
              &ipatb,&ipaix,&ipami,&iparg,&ipapk);
      pdf_mul(xr0tb,xr0ix,xr0mi,xr0rg,xr0pk,yi0tb,yi0ix,yi0mi,yi0rg,yi0pk,
              &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_inc(&ipatb,&ipaix,&ipami,&iparg,&ipapk,
              tmptb,tmpix,tmpmi,tmprg,tmppk);

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0tb = xretb[i]; xr0ix = xreix[i]; xr0mi = xremi[i];
         xr0rg = xrerg[i]; xr0pk = xrepk[i];
         xi0tb = ximtb[i]; xi0ix = ximix[i]; xi0mi = ximmi[i];
         xi0rg = ximrg[i]; xi0pk = ximpk[i];
         yr0tb = yretb[idx]; yr0ix = yreix[idx]; yr0mi = yremi[idx];
         yr0rg = yrerg[idx]; yr0pk = yrepk[idx];
         yi0tb = yimtb[idx]; yi0ix = yimix[idx]; yi0mi = yimmi[idx];
         yi0rg = yimrg[idx]; yi0pk = yimpk[idx];
         // rpa = rpa + xr0*yr0 - xi0*yi0;
         pdf_mul(xr0tb,xr0ix,xr0mi,xr0rg,xr0pk,yr0tb,yr0ix,yr0mi,yr0rg,yr0pk,
                 &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
         pdf_inc(&rpatb,&rpaix,&rpami,&rparg,&rpapk,
                 tmptb,tmpix,tmpmi,tmprg,tmppk);
         pdf_mul(xi0tb,xi0ix,xi0mi,xi0rg,xi0pk,yi0tb,yi0ix,yi0mi,yi0rg,yi0pk,
                 &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
         pdf_minus(&tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
         pdf_inc(&rpatb,&rpaix,&rpami,&rparg,&rpapk,
                 tmptb,tmpix,tmpmi,tmprg,tmppk);
         // ipa = ipa + xi0*yr0 + xr0*yi0;
         pdf_mul(xi0tb,xi0ix,xi0mi,xi0rg,xi0pk,yr0tb,yr0ix,yr0mi,yr0rg,yr0pk,
                 &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
         pdf_inc(&ipatb,&ipaix,&ipami,&iparg,&ipapk,
                 tmptb,tmpix,tmpmi,tmprg,tmppk);
         pdf_mul(xr0tb,xr0ix,xr0mi,xr0rg,xr0pk,yi0tb,yi0ix,yi0mi,yi0rg,yi0pk,
                 &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
         pdf_inc(&ipatb,&ipaix,&ipami,&iparg,&ipapk,
                 tmptb,tmpix,tmpmi,tmprg,tmppk);
      }
      zretb[k] = rpatb; zreix[k] = rpaix; zremi[k] = rpami;
      zrerg[k] = rparg; zrepk[k] = rpapk;
      zimtb[k] = ipatb; zimix[k] = ipaix; zimmi[k] = ipami;
      zimrg[k] = iparg; zimpk[k] = ipapk;
   }
}
