// The file random5_series.cpp defines functions specified
// in random5_series.h.

#include "random5_vectors.h"
#include "penta_double_functions.h"

void dbl5_exponentials
 ( int deg, double xtb, double xix, double xmi, double xrg, double xpk, 
   double *pluxtb, double *pluxix, double *pluxmi, double *pluxrg,
   double *pluxpk, double *minxtb, double *minxix, double *minxmi,
   double *minxrg, double *minxpk )
{
   double ftb,fix,fmi,frg,fpk;

   pluxtb[0] = 1.0;  pluxix[0] = 0.0;  pluxmi[0] = 0.0;
   pluxrg[0] = 0.0;  pluxpk[0] = 0.0;
   minxtb[0] = 1.0;  minxix[0] = 0.0;  minxmi[0] = 0.0;
   minxrg[0] = 0.0;  minxpk[0] = 0.0;
   pluxtb[1] = xtb;  pluxix[1] = xix;  pluxmi[1] = xmi;
   pluxrg[1] = xrg;  pluxpk[1] = xpk;
   minxtb[1] = -xtb; minxix[1] = -xix; minxmi[1] = -xmi;
   minxrg[1] = -xrg; minxpk[1] = -xpk;

   for(int k=2; k<=deg; k++)
   {
      pdf_mul(pluxtb[k-1],pluxix[k-1],pluxmi[k-1],pluxrg[k-1],pluxpk[k-1],
              xtb,xix,xmi,xrg,xpk,
              &pluxtb[k],&pluxix[k],&pluxmi[k],&pluxrg[k],&pluxpk[k]);
      // x[k] = x[k-1]*r
      pdf_mul(minxtb[k-1],minxix[k-1],minxmi[k-1],minxrg[k-1],minxpk[k-1],
              -xtb,-xix,-xmi,-xrg,-xpk,
              &minxtb[k],&minxix[k],&minxmi[k],&minxrg[k],&minxpk[k]); 
      // y[k] = y[k-1]*(-r);
      ftb = (double) k; fix = 0.0; fmi = 0.0; frg = 0.0; fpk = 0.0;
      pdf_div(pluxtb[k],pluxix[k],pluxmi[k],pluxrg[k],pluxpk[k],
              ftb,fix,fmi,frg,fpk,
              &pluxtb[k],&pluxix[k],&pluxmi[k],&pluxrg[k],&pluxpk[k]);
      pdf_div(minxtb[k],minxix[k],minxmi[k],minxrg[k],minxpk[k],
              ftb,fix,fmi,frg,fpk,
              &minxtb[k],&minxix[k],&minxmi[k],&minxrg[k],&minxpk[k]);
   }
}

void random_dbl5_exponentials
 ( int deg, double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *pluxtb, double *pluxix, double *pluxmi, double *pluxrg,
   double *pluxpk, double *minxtb, double *minxix, double *minxmi,
   double *minxrg, double *minxpk )
{
   random_penta_double(xtb,xix,xmi,xrg,xpk);

   dbl5_exponentials
      (deg,*xtb,*xix,*xmi,*xrg,*xpk,
       pluxtb,pluxix,pluxmi,pluxrg,pluxpk,
       minxtb,minxix,minxmi,minxrg,minxpk);
}

void cmplx5_exponentials
 ( int deg,
   double xretb, double xreix, double xremi, double xrerg, double xrepk,
   double ximtb, double ximix, double ximmi, double ximrg, double ximpk,
   double *pluxretb, double *pluxreix, double *pluxremi, double *pluxrerg,
   double *pluxrepk, double *pluximtb, double *pluximix, double *pluximmi,
   double *pluximrg, double *pluximpk,
   double *minxretb, double *minxreix, double *minxremi, double *minxrerg,
   double *minxrepk, double *minximtb, double *minximix, double *minximmi,
   double *minximrg, double *minximpk )
{
   double tmptb,tmpix,tmpmi,tmprg,tmppk;

   pluxretb[0] = 1.0; pluxreix[0] = 0.0; pluxremi[0] = 0.0;
   pluxrerg[0] = 0.0; pluxrepk[0] = 0.0;
   minxretb[0] = 1.0; minxreix[0] = 0.0; minxremi[0] = 0.0;
   minxrerg[0] = 0.0; minxrepk[0] = 0.0;
   pluximtb[0] = 0.0; pluximix[0] = 0.0; pluximmi[0] = 0.0;
   pluximrg[0] = 0.0; pluximpk[0] = 0.0;
   minximtb[0] = 0.0; minximix[0] = 0.0; minximmi[0] = 0.0;
   pluximrg[0] = 0.0; minximpk[0] = 0.0;
   pluxretb[1] = xretb; pluxreix[1] = xreix; pluxremi[1] = xremi;
   pluxrerg[1] = xrerg; pluxrepk[1] = xrepk;
   pluximtb[1] = ximtb; pluximix[1] = ximix; pluximmi[1] = ximmi;
   pluximrg[1] = ximrg; pluximpk[1] = ximpk;
   minxretb[1] = -xretb; minxreix[1] = -xreix; minxremi[1] = -xremi;
   minxrerg[1] = -xrerg; minxrepk[1] = -xrepk;
   minximtb[1] = -ximtb; minximix[1] = -ximix; minximmi[1] = -ximmi;
   minximrg[1] = -ximrg; minximpk[1] = -ximpk;

   for(int k=2; k<=deg; k++)
   {
      // pluxre[k] = (pluxre[k-1]*cr - pluxim[k-1]*sr)/k;
      pdf_mul(pluxretb[k-1],pluxreix[k-1],pluxremi[k-1],
              pluxrerg[k-1],pluxrepk[k-1],
              xretb,xreix,xremi,xrerg,xrepk,
              &pluxretb[k],&pluxreix[k],&pluxremi[k],
              &pluxrerg[k],&pluxrepk[k]);
      pdf_mul(pluximtb[k-1],pluximix[k-1],pluximmi[k-1],
              pluximrg[k-1],pluximpk[k-1],
              ximtb,ximix,ximmi,ximrg,ximpk,
              &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_minus(&tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_inc(&pluxretb[k],&pluxreix[k],&pluxremi[k],&pluxrerg[k],
              &pluxrepk[k],tmptb,tmpix,tmpmi,tmprg,tmppk);
      tmptb = (double) k; tmpix = 0.0; tmpmi = 0.0; tmprg = 0.0; tmppk = 0.0;
      pdf_div(pluxretb[k],pluxreix[k],pluxremi[k],pluxrerg[k],pluxrepk[k],
              tmptb,tmpix,tmpmi,tmprg,tmppk,
              &pluxretb[k],&pluxreix[k],&pluxremi[k],&pluxrerg[k],
              &pluxrepk[k]);
      // pluxim[k] = (pluxre[k-1]*sr + pluxim[k-1]*cr)/k;
      pdf_mul(pluxretb[k-1],pluxreix[k-1],pluxremi[k-1],pluxrerg[k-1],
              pluxrepk[k-1],ximtb,ximix,ximmi,ximrg,ximpk,
              &pluximtb[k],&pluximix[k],&pluximmi[k],&pluximrg[k],
              &pluximpk[k]);
      pdf_mul(pluximtb[k-1],pluximix[k-1],pluximmi[k-1],pluximrg[k-1],
              pluximpk[k-1],xretb,xreix,xremi,xrerg,xrepk,
              &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_inc(&pluximtb[k],&pluximix[k],&pluximmi[k],&pluximrg[k],
              &pluximpk[k],tmptb,tmpix,tmpmi,tmprg,tmppk);
      tmptb = (double) k; tmpix = 0.0; tmpmi = 0.0; tmprg = 0.0; tmppk = 0.0;
      pdf_div(pluximtb[k],pluximix[k],pluximmi[k],pluximrg[k],pluximpk[k],
              tmptb,tmpix,tmpmi,tmprg,tmppk,&pluximtb[k],&pluximix[k],
              &pluximmi[k],&pluximrg[k],&pluximpk[k]);
      // minxre[k] = (minxre[k-1]*(-cr) - minxim[k-1]*(-sr))/k;
      pdf_mul(minxretb[k-1],minxreix[k-1],minxremi[k-1],minxrerg[k-1],
              minxrepk[k-1],-xretb,-xreix,-xremi,-xrerg,-xrepk,
              &minxretb[k],&minxreix[k],&minxremi[k],&minxrerg[k],
              &minxrepk[k]);
      pdf_mul(minximtb[k-1],minximix[k-1],minximmi[k-1],minximrg[k-1],
              minximpk[k-1],-ximtb,-ximix,-ximmi,-ximrg,-ximpk,
              &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_minus(&tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_inc(&minxretb[k],&minxreix[k],&minxremi[k],&minxrerg[k],
              &minxrepk[k],tmptb,tmpix,tmpmi,tmprg,tmppk);
      tmptb = (double) k; tmpix = 0.0; tmpmi = 0.0; tmprg = 0.0; tmppk = 0.0;
      pdf_div(minxretb[k],minxreix[k],minxremi[k],minxrerg[k],minxrepk[k],
              tmptb,tmpix,tmpmi,tmprg,tmppk,&minxretb[k],&minxreix[k],
              &minxremi[k],&minxrerg[k],&minxrepk[k]);
      // minxim[k] = (minxre[k-1]*(-sr) + minxim[k-1]*(-cr))/k;
      pdf_mul(minxretb[k-1],minxreix[k-1],minxremi[k-1],minxrerg[k-1],
              minxrepk[k-1],-ximtb,-ximix,-ximmi,-ximrg,-ximpk,
              &minximtb[k],&minximix[k],&minximmi[k],&minximrg[k],
              &minximpk[k]);
      pdf_mul(minximtb[k-1],minximix[k-1],minximmi[k-1],minximrg[k-1],
              minximpk[k-1],-xretb,-xreix,-xremi,-xrerg,-xrepk,
              &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);
      pdf_inc(&minximtb[k],&minximix[k],&minximmi[k],&minximrg[k],
              &minximpk[k],tmptb,tmpix,tmpmi,tmprg,tmppk);
      tmptb = (double) k; tmpix = 0.0; tmpmi = 0.0; tmprg = 0.0; tmppk = 0.0;
      pdf_div(minximtb[k],minximix[k],minximmi[k],minximrg[k],minximpk[k],
              tmptb,tmpix,tmpmi,tmprg,tmppk,
              &minximtb[k],&minximix[k],&minximmi[k],&minximrg[k],
              &minximpk[k]);
   }
}

void random_cmplx5_exponentials
 ( int deg,
   double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *pluxretb, double *pluxreix, double *pluxremi, double *pluxrerg,
   double *pluxrepk, double *pluximtb, double *pluximix, double *pluximmi,
   double *pluximrg, double *pluximpk,
   double *minxretb, double *minxreix, double *minxremi, double *minxrerg,
   double *minxrepk, double *minximtb, double *minximix, double *minximmi,
   double *minximrg, double *minximpk )
{
   double tmptb,tmpix,tmpmi,tmprg,tmppk;

   random_penta_double(xretb,xreix,xremi,xrerg,xrepk);      // cos(a)

   pdf_sqr(*xretb,*xreix,*xremi,*xrerg,*xrepk,
           &tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);             // cos^2(a)
   pdf_minus(&tmptb,&tmpix,&tmpmi,&tmprg,&tmppk);           // -cos^2(a)
   pdf_inc_d(&tmptb,&tmpix,&tmpmi,&tmprg,&tmppk,1.0);       // 1-cos^2(a)
   pdf_sqrt(tmptb,tmpix,tmpmi,tmprg,tmppk,
            ximtb,ximix,ximmi,ximrg,ximpk);                 // sin is sqrt

   cmplx5_exponentials
      (deg,*xretb,*xreix,*xremi,*xrerg,*xrepk,
           *ximtb,*ximix,*ximmi,*ximrg,*ximpk,
           pluxretb,pluxreix,pluxremi,pluxrerg,pluxrepk,
           pluximtb,pluximix,pluximmi,pluximrg,pluximpk,
           minxretb,minxreix,minxremi,minxrerg,minxrepk,
           minximtb,minximix,minximmi,minximrg,minximpk);
}
