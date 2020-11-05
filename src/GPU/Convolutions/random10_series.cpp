// The file random10_series.cpp defines functions specified
// in random10_series.h.

#include "deca_double_functions.h"
#include "random10_vectors.h"
#include "random10_series.h"

void dbl10_exponentials
 ( int deg,
   double xrtb, double xrix, double xrmi, double xrrg, double xrpk, 
   double xltb, double xlix, double xlmi, double xlrg, double xlpk, 
   double *pluxrtb, double *pluxrix, double *pluxrmi, double *pluxrrg,
   double *pluxrpk, double *pluxltb, double *pluxlix, double *pluxlmi,
   double *pluxlrg, double *pluxlpk,
   double *minxrtb, double *minxrix, double *minxrmi, double *minxrrg,
   double *minxrpk, double *minxltb, double *minxlix, double *minxlmi,
   double *minxlrg, double *minxlpk )
{
   double frtb,frix,frmi,frrg,frpk,fltb,flix,flmi,flrg,flpk;

   pluxrtb[0] = 1.0; pluxrix[0] = 0.0; pluxrmi[0] = 0.0;
   pluxrrg[0] = 0.0; pluxrpk[0] = 0.0;
   pluxltb[0] = 0.0; pluxlix[0] = 0.0; pluxlmi[0] = 0.0;
   pluxlrg[0] = 0.0; pluxlpk[0] = 0.0;
   minxrtb[0] = 1.0; minxrix[0] = 0.0; minxrmi[0] = 0.0;
   minxrrg[0] = 0.0; minxrpk[0] = 0.0;
   minxltb[0] = 0.0; minxlix[0] = 0.0; minxlmi[0] = 0.0;
   minxlrg[0] = 0.0; minxlpk[0] = 0.0;
   pluxrtb[1] = xrtb; minxrtb[1] = -xrtb;
   pluxrix[1] = xrix; minxrix[1] = -xrix;
   pluxrmi[1] = xrmi; minxrmi[1] = -xrmi;
   pluxrrg[1] = xrrg; minxrrg[1] = -xrrg;
   pluxrpk[1] = xrpk; minxrpk[1] = -xrpk;
   pluxltb[1] = xltb; minxltb[1] = -xltb;
   pluxlix[1] = xlix; minxlix[1] = -xlix;
   pluxlmi[1] = xlmi; minxlmi[1] = -xlmi;
   pluxlrg[1] = xlrg; minxlrg[1] = -xlrg;
   pluxlpk[1] = xlpk; minxlpk[1] = -xlpk;

   for(int k=2; k<=deg; k++)
   {
      daf_mul(pluxrtb[k-1],pluxrix[k-1],pluxrmi[k-1],pluxrrg[k-1],pluxrpk[k-1],
              pluxltb[k-1],pluxlix[k-1],pluxlmi[k-1],pluxlrg[k-1],pluxlpk[k-1],
              xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
              &pluxrtb[k],&pluxrix[k],&pluxrmi[k],&pluxrrg[k],&pluxrpk[k],
              &pluxltb[k],&pluxlix[k],&pluxlmi[k],&pluxlrg[k],&pluxlpk[k]);
      // x[k] = x[k-1]*r
      daf_mul(minxrtb[k-1],minxrix[k-1],minxrmi[k-1],minxrrg[k-1],minxrpk[k-1],
              minxltb[k-1],minxlix[k-1],minxlmi[k-1],minxlrg[k-1],minxlpk[k-1],
              -xrtb,-xrix,-xrmi,-xrrg,-xrpk,
              -xltb,-xlix,-xlmi,-xlrg,-xlpk,
              &minxrtb[k],&minxrix[k],&minxrmi[k],&minxrrg[k],&minxrpk[k],
              &minxltb[k],&minxlix[k],&minxlmi[k],&minxlrg[k],&minxlpk[k]); 
      // y[k] = y[k-1]*(-r);
      frtb = (double) k;
                  frix = 0.0; frmi = 0.0; frrg = 0.0; frpk = 0.0;
      fltb = 0.0; flix = 0.0; flmi = 0.0; flrg = 0.0; flpk = 0.0;
      daf_div(pluxrtb[k],pluxrix[k],pluxrmi[k],pluxrrg[k],pluxrpk[k],
              pluxltb[k],pluxlix[k],pluxlmi[k],pluxlrg[k],pluxlpk[k],
              frtb,frix,frmi,frrg,frpk,fltb,flix,flmi,flrg,flpk,
              &pluxrtb[k],&pluxrix[k],&pluxrmi[k],&pluxrrg[k],&pluxrpk[k],
              &pluxltb[k],&pluxlix[k],&pluxlmi[k],&pluxlrg[k],&pluxlpk[k]);
      daf_div(minxrtb[k],minxrix[k],minxrmi[k],minxrrg[k],minxrpk[k],
              minxltb[k],minxlix[k],minxlmi[k],minxlrg[k],minxlpk[k],
              frtb,frix,frmi,frrg,frpk,fltb,flix,flmi,flrg,flpk,
              &minxrtb[k],&minxrix[k],&minxrmi[k],&minxrrg[k],&minxrpk[k],
              &minxltb[k],&minxlix[k],&minxlmi[k],&minxlrg[k],&minxlpk[k]);
   }
}

void random_dbl10_exponentials
 ( int deg,
   double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *pluxrtb, double *pluxrix, double *pluxrmi, double *pluxrrg,
   double *pluxrpk, double *pluxltb, double *pluxlix, double *pluxlmi,
   double *pluxlrg, double *pluxlpk,
   double *minxrtb, double *minxrix, double *minxrmi, double *minxrrg,
   double *minxrpk, double *minxltb, double *minxlix, double *minxlmi,
   double *minxlrg, double *minxlpk )
{
   random_deca_double
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk);

   dbl10_exponentials
      (deg,*xrtb,*xrix,*xrmi,*xrrg,*xrpk,*xltb,*xlix,*xlmi,*xlrg,*xlpk,
           pluxrtb,pluxrix,pluxrmi,pluxrrg,pluxrpk,
           pluxltb,pluxlix,pluxlmi,pluxlrg,pluxlpk,
           minxrtb,minxrix,minxrmi,minxrrg,minxrpk,
           minxltb,minxlix,minxlmi,minxlrg,minxlpk);
}

void cmplx10_exponentials
 ( int deg,
   double xrertb, double xrerix, double xrermi, double xrerrg, double xrerpk,
   double xreltb, double xrelix, double xrelmi, double xrelrg, double xrelpk,
   double ximrtb, double ximrix, double ximrmi, double ximrrg, double ximrpk,
   double ximltb, double ximlix, double ximlmi, double ximlrg, double ximlpk,
   double *pluxrertb, double *pluxrerix, double *pluxrermi, double *pluxrerrg,
   double *pluxrerpk, double *pluxreltb, double *pluxrelix, double *pluxrelmi,
   double *pluxrelrg, double *pluxrelpk,
   double *pluximrtb, double *pluximrix, double *pluximrmi, double *pluximrrg,
   double *pluximrpk, double *pluximltb, double *pluximlix, double *pluximlmi,
   double *pluximlrg, double *pluximlpk,
   double *minxrertb, double *minxrerix, double *minxrermi, double *minxrerrg,
   double *minxrerpk, double *minxreltb, double *minxrelix, double *minxrelmi,
   double *minxrelrg, double *minxrelpk,
   double *minximrtb, double *minximrix, double *minximrmi, double *minximrrg,
   double *minximrpk, double *minximltb, double *minximlix, double *minximlmi,
   double *minximlrg, double *minximlpk )
{
   double tmprtb,tmprix,tmprmi,tmprrg,tmprpk;
   double tmpltb,tmplix,tmplmi,tmplrg,tmplpk;

   pluxrertb[0] = 1.0; pluxrerix[0] = 0.0; pluxrermi[0] = 0.0;
   pluxrerrg[0] = 0.0; pluxrerpk[0] = 0.0;
   pluxreltb[0] = 0.0; pluxrelix[0] = 0.0; pluxrelmi[0] = 0.0;
   pluxrelrg[0] = 0.0; pluxrelpk[0] = 0.0;
   minxrertb[0] = 1.0; minxrerix[0] = 0.0; minxrermi[0] = 0.0;
   minxrerrg[0] = 0.0; minxrerpk[0] = 0.0;
   minxreltb[0] = 0.0; minxrelix[0] = 0.0; minxrelmi[0] = 0.0;
   minxrelrg[0] = 0.0; minxrelpk[0] = 0.0;
   pluximrtb[0] = 0.0; pluximrix[0] = 0.0; pluximrmi[0] = 0.0;
   pluximrrg[0] = 0.0; pluximrpk[0] = 0.0;
   pluximltb[0] = 0.0; pluximlix[0] = 0.0; pluximlmi[0] = 0.0;
   pluximlrg[0] = 0.0; pluximlpk[0] = 0.0;
   minximrtb[0] = 0.0; minximrix[0] = 0.0; minximrmi[0] = 0.0;
   minximrrg[0] = 0.0; minximrpk[0] = 0.0;
   minximltb[0] = 0.0; minximlix[0] = 0.0; minximlmi[0] = 0.0;
   minximlrg[0] = 0.0; minximlpk[0] = 0.0;
   pluxrertb[1] = xrertb; pluxrerix[1] = xrerix; pluxrermi[1] = xrermi;
   pluxrerrg[1] = xrerrg; pluxrerpk[1] = xrerpk;
   pluxreltb[1] = xreltb; pluxrelix[1] = xrelix; pluxrelmi[1] = xrelmi;
   pluxrelrg[1] = xrelrg; pluxrelpk[1] = xrelpk;
   pluximrtb[1] = ximrtb; pluximrix[1] = ximrix; pluximrmi[1] = ximrmi;
   pluximrrg[1] = ximrrg; pluximrpk[1] = ximrpk;
   pluximltb[1] = ximltb; pluximlix[1] = ximlix; pluximlmi[1] = ximlmi;
   pluximlrg[1] = ximlrg; pluximlpk[1] = ximlpk;
   minxrertb[1] = -xrertb; minxrerix[1] = -xrerix; minxrermi[1] = -xrermi;
   minxrerrg[1] = -xrerrg; minxrerpk[1] = -xrerpk;
   minxreltb[1] = -xreltb; minxrelix[1] = -xrelix; minxrelmi[1] = -xrelmi;
   minxrelrg[1] = -xrelrg; minxrelpk[1] = -xrelpk;
   minximrtb[1] = -ximrtb; minximrix[1] = -ximrix; minximrmi[1] = -ximrmi;
   minximrrg[1] = -ximrrg; minximrpk[1] = -ximrpk;
   minximltb[1] = -ximltb; minximlix[1] = -ximlix; minximlmi[1] = -ximlmi;
   minximlrg[1] = -ximlrg; minximlpk[1] = -ximlpk;

   for(int k=2; k<=deg; k++)
   {
      // xre[k] = (xre[k-1]*cr - xim[k-1]*sr)/k;
      daf_mul(pluxrertb[k-1],pluxrerix[k-1],pluxrermi[k-1],pluxrerrg[k-1],
              pluxrerpk[k-1],pluxreltb[k-1],pluxrelix[k-1],pluxrelmi[k-1],
              pluxrelrg[k-1],pluxrelpk[k-1],
              xrertb,xrerix,xrermi,xrerrg,xrerpk,
              xreltb,xrelix,xrelmi,xrelrg,xrelpk,
              &pluxrertb[k],&pluxrerix[k],&pluxrermi[k],&pluxrerrg[k],
              &pluxrerpk[k],&pluxreltb[k],&pluxrelix[k],&pluxrelmi[k],
              &pluxrelrg[k],&pluxrelpk[k]);
      daf_mul(pluximrtb[k-1],pluximrix[k-1],pluximrmi[k-1],pluximrrg[k-1],
              pluximrpk[k-1],pluximltb[k-1],pluximlix[k-1],pluximlmi[k-1],
              pluximlrg[k-1],pluximlpk[k-1],
              ximrtb,ximrix,ximrmi,ximrrg,ximrpk,
              ximltb,ximlix,ximlmi,ximlrg,ximlpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&pluxrertb[k],&pluxrerix[k],&pluxrermi[k],&pluxrerrg[k],
              &pluxrerpk[k],&pluxreltb[k],&pluxrelix[k],&pluxrelmi[k],
              &pluxrelrg[k],&pluxrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(pluxrertb[k],pluxrerix[k],pluxrermi[k],pluxrerrg[k],
              pluxrerpk[k],pluxreltb[k],pluxrelix[k],pluxrelmi[k],
              pluxrelrg[k],pluxrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &pluxrertb[k],&pluxrerix[k],&pluxrermi[k],&pluxrerrg[k],
              &pluxrerpk[k],&pluxreltb[k],&pluxrelix[k],&pluxrelmi[k],
              &pluxrelrg[k],&pluxrelpk[k]);
      // xim[k] = (xre[k-1]*sr + xim[k-1]*cr)/k;
      daf_mul(pluxrertb[k-1],pluxrerix[k-1],pluxrermi[k-1],pluxrerrg[k-1],
              pluxrerpk[k-1],pluxreltb[k-1],pluxrelix[k-1],pluxrelmi[k-1],
              pluxrelrg[k-1],pluxrelpk[k-1],
              ximrtb,ximrix,ximrmi,ximrrg,ximrpk,
              ximltb,ximlix,ximlmi,ximlrg,ximlpk,
              &pluximrtb[k],&pluximrix[k],&pluximrmi[k],&pluximrrg[k],
              &pluximrpk[k],&pluximltb[k],&pluximlix[k],&pluximlmi[k],
              &pluximlrg[k],&pluximlpk[k]);
      daf_mul(pluximrtb[k-1],pluximrix[k-1],pluximrmi[k-1],pluximrrg[k-1],
              pluximrpk[k-1],pluximltb[k-1],pluximlix[k-1],pluximlmi[k-1],
              pluximlrg[k-1],pluximlpk[k-1],
              xrertb,xrerix,xrermi,xrerrg,xrerpk,
              xreltb,xrelix,xrelmi,xrelrg,xrelpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&pluximrtb[k],&pluximrix[k],&pluximrmi[k],&pluximrrg[k],
              &pluximrpk[k],&pluximltb[k],&pluximlix[k],&pluximlmi[k],
              &pluximlrg[k],&pluximlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(pluximrtb[k],pluximrix[k],pluximrmi[k],pluximrrg[k],
              pluximrpk[k],pluximltb[k],pluximlix[k],pluximlmi[k],
              pluximlrg[k],pluximlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &pluximrtb[k],&pluximrix[k],&pluximrmi[k],&pluximrrg[k],
              &pluximrpk[k],&pluximltb[k],&pluximlix[k],&pluximlmi[k],
              &pluximlrg[k],&pluximlpk[k]);
      // minxre[k] = (minxre[k-1]*(-cr) - minxim[k-1]*(-sr))/k;
      daf_mul(minxrertb[k-1],minxrerix[k-1],minxrermi[k-1],minxrerrg[k-1],
              minxrerpk[k-1],minxreltb[k-1],minxrelix[k-1],minxrelmi[k-1],
              minxrelrg[k-1],minxrelpk[k-1],
              -xrertb,-xrerix,-xrermi,-xrerrg,-xrerpk,
              -xreltb,-xrelix,-xrelmi,-xrelrg,-xrelpk,
              &minxrertb[k],&minxrerix[k],&minxrermi[k],&minxrerrg[k],
              &minxrerpk[k],&minxreltb[k],&minxrelix[k],&minxrelmi[k],
              &minxrelrg[k],&minxrelpk[k]);
      daf_mul(minximrtb[k-1],minximrix[k-1],minximrmi[k-1],minximrrg[k-1],
              minximrpk[k-1],minximltb[k-1],minximlix[k-1],minximlmi[k-1],
              minximlrg[k-1],minximlpk[k-1],
              -ximrtb,-ximrix,-ximrmi,-ximrrg,-ximrpk,
              -ximltb,-ximlix,-ximlmi,-ximlrg,-ximlpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&minxrertb[k],&minxrerix[k],&minxrermi[k],&minxrerrg[k],
              &minxrerpk[k],&minxreltb[k],&minxrelix[k],&minxrelmi[k],
              &minxrelrg[k],&minxrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(minxrertb[k],minxrerix[k],minxrermi[k],minxrerrg[k],
              minxrerpk[k],minxreltb[k],minxrelix[k],minxrelmi[k],
              minxrelrg[k],minxrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &minxrertb[k],&minxrerix[k],&minxrermi[k],&minxrerrg[k],
              &minxrerpk[k],&minxreltb[k],&minxrelix[k],&minxrelmi[k],
              &minxrelrg[k],&minxrelpk[k]);
      // minxim[k] = (minxre[k-1]*(-sr) + minxim[k-1]*(-cr))/k;
      daf_mul(minxrertb[k-1],minxrerix[k-1],minxrermi[k-1],minxrerrg[k-1],
              minxrerpk[k-1],minxreltb[k-1],minxrelix[k-1],minxrelmi[k-1],
              minxrelrg[k-1],minxrelpk[k-1],
              -ximrtb,-ximrix,-ximrmi,-ximrrg,-ximrpk,
              -ximltb,-ximlix,-ximlmi,-ximlrg,-ximlpk,
              &minximrtb[k],&minximrix[k],&minximrmi[k],&minximrrg[k],
              &minximrpk[k],&minximltb[k],&minximlix[k],&minximlmi[k],
              &minximlrg[k],&minximlpk[k]);
      daf_mul(minximrtb[k-1],minximrix[k-1],minximrmi[k-1],minximrrg[k-1],
              minximrpk[k-1],minximltb[k-1],minximlix[k-1],minximlmi[k-1],
              minximlrg[k-1],minximlpk[k-1],
              -xrertb,-xrerix,-xrermi,-xrerrg,-xrerpk,
              -xreltb,-xrelix,-xrelmi,-xrelrg,-xrelpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&minximrtb[k],&minximrix[k],&minximrmi[k],&minximrrg[k],
              &minximrpk[k],&minximltb[k],&minximlix[k],&minximlmi[k],
              &minximlrg[k],&minximlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(minximrtb[k],minximrix[k],minximrmi[k],minximrrg[k],
              minximrpk[k],minximltb[k],minximlix[k],minximlmi[k],
              minximlrg[k],minximlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &minximrtb[k],&minximrix[k],&minximrmi[k],&minximrrg[k],
              &minximrpk[k],&minximltb[k],&minximlix[k],&minximlmi[k],
              &minximlrg[k],&minximlpk[k]);
   }
}

void random_cmplx10_exponentials
 ( int deg,
   double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk,
   double *ximrtb, double *ximrix, double *ximrmi, double *ximrrg,
   double *ximrpk, double *ximltb, double *ximlix, double *ximlmi,
   double *ximlrg, double *ximlpk,
   double *pluxrertb, double *pluxrerix, double *pluxrermi, double *pluxrerrg,
   double *pluxrerpk, double *pluxreltb, double *pluxrelix, double *pluxrelmi,
   double *pluxrelrg, double *pluxrelpk,
   double *pluximrtb, double *pluximrix, double *pluximrmi, double *pluximrrg,
   double *pluximrpk, double *pluximltb, double *pluximlix, double *pluximlmi,
   double *pluximlrg, double *pluximlpk,
   double *minxrertb, double *minxrerix, double *minxrermi, double *minxrerrg,
   double *minxrerpk, double *minxreltb, double *minxrelix, double *minxrelmi,
   double *minxrelrg, double *minxrelpk,
   double *minximrtb, double *minximrix, double *minximrmi, double *minximrrg,
   double *minximrpk, double *minximltb, double *minximlix, double *minximlmi,
   double *minximlrg, double *minximlpk )
{
   double tmprtb,tmprix,tmprmi,tmprrg,tmprpk;
   double tmpltb,tmplix,tmplmi,tmplrg,tmplpk;

   random_deca_double
      (xrertb,xrerix,xrermi,xrerrg,xrerpk,
       xreltb,xrelix,xrelmi,xrelrg,xrelpk);                 // cos(a)

   daf_sqr(*xrertb,*xrerix,*xrermi,*xrerrg,*xrerpk,
           *xreltb,*xrelix,*xrelmi,*xrelrg,*xrelpk,
           &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
           &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);        // cos^2(a)
   daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
             &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);      // -cos^2(a)
   daf_inc_d(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
             &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk,1.0);  // 1-cos^2(a)
   daf_sqrt(tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
            tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
            ximrtb,ximrix,ximrmi,ximrrg,ximrpk,
            ximltb,ximlix,ximlmi,ximlrg,ximlpk);            // sin is sqrt

   cmplx10_exponentials
      (deg,*xrertb,*xrerix,*xrermi,*xrerrg,*xrerpk,
           *xreltb,*xrelix,*xrelmi,*xrelrg,*xrelpk,
           *ximrtb,*ximrix,*ximrmi,*ximrrg,*ximrpk,
           *ximltb,*ximlix,*ximlmi,*ximlrg,*ximlpk,
           pluxrertb,pluxrerix,pluxrermi,pluxrerrg,pluxrerpk,
           pluxreltb,pluxrelix,pluxrelmi,pluxrelrg,pluxrelpk,
           pluximrtb,pluximrix,pluximrmi,pluximrrg,pluximrpk,
           pluximltb,pluximlix,pluximlmi,pluximlrg,pluximlpk,
           minxrertb,minxrerix,minxrermi,minxrerrg,minxrerpk,
           minxreltb,minxrelix,minxrelmi,minxrelrg,minxrelpk,
           minximrtb,minximrix,minximrmi,minximrrg,minximrpk,
           minximltb,minximlix,minximlmi,minximlrg,minximlpk);
}
