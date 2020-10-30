/* The file dbl10_convolutions_host.cpp defines functions
 * specified in dbl10_convolutions_host.h. */

#include "dbl10_convolutions_host.h"
#include "deca_double_functions.h"

void CPU_dbl10_product
 ( int deg,
   double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk )
{
   int idx;
   double prtb,prix,prmi,prrg,prpk;
   double pltb,plix,plmi,plrg,plpk;

   daf_mul(xrtb[0],xrix[0],xrmi[0],xrrg[0],xrpk[0],
           xltb[0],xlix[0],xlmi[0],xlrg[0],xlpk[0],
           yrtb[0],yrix[0],yrmi[0],yrrg[0],yrpk[0],
           yltb[0],ylix[0],ylmi[0],ylrg[0],ylpk[0],
           &zrtb[0],&zrix[0],&zrmi[0],&zrrg[0],&zrpk[0],
           &zltb[0],&zlix[0],&zlmi[0],&zlrg[0],&zlpk[0]); // z[0] = x[0]*y[0]

   for(int k=1; k<=deg; k++)
   {
      daf_mul(xrtb[0],xrix[0],xrmi[0],xrrg[0],xrpk[0],
              xltb[0],xlix[0],xlmi[0],xlrg[0],xlpk[0],
              yrtb[k],yrix[k],yrmi[k],yrrg[k],yrpk[k],
              yltb[k],ylix[k],ylmi[k],ylrg[k],ylpk[k],
              &zrtb[k],&zrix[k],&zrmi[k],&zrrg[k],&zrpk[k],
              &zltb[k],&zlix[k],&zlmi[k],&zlrg[k],&zlpk[k]);
      // z[k] = x[0]*y[k]
      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         daf_mul(xrtb[i],xrix[i],xrmi[i],xrrg[i],xrpk[i],
                 xltb[i],xlix[i],xlmi[i],xlrg[i],xlpk[i],
                 yrtb[idx],yrix[idx],yrmi[idx],yrrg[idx],yrpk[idx],
                 yltb[idx],ylix[idx],ylmi[idx],ylrg[idx],ylpk[idx],
                 &prtb,&prix,&prmi,&prrg,&prpk,
                 &pltb,&plix,&plmi,&plrg,&plpk);             // x[i]*y[k-i]
         // z[k] = z[k] + x[i]*y[k-i]
         daf_inc(&zrtb[k],&zrix[k],&zrmi[k],&zrrg[k],&zrpk[k],
                 &zltb[k],&zlix[k],&zlmi[k],&zlrg[k],&zlpk[k],
                 prtb,prix,prmi,prrg,prpk,pltb,plix,plmi,plrg,plpk); 
      }
   }
}

void CPU_cmplx10_product
 ( int deg,
   double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk, double *ximrtb, double *ximrix,
   double *ximrmi, double *ximrrg, double *ximrpk, double *ximltb,
   double *ximlix, double *ximlmi, double *ximlrg, double *ximlpk,
   double *yrertb, double *yrerix, double *yrermi, double *yrerrg,
   double *yrerpk, double *yreltb, double *yrelix, double *yrelmi,
   double *yrelrg, double *yrelpk, double *yimrtb, double *yimrix,
   double *yimrmi, double *yimrrg, double *yimrpk, double *yimltb,
   double *yimlix, double *yimlmi, double *yimlrg, double *yimlpk,
   double *zrertb, double *zrerix, double *zrermi, double *zrerrg,
   double *zrerpk, double *zreltb, double *zrelix, double *zrelmi,
   double *zrelrg, double *zrelpk, double *zimrtb, double *zimrix,
   double *zimrmi, double *zimrrg, double *zimrpk, double *zimltb,
   double *zimlix, double *zimlmi, double *zimlrg, double *zimlpk )
{
   double rpartb,rparix,rparmi,rparrg,rparpk;   // high real parts
   double rpaltb,rpalix,rpalmi,rpalrg,rpalpk;   // low real parts
   double ipartb,iparix,iparmi,iparrg,iparpk;   // high imaginary parts
   double ipaltb,ipalix,ipalmi,ipalrg,ipalpk;   // low imaginary parts
   double tmprtb,tmprix,tmprmi,tmprrg,tmprpk;   // high temporary deca double
   double tmpltb,tmplix,tmplmi,tmplrg,tmplpk;   // low temporary deca double
   double xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk;   // high real values in xr 
   double xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk;   // low real values in xr 
   double xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk;   // high imaginary values in xi
   double xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk;   // low imaginary values in xi
   double yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk;   // high real values in yr
   double yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk;   // low real values in yr
   double yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk;   // high imaginary values in yi
   double yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk;   // low imaginary values in yi
   int idx;

   xr0rtb = xrertb[0]; xr0rix = xrerix[0]; xr0rmi = xrermi[0];
   xr0rrg = xrerrg[0]; xr0rpk = xrerpk[0];
   xr0ltb = xreltb[0]; xr0lix = xrelix[0]; xr0lmi = xrelmi[0];
   xr0lrg = xrelrg[0]; xr0lpk = xrelpk[0];
   xi0rtb = ximrtb[0]; xi0rix = ximrix[0]; xi0rmi = ximrmi[0];
   xi0rrg = ximrrg[0]; xi0rpk = ximrpk[0];
   xi0ltb = ximltb[0]; xi0lix = ximlix[0]; xi0lmi = ximlmi[0];
   xi0lrg = ximlrg[0]; xi0lpk = ximlpk[0];
   yr0rtb = yrertb[0]; yr0rix = yrerix[0]; yr0rmi = yrermi[0];
   yr0rrg = yrerrg[0]; yr0rpk = yrerpk[0];
   yr0ltb = yreltb[0]; yr0lix = yrelix[0]; yr0lmi = yrelmi[0];
   yr0lrg = yrelrg[0]; yr0lpk = yrelpk[0];
   yi0rtb = yimrtb[0]; yi0rix = yimrix[0]; yi0rmi = yimrmi[0];
   yi0rrg = yimrrg[0]; yi0rpk = yimrpk[0];
   yi0ltb = yimltb[0]; yi0lix = yimlix[0]; yi0lmi = yimlmi[0];
   yi0lrg = yimlrg[0]; yi0lpk = yimlpk[0];
   // zre[0] = xr0*yr0 - xi0*yi0;
   daf_mul(xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk,
           xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk,
           yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk,
           yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk,
           &zrertb[0],&zrerix[0],&zrermi[0],&zrerrg[0],&zrerpk[0],
           &zreltb[0],&zrelix[0],&zrelmi[0],&zrelrg[0],&zrelpk[0]);
   daf_mul(xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk,
           xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk,
           yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk,
           yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk,
           &rpartb,&rparix,&rparmi,&rparrg,&rparpk,
           &rpaltb,&rpalix,&rpalmi,&rpalrg,&rpalpk);
   daf_minus(&rpartb,&rparix,&rparmi,&rparrg,&rparpk,
             &rpaltb,&rpalix,&rpalmi,&rpalrg,&rpalpk);
   daf_inc(&zrertb[0],&zrerix[0],&zrermi[0],&zrerrg[0],&zrerpk[0],
           &zreltb[0],&zrelix[0],&zrelmi[0],&zrelrg[0],&zrelpk[0],
           rpartb,rparix,rparmi,rparrg,rparpk,
           rpaltb,rpalix,rpalmi,rpalrg,rpalpk);
   // zim[0] = xi0*yr0 + xr0*yi0;
   daf_mul(xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk,
           xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk,
           yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk,
           yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk,
           &zimrtb[0],&zimrix[0],&zimrmi[0],&zimrrg[0],&zimrpk[0],
           &zimltb[0],&zimlix[0],&zimlmi[0],&zimlrg[0],&zimlpk[0]);
   daf_mul(xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk,
           xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk,
           yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk,
           yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk,
           &ipartb,&iparix,&iparmi,&iparrg,&iparpk,
           &ipaltb,&ipalix,&ipalmi,&ipalrg,&ipalpk);
   daf_inc(&zimrtb[0],&zimrix[0],&zimrmi[0],&zimrrg[0],&zimrpk[0],
           &zimltb[0],&zimlix[0],&zimlmi[0],&zimlrg[0],&zimlpk[0],
           ipartb,iparix,iparmi,iparrg,iparpk,
           ipaltb,ipalix,ipalmi,ipalrg,ipalpk);

   for(int k=1; k<=deg; k++)
   {
      xr0rtb = xrertb[0]; xr0rix = xrerix[0]; xr0rmi = xrermi[0];
      xr0rrg = xrerrg[0]; xr0rpk = xrerpk[0];
      xr0ltb = xreltb[0]; xr0lix = xrelix[0]; xr0lmi = xrelmi[0];
      xr0lrg = xrelrg[0]; xr0lpk = xrelpk[0];
      xi0rtb = ximrtb[0]; xi0rix = ximrix[0]; xi0rmi = ximrmi[0];
      xi0rrg = ximrrg[0]; xi0rpk = ximrpk[0];
      xi0ltb = ximltb[0]; xi0lix = ximlix[0]; xi0lmi = ximlmi[0];
      xi0lrg = ximlrg[0]; xi0lpk = ximlpk[0];
      yr0rtb = yrertb[k]; yr0rix = yrerix[k]; yr0rmi = yrermi[k];
      yr0rrg = yrerrg[k]; yr0rpk = yrerpk[k];
      yr0ltb = yreltb[k]; yr0lix = yrelix[k]; yr0lmi = yrelmi[k];
      yr0lrg = yrelrg[k]; yr0lpk = yrelpk[k];
      yi0rtb = yimrtb[k]; yi0rix = yimrix[k]; yi0rmi = yimrmi[k];
      yi0rrg = yimrrg[k]; yi0rpk = yimrpk[k];
      yi0ltb = yimltb[k]; yi0lix = yimlix[k]; yi0lmi = yimlmi[k];
      yi0lrg = yimlrg[k]; yi0lpk = yimlpk[k];
      // rpa = xr0*yr0 - xi0*yi0;
      daf_mul(xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk,
              xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk,
              yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk,
              yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk,
              &rpartb,&rparix,&rparmi,&rparrg,&rparpk,
              &rpaltb,&rpalix,&rpalmi,&rpalrg,&rpalpk);
      daf_mul(xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk,
              xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk,
              yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk,
              yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&rpartb,&rparix,&rparmi,&rparrg,&rparpk,
              &rpaltb,&rpalix,&rpalmi,&rpalrg,&rpalpk,
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      // ipa = xi0*yr0 + xr0*yi0;
      daf_mul(xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk,
              xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk,
              yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk,
              yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk,
              &ipartb,&iparix,&iparmi,&iparrg,&iparpk,
              &ipaltb,&ipalix,&ipalmi,&ipalrg,&ipalpk);
      daf_mul(xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk,
              xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk,
              yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk,
              yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&ipartb,&iparix,&iparmi,&iparrg,&iparpk,
              &ipaltb,&ipalix,&ipalmi,&ipalrg,&ipalpk,
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0rtb = xrertb[i]; xr0rix = xrerix[i]; xr0rmi = xrermi[i];
         xr0rrg = xrerrg[i]; xr0rpk = xrerpk[i];
         xr0ltb = xreltb[i]; xr0lix = xrelix[i]; xr0lmi = xrelmi[i];
         xr0lrg = xrelrg[i]; xr0lpk = xrelpk[i];
         xi0rtb = ximrtb[i]; xi0rix = ximrix[i]; xi0rmi = ximrmi[i];
         xi0rrg = ximrrg[i]; xi0rpk = ximrpk[i];
         xi0ltb = ximltb[i]; xi0lix = ximlix[i]; xi0lmi = ximlmi[i];
         xi0lrg = ximlrg[i]; xi0lpk = ximlpk[i];
         yr0rtb = yrertb[idx]; yr0rix = yrerix[idx]; yr0rmi = yrermi[idx];
         yr0rrg = yrerrg[idx]; yr0rpk = yrerpk[idx];
         yr0ltb = yreltb[idx]; yr0lix = yrelix[idx]; yr0lmi = yrelmi[idx];
         yr0lrg = yrelrg[idx]; yr0lpk = yrelpk[idx];
         yi0rtb = yimrtb[idx]; yi0rix = yimrix[idx]; yi0rmi = yimrmi[idx];
         yi0rrg = yimrrg[idx]; yi0rpk = yimrpk[idx];
         yi0ltb = yimltb[idx]; yi0lix = yimlix[idx]; yi0lmi = yimlmi[idx];
         yi0lrg = yimlrg[idx]; yi0lpk = yimlpk[idx];
         // rpa = rpa + xr0*yr0 - xi0*yi0;
         daf_mul(xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk,
                 xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk,
                 yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk,
                 yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk,
                 &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                 &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
         daf_inc(&rpartb,&rparix,&rparmi,&rparrg,&rparpk,
                 &rpaltb,&rpalix,&rpalmi,&rpalrg,&rpalpk,
                 tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
                 tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
         daf_mul(xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk,
                 xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk,
                 yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk,
                 yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk,
                 &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                 &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
         daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                   &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
         daf_inc(&rpartb,&rparix,&rparmi,&rparrg,&rparpk,
                 &rpaltb,&rpalix,&rpalmi,&rpalrg,&rpalpk,
                 tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
                 tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
         // ipa = ipa + xi0*yr0 + xr0*yi0;
         daf_mul(xi0rtb,xi0rix,xi0rmi,xi0rrg,xi0rpk,
                 xi0ltb,xi0lix,xi0lmi,xi0lrg,xi0lpk,
                 yr0rtb,yr0rix,yr0rmi,yr0rrg,yr0rpk,
                 yr0ltb,yr0lix,yr0lmi,yr0lrg,yr0lpk,
                 &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                 &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
         daf_inc(&ipartb,&iparix,&iparmi,&iparrg,&iparpk,
                 &ipaltb,&ipalix,&ipalmi,&ipalrg,&ipalpk,
                 tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
                 tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
         daf_mul(xr0rtb,xr0rix,xr0rmi,xr0rrg,xr0rpk,
                 xr0ltb,xr0lix,xr0lmi,xr0lrg,xr0lpk,
                 yi0rtb,yi0rix,yi0rmi,yi0rrg,yi0rpk,
                 yi0ltb,yi0lix,yi0lmi,yi0lrg,yi0lpk,
                 &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                 &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
         daf_inc(&ipartb,&iparix,&iparmi,&iparrg,&iparpk,
                 &ipaltb,&ipalix,&ipalmi,&ipalrg,&ipalpk,
                 tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
                 tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      }
      zrertb[k] = rpartb; zrerix[k] = rparix; zrermi[k] = rparmi;
      zrerrg[k] = rparrg; zrerpk[k] = rparpk;
      zreltb[k] = rpaltb; zrelix[k] = rpalix; zrelmi[k] = rpalmi;
      zrelrg[k] = rpalrg; zrelpk[k] = rpalpk;
      zimrtb[k] = ipartb; zimrix[k] = iparix; zimrmi[k] = iparmi;
      zimrrg[k] = iparrg; zimrpk[k] = iparpk;
      zimltb[k] = ipaltb; zimlix[k] = ipalix; zimlmi[k] = ipalmi;
      zimlrg[k] = ipalrg; zimlpk[k] = ipalpk;
   }
}
