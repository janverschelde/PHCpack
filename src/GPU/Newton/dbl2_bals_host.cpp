// The file dbl2_bals_host.cpp defines functions with prototypes in
// the file dbl2_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "dbl2_factorizations.h"
#include "dbl_onenorms_host.h"
#include "dbl2_bals_host.h"

using namespace std;

void CPU_dbl2_qrbs_head
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo,
   bool *zeroQ, bool *noqr, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "  " << mathi[0][i][j] << "  " << matlo[0][i][j];
         }
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhshi[0][i] << "  " << rhslo[0][i] << endl;
   }
   double nrm;
   CPU_dbl_onenorm(dim,rhshi[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm < 1.0e-28)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping of CPU_dbl2_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_dbl2_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl2_factors_houseqr ..." << endl;

      CPU_dbl2_factors_houseqr(dim,dim,mathi[0],matlo[0],Qhi,Qlo,Rhi,Rlo);

      *zeroQ = false;

      CPU_dbl2_factors_qrbs
         (dim,dim,Qhi,Qlo,Rhi,Rlo,rhshi[0],rhslo[0],solhi[0],sollo[0],
          wrkvechi,wrkveclo);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_dbl_onenorm(dim,wrkvechi,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_dbl_onenorm(dim,solhi[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double acchi,acclo;

         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solhi[0][i] << "  " << sollo[0][i] << endl;
         for(int i=0; i<dim; i++)
         {
            wrkvechi[i] = rhshi[0][i];
            wrkveclo[i] = rhslo[0][i];
            for(int j=0; j<dim; j++)
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
            {
               ddf_mul(mathi[0][i][j],matlo[0][i][j],solhi[0][j],sollo[0][j],
                       &acchi,&acclo);
               ddf_dec(&wrkvechi[i],&wrkveclo[i],acchi,acclo);
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
           cout << wrkvechi[i] << "  " << wrkveclo[i] << endl;
      }
   }
}

void CPU_cmplx2_qrbs_head
 ( int dim, int degp1,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo,
   bool *zeroQ, bool *noqr, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            cout << "  " << matrehi[0][i][j]
                 << "  " << matrelo[0][i][j]
                 << "  " << matimhi[0][i][j]
                 << "  " << matimlo[0][i][j];
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhsrehi[0][i] << "  " << rhsrelo[0][i] << endl << "  " 
              << rhsimhi[0][i] << "  " << rhsimlo[0][i] << endl;
   }
   double nrm;
   CPU_cmplx_onenorm(dim,rhsrehi[0],rhsimhi[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm < 1.0e-28)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping of CPU_cmplx2_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_cmplx2_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx2_factors_houseqr ..." << endl;

      CPU_cmplx2_factors_houseqr
         (dim,dim,matrehi[0],matrelo[0],matimhi[0],matimlo[0],
          Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo);

      *zeroQ = false;

      CPU_cmplx2_factors_qrbs
         (dim,dim,Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo,
          rhsrehi[0],rhsrelo[0],rhsimhi[0],rhsimlo[0],
          solrehi[0],solrelo[0],solimhi[0],solimlo[0],
          wrkvecrehi,wrkvecrelo,wrkvecimhi,wrkvecimlo);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_cmplx_onenorm(dim,wrkvecrehi,wrkvecimhi,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_cmplx_onenorm(dim,solrehi[0],solimhi[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double acchi,acclo;

         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solrehi[0][i] << "  " << solrelo[0][i] << endl << "  "
                 << solimhi[0][i] << "  " << solimlo[0][i] << endl;

         for(int i=0; i<dim; i++)
         {
            wrkvecrehi[i] = rhsrehi[0][i]; wrkvecrelo[i] = rhsrelo[0][i];
            wrkvecimhi[i] = rhsimhi[0][i]; wrkvecimlo[i] = rhsimlo[0][i];

            for(int j=0; j<dim; j++)
            {
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
               // zre = matre[0][i][j]*solre[0][j] - matim[0][i][j]*solim[0][j];
               // wrkvecre[i] = wrkvecre[i] + zre;
               ddf_mul(matrehi[0][i][j],matrelo[0][i][j],
                       solrehi[0][j],solrelo[0][j],&acchi,&acclo);
               ddf_inc(&wrkvecrehi[i],&wrkvecrelo[i],acchi,acclo);
               ddf_mul(matimhi[0][i][j],matimlo[0][i][j],
                       solimhi[0][j],solimlo[0][j],&acchi,&acclo);
               ddf_dec(&wrkvecrehi[i],&wrkvecrelo[i],acchi,acclo);
               // zim = matre[0][i][j]*solim[0][j] + matim[0][i][j]*solre[0][j];
               // wrkvecim[i] = wrkvecim[i] + zim;
               ddf_mul(matrehi[0][i][j],matrelo[0][i][j],
                       solimhi[0][j],solimlo[0][j],&acchi,&acclo);
               ddf_inc(&wrkvecimhi[i],&wrkvecimlo[i],acchi,acclo);
               ddf_mul(matimhi[0][i][j],matimlo[0][i][j],
                       solrehi[0][j],solrelo[0][j],&acchi,&acclo);
               ddf_inc(&wrkvecimhi[i],&wrkvecimlo[i],acchi,acclo);
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
            cout << wrkvecrehi[i] << "  " << wrkvecrelo[i] << endl << "  "
                 << wrkvecimhi[i] << "  " << wrkvecimlo[i] << endl;
      }
   }
}

void CPU_dbl2_qrbs_tail
 ( int dim, int degp1, int tailidx, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   double nrm,acchi,acclo;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      double *xihi = solhi[i-1];    // solution to do the update with
      double *xilo = sollo[i-1];

      CPU_dbl_onenorm(dim,xihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-28)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << i-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << i-1 << "] ..." << endl;

         for(int j=i; j<degp1; j++)
         {
            double **Ajhi = mathi[j-i+1]; // always start with A[1]
            double **Ajlo = matlo[j-i+1]; 
            double *wjhi = rhshi[j];      // current right hand side vector
            double *wjlo = rhslo[j];

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  ddf_mul(Ajhi[k][L],Ajlo[k][L],xihi[L],xilo[L],&acchi,&acclo);
                  ddf_dec(&wjhi[k],&wjlo[k],acchi,acclo);
               }
         }
      }
      // compute sol[i] with back substitution
      double *xhi = solhi[i];
      double *xlo = sollo[i];
      double *bhi = rhshi[i];
      double *blo = rhslo[i];

      CPU_dbl_onenorm(dim,bhi,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-28) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++)
         {
            xhi[j] = 0.0;
            xlo[j] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "-> run backsubstitution for x[" << i << "] ..." << endl;

         if(firstbs)
         {
            *newtail = i;
            firstbs = false;
         }
         CPU_dbl2_factors_qrbs
            (dim,dim,Qhi,Qlo,Rhi,Rlo,bhi,blo,xhi,xlo,wrkvechi,wrkveclo);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "Qtb[" << i << "] : "
                    << wrkvechi[i] << "  " << wrkveclo[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xhi[j] << "  " << xlo[j] << endl;
         }
      }
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;
}

void CPU_cmplx2_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo,
   int *upidx, int *bsidx, int *newtail, int vrblvl )
{
   double nrm,acchi,acclo;
   int skipupcnt = 0;
   int skipbscnt = 0;
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1

      double *xirehi = solrehi[i-1]; // solution to do the update with
      double *xirelo = solrelo[i-1]; 
      double *xiimhi = solimhi[i-1]; 
      double *xiimlo = solimlo[i-1]; 

      CPU_cmplx_onenorm(dim,xirehi,xiimhi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-28)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << i-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << i-1 << "] ..." << endl;

         for(int j=i; j<degp1; j++)
         {
            double **Ajrehi = matrehi[j-i+1]; // always start with A[1]
            double **Ajrelo = matrelo[j-i+1];
            double **Ajimhi = matimhi[j-i+1];
            double **Ajimlo = matimlo[j-i+1];
            double *wjrehi = rhsrehi[j]; // current right hand side vector
            double *wjrelo = rhsrelo[j];
            double *wjimhi = rhsimhi[j]; 
            double *wjimlo = rhsimlo[j]; 

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  // zre = Ajre[k][L]*xire[L] - Ajim[k][L]*xiim[L];
                  // wjre[k] = wjre[k] - zre;
                  ddf_mul(Ajrehi[k][L],Ajrelo[k][L],xirehi[L],xirelo[L],
                          &acchi,&acclo);
                  ddf_dec(&wjrehi[k],&wjrelo[k],acchi,acclo);
                  ddf_mul(Ajimhi[k][L],Ajimlo[k][L],xiimhi[L],xiimlo[L],
                          &acchi,&acclo);
                  ddf_inc(&wjrehi[k],&wjrelo[k],acchi,acclo);
                  // zim = Ajre[k][L]*xiim[L] + Ajim[k][L]*xire[L];
                  // wjim[k] = wjim[k] - zim;
                  ddf_mul(Ajrehi[k][L],Ajrelo[k][L],xiimhi[L],xiimlo[L],
                          &acchi,&acclo);
                  ddf_dec(&wjimhi[k],&wjimlo[k],acchi,acclo);
                  ddf_mul(Ajimhi[k][L],Ajimlo[k][L],xirehi[L],xirelo[L],
                          &acchi,&acclo);
                  ddf_dec(&wjimhi[k],&wjimlo[k],acchi,acclo);
               }
         }
      }
      // compute sol[i] with back substitution
      double *xrehi = solrehi[i];
      double *xrelo = solrelo[i];
      double *ximhi = solimhi[i];
      double *ximlo = solimlo[i];
      double *brehi = rhsrehi[i];
      double *brelo = rhsrelo[i];
      double *bimhi = rhsimhi[i];
      double *bimlo = rhsimlo[i];

      CPU_cmplx_onenorm(dim,brehi,bimhi,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-28) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++)
         {
            xrehi[j] = 0.0; xrelo[j] = 0.0;
            ximhi[j] = 0.0; ximlo[j] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "-> run backsubstitution for x[" << i << "] ..." << endl;

         if(firstbs)
         {
            *newtail = i;
            firstbs = false;
         }
         CPU_cmplx2_factors_qrbs
            (dim,dim,Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo,
             brehi,brelo,bimhi,bimlo,xrehi,xrelo,ximhi,ximlo,
             wrkvecrehi,wrkvecrelo,wrkvecimhi,wrkvecimlo);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "QHb[" << i << "] : "
                    << wrkvecrehi[i] << "  " << wrkvecrehi[i] << endl << "  "
                    << wrkvecimhi[i] << "  " << wrkvecimlo[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xrehi[j] << "  " << xrelo[j] << endl << "  "
                    << ximhi[j] << "  " << ximlo[j] << endl;
         }
      }
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;
}

void CPU_dbl2_qrbs_solve
 ( int dim, int degp1, int tailidx, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_dbl2_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_head ..." << endl;

      CPU_dbl2_qrbs_head
         (dim,degp1,mathi,matlo,rhshi,rhslo,solhi,sollo,
          Qhi,Qlo,Rhi,Rlo,wrkvechi,wrkveclo,zeroQ,noqr,vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_tail ..." << endl;

      CPU_dbl2_qrbs_tail
         (dim,degp1,tailidx,mathi,matlo,rhshi,rhslo,solhi,sollo,
          Qhi,Qlo,Rhi,Rlo,wrkvechi,wrkveclo,upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_cmplx2_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo, 
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_cmplx2_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx2_qrbs_head ..." << endl;

      CPU_cmplx2_qrbs_head
         (dim,degp1,matrehi,matrelo,matimhi,matimlo,
          rhsrehi,rhsrelo,rhsimhi,rhsimlo,solrehi,solrelo,solimhi,solimlo,
          Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo,
          wrkvecrehi,wrkvecrelo,wrkvecimhi,wrkvecimlo,zeroQ,noqr,vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx2_qrbs_tail ..." << endl;

      CPU_cmplx2_qrbs_tail
         (dim,degp1,tailidx,matrehi,matrelo,matimhi,matimlo,
          rhsrehi,rhsrelo,rhsimhi,rhsimlo,solrehi,solrelo,solimhi,solimlo,
          Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo,
          wrkvecrehi,wrkvecrelo,wrkvecimhi,wrkvecimlo,
          upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_dbl2_linear_residue
 ( int dim, int degp1, int tailidx, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   int vrblvl )
{
   *resmaxhi = 0.0;
   *resmaxlo = 0.0;
   double acchi,acclo;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rihi = resvechi[i];
      double *rilo = resveclo[i];
      for(int j=0; j<dim; j++)
      {
         rihi[j] = rhshi[i][j];
         rilo[j] = rhslo[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajhi = mathi[j];
         double **Ajlo = matlo[j];
         double *xhi = solhi[i-j];
         double *xlo = sollo[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               ddf_mul(Ajhi[L][k],Ajlo[L][k],xhi[k],xlo[k],&acchi,&acclo);
               ddf_dec(&rihi[L],&rilo[L],acchi,acclo);
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solhi[i][j] << "  " << sollo[i][j] << endl;
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rihi[j] << "  " << rilo[j] << endl;
      }
      for(int j=0; j<dim; j++)
         if(abs(rihi[j]) > *resmaxhi)
         {
            *resmaxhi = abs(rihi[j]);
            *resmaxlo = abs(rilo[j]);
         }
   }
}

void CPU_cmplx2_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo, 
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo, int vrblvl )
{
   *resmaxhi = 0.0;
   *resmaxlo = 0.0;
   double acchi,acclo;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rirehi = resvecrehi[i];
      double *rirelo = resvecrelo[i];
      double *riimhi = resvecimhi[i];
      double *riimlo = resvecimlo[i];

      for(int j=0; j<dim; j++)
      {
         rirehi[j] = rhsrehi[i][j]; rirelo[j] = rhsrelo[i][j];
         riimhi[j] = rhsimhi[i][j]; riimlo[j] = rhsimlo[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajrehi = matrehi[j];
         double **Ajrelo = matrelo[j];
         double **Ajimhi = matimhi[j];
         double **Ajimlo = matimlo[j];
         double *xrehi = solrehi[i-j];
         double *xrelo = solrelo[i-j];
         double *ximhi = solimhi[i-j];
         double *ximlo = solimlo[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               // zre = Ajre[L][k]*xre[k] - Ajim[L][k]*xim[k];
               // rire[L] = rire[L] - zre;
               ddf_mul(Ajrehi[L][k],Ajrelo[L][k],xrehi[k],xrelo[k],
                       &acchi,&acclo);
               ddf_dec(&rirehi[L],&rirelo[L],acchi,acclo);
               ddf_mul(Ajimhi[L][k],Ajimlo[L][k],ximhi[k],ximlo[k],
                       &acchi,&acclo);
               ddf_inc(&rirehi[L],&rirelo[L],acchi,acclo);
               // zim = Ajre[L][k]*xim[k] + Ajim[L][k]*xre[k];
               // riim[L] = riim[L] - zim;
               ddf_mul(Ajrehi[L][k],Ajrelo[L][k],ximhi[k],ximlo[k],
                       &acchi,&acclo);
               ddf_dec(&riimhi[L],&riimlo[L],acchi,acclo);
               ddf_mul(Ajimhi[L][k],Ajimlo[L][k],xrehi[k],xrelo[k],
                       &acchi,&acclo);
               ddf_dec(&riimhi[L],&riimlo[L],acchi,acclo);
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solrehi[i][j] << "  " << solrelo[i][j] << endl << "  "
                 << solimhi[i][j] << "  " << solimlo[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rirehi[j] << "  " << rirelo[j] << endl << "  "
                 << riimhi[j] << "  " << riimlo[j] << endl;
      }
      for(int j=0; j<dim; j++)
         if(abs(rirehi[j]) + abs(riimhi[j]) > *resmaxhi)
         {
            *resmaxhi = abs(rirehi[j]) + abs(riimhi[j]);
            *resmaxlo = abs(rirelo[j]) + abs(riimlo[j]);
         }
   }
}
