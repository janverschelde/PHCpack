// The file dbl_bals_host.cpp defines functions with prototypes in
// the file dbl_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "dbl_factorizations.h"
#include "dbl_onenorms_host.h"
#include "dbl_bals_host.h"

using namespace std;

void CPU_dbl_qrbs_head
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec,
   bool *zeroQ, bool *noqr, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++) cout << " " << mat[0][i][j];
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++) cout << rhs[0][i] << endl;
   }
   double nrm;
   CPU_dbl_onenorm(dim,rhs[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm + 1.0 == 1.0)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping of CPU_dbl_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_dbl_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl_factors_houseqr ..." << endl;

      CPU_dbl_factors_houseqr(dim,dim,mat[0],Q,R); *zeroQ = false;
      CPU_dbl_factors_qrbs(dim,dim,Q,R,rhs[0],sol[0],wrkvec);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_dbl_onenorm(dim,wrkvec,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_dbl_onenorm(dim,sol[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++) cout << sol[0][i] << endl;
         for(int i=0; i<dim; i++)
         {
            wrkvec[i] = rhs[0][i];
            for(int j=0; j<dim; j++)
               wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++) cout << wrkvec[i] << endl;
      }
   }
}

void CPU_cmplx_qrbs_head
 ( int dim, int degp1, double ***matre, double ***matim,
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *wrkvecre, double *wrkvecim,
   bool *zeroQ, bool *noqr, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            cout << "  " << matre[0][i][j] << "  " << matim[0][i][j];
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhsre[0][i] << "  " << rhsim[0][i] << endl;
   }
   double nrm;
   CPU_cmplx_onenorm(dim,rhsre[0],rhsim[0],&nrm);

   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm + 1.0 == 1.0)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping of CPU_cmplx_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_cmplx_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx_factors_houseqr ..." << endl;

      CPU_cmplx_factors_houseqr(dim,dim,matre[0],matim[0],Qre,Qim,Rre,Rim);
      *zeroQ = false;
      CPU_cmplx_factors_qrbs
         (dim,dim,Qre,Qim,Rre,Rim,rhsre[0],rhsim[0],solre[0],solim[0],
          wrkvecre,wrkvecim);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_cmplx_onenorm(dim,wrkvecre,wrkvecim,&nrm);
         cout << "1-norm of Q^H*b : " << nrm << endl;
         CPU_cmplx_onenorm(dim,solre[0],solim[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double zre,zim;

         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solre[0][i] << "  " << solim[0][i] << endl;

         for(int i=0; i<dim; i++)
         {
            wrkvecre[i] = rhsre[0][i];
            wrkvecim[i] = rhsim[0][i];

            for(int j=0; j<dim; j++)
            {
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
               zre = matre[0][i][j]*solre[0][j] - matim[0][i][j]*solim[0][j];
               zim = matre[0][i][j]*solim[0][j] + matim[0][i][j]*solre[0][j];
               wrkvecre[i] = wrkvecre[i] - zre;
               wrkvecim[i] = wrkvecim[i] - zim;
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
            cout << wrkvecre[i] << "  " << wrkvecim[i] << endl;
      }
   }
}

void CPU_dbl_qrbs_tail
 ( int dim, int degp1, int tailidx, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec, int *upidx, int *bsidx,
   int *newtail, int vrblvl )
{
   double nrm;
   int skipupcnt = 0; // counts updates skipped
   int skipbscnt = 0; // counts backsubstitutions skipped
   double prevnorm = 1.0e+99; // previous norm
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      double *xi = sol[i-1];       // solution to do the update with
      CPU_dbl_onenorm(dim,xi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-15)
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
            double **Aj = mat[j-i+1]; // always start with A[1]
            double *wj = rhs[j];      // current right hand side vector

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) wj[k] = wj[k] - Aj[k][L]*xi[L];
         }
      }
      // compute sol[i] with back substitution
      double *x = sol[i];
      double *b = rhs[i];

      CPU_dbl_onenorm(dim,b,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-15) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++) x[j] = 0.0;
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
         CPU_dbl_factors_qrbs(dim,dim,Q,R,b,x,wrkvec);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "Qtb[" << i << "] : " << wrkvec[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++) cout << x[j] << endl;
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

void CPU_cmplx_qrbs_tail
 ( int dim, int degp1, int tailidx, double ***matre, double ***matim,
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *wrkvecre, double *wrkvecim, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   double nrm,zre,zim;
   int skipupcnt = 0; // counts skipped updates
   int skipbscnt = 0; // counts skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1

      double *xire = solre[i-1];    // solution to do the update with
      double *xiim = solim[i-1]; 

      CPU_cmplx_onenorm(dim,xire,xiim,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-15)
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
            double **Ajre = matre[j-i+1]; // always start with A[1]
            double **Ajim = matim[j-i+1];
            double *wjre = rhsre[j];      // current right hand side vector
            double *wjim = rhsim[j]; 

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  zre = Ajre[k][L]*xire[L] - Ajim[k][L]*xiim[L];
                  zim = Ajre[k][L]*xiim[L] + Ajim[k][L]*xire[L];
                  wjre[k] = wjre[k] - zre;
                  wjim[k] = wjim[k] - zim;
               }
         }
      }
      // compute sol[i] with back substitution
      double *xre = solre[i];
      double *xim = solim[i];
      double *bre = rhsre[i];
      double *bim = rhsim[i];

      CPU_cmplx_onenorm(dim,bre,bim,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-15) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..."
                 << endl;

         for(int j=0; j<dim; j++)
         {
            xre[j] = 0.0;
            xim[j] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;;

         if(vrblvl > 0)
            cout << "-> run backsubstitution for x[" << i << "] ..." << endl;

         if(firstbs)
         {
            *newtail = i;
            firstbs = false;
         }
         CPU_cmplx_factors_qrbs
            (dim,dim,Qre,Qim,Rre,Rim,bre,bim,xre,xim,wrkvecre,wrkvecim);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "QHb[" << i << "] : "
                    << wrkvecre[i] << "  " << wrkvecim[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xre[j] << "  " << xim[j] << endl;
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

void CPU_dbl_qrbs_solve
 ( int dim, int degp1, int tailidx, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_dbl_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_dbl_qrbs_head ..." << endl;

      CPU_dbl_qrbs_head
         (dim,degp1,mat,rhs,sol,Q,R,wrkvec,zeroQ,noqr,vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl_qrbs_tail ..." << endl;

      CPU_dbl_qrbs_tail
         (dim,degp1,tailidx,mat,rhs,sol,Q,R,wrkvec,upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_cmplx_qrbs_solve
 ( int dim, int degp1, int tailidx, double ***matre, double ***matim, 
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *wrkvecre, double *wrkvecim,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_cmplx_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx_qrbs_head ..." << endl;

      CPU_cmplx_qrbs_head
         (dim,degp1,matre,matim,rhsre,rhsim,solre,solim,
          Qre,Qim,Rre,Rim,wrkvecre,wrkvecim,zeroQ,noqr,vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx_qrbs_tail ..." << endl;

      CPU_cmplx_qrbs_tail
         (dim,degp1,tailidx,matre,matim,rhsre,rhsim,solre,solim,
          Qre,Qim,Rre,Rim,wrkvecre,wrkvecim,upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_dbl_linear_residue
 ( int dim, int degp1, int tailidx, double ***mat, double **rhs, double **sol,
   double **resvec, double *resmax, int vrblvl )
{
   *resmax = 0.0;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *ri = resvec[i];
      for(int j=0; j<dim; j++) ri[j] = rhs[i][j];
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Aj = mat[j];
         double *x = sol[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) ri[L] = ri[L] - Aj[L][k]*x[k];
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << sol[i][j] << endl;
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << ri[j] << endl;
      }
      for(int j=0; j<dim; j++)
         if(abs(ri[j]) > *resmax) *resmax = abs(ri[j]);
   }
}

void CPU_cmplx_linear_residue
 ( int dim, int degp1, int tailidx, double ***matre, double ***matim,
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **resvecre, double **resvecim, double *resmax, int vrblvl )
{
   *resmax = 0.0;
   double zre,zim;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rire = resvecre[i];
      double *riim = resvecim[i];

      for(int j=0; j<dim; j++)
      {
         rire[j] = rhsre[i][j];
         riim[j] = rhsim[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajre = matre[j];
         double **Ajim = matim[j];
         double *xre = solre[i-j];
         double *xim = solim[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               zre = Ajre[L][k]*xre[k] - Ajim[L][k]*xim[k];
               zim = Ajre[L][k]*xim[k] + Ajim[L][k]*xre[k];
               rire[L] = rire[L] - zre;
               riim[L] = riim[L] - zim;
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solre[i][j] << "  " << solim[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rire[j] << "  " << riim[j] << endl;
      }
      for(int j=0; j<dim; j++)
         if(abs(rire[j]) + abs(riim[j]) > *resmax)
            *resmax = abs(rire[j]) + abs(riim[j]);
   }
}
