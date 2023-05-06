// The file dbl4_bals_host.cpp defines functions with prototypes in
// the file dbl4_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"
#include "dbl_onenorms_host.h"
#include "dbl4_bals_host.h"

using namespace std;

void CPU_dbl4_qrbs_head
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo, 
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo,
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
            cout << "  " << mathihi[i][j]
                 << "  " << matlohi[i][j]
                 << "  " << mathilo[i][j]
                 << "  " << matlolo[i][j];
         }
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhshihi[0][i] << "  " << rhslohi[0][i] << "  "
              << rhshilo[0][i] << "  " << rhslolo[0][i] << endl;
   }
   double nrm;
   CPU_dbl_onenorm(dim,rhshihi[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm < 1.0e-56)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping CPU_dbl4_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_dbl4_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl4_factors_houseqr ..." << endl;

      CPU_dbl4_factors_houseqr
         (dim,dim,mathihi[0],matlohi[0],mathilo[0],matlolo[0],
          Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo);

      *zeroQ = false;

      CPU_dbl4_factors_qrbs
         (dim,dim,Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
          rhshihi[0],rhslohi[0],rhshilo[0],rhslolo[0],
          solhihi[0],sollohi[0],solhilo[0],sollolo[0],
          wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_dbl_onenorm(dim,wrkvechihi,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_dbl_onenorm(dim,solhihi[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double acchihi,acclohi,acchilo,acclolo;
   
         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solhihi[0][i] << "  " << sollohi[0][i] << "  "
                 << solhilo[0][i] << "  " << sollolo[0][i] << endl;

         for(int i=0; i<dim; i++)
         {
            wrkvechihi[i] = rhshihi[0][i]; wrkveclohi[i] = rhslohi[0][i];
            wrkvechilo[i] = rhshilo[0][i]; wrkveclolo[i] = rhslolo[0][i];
   
            for(int j=0; j<dim; j++)
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
            {
               qdf_mul(mathihi[0][i][j],matlohi[0][i][j],
                       mathilo[0][i][j],matlolo[0][i][j],
                       solhihi[0][j],sollohi[0][j],
                       solhilo[0][j],sollolo[0][j],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&wrkvechihi[i],&wrkveclohi[i],
                       &wrkvechilo[i],&wrkveclolo[i],
                       acchihi,acclohi,acchilo,acclolo);
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
           cout << wrkvechihi[i] << "  " << wrkveclohi[i] << "  "
                << wrkvechilo[i] << "  " << wrkveclolo[i] << endl;
      }
   }
}

void CPU_cmplx4_qrbs_head
 ( int dim, int degp1,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo,
   bool *zeroQ, bool *noqr, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix : " << endl;
      // cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            cout << "  " << matrehihi[0][i][j]
                 << "  " << matrelohi[0][i][j] << endl
                 << "  " << matrehilo[0][i][j]
                 << "  " << matrelolo[0][i][j] << endl
                 << "  " << matimhihi[0][i][j]
                 << "  " << matimlohi[0][i][j] << endl
                 << "  " << matimhilo[0][i][j]
                 << "  " << matimlolo[0][i][j] << endl;
         cout << endl;
      }
      // cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhsrehihi[0][i] << "  " << rhsrelohi[0][i] << endl << "  " 
              << rhsrehilo[0][i] << "  " << rhsrelolo[0][i] << endl << "  " 
              << rhsimhihi[0][i] << "  " << rhsimlohi[0][i] << endl << "  "
              << rhsimhilo[0][i] << "  " << rhsimlolo[0][i] << endl;
   }
   double nrm;
   CPU_cmplx_onenorm(dim,rhsrehihi[0],rhsimhihi[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm < 1.0e-56)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping CPU_cmplx4_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_cmplx4_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx4_factors_houseqr ..." << endl;

      CPU_cmplx4_factors_houseqr
         (dim,dim,matrehihi[0],matrelohi[0],matrehilo[0],matrelolo[0],
                  matimhihi[0],matimlohi[0],matimhilo[0],matimlolo[0],
          Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
          Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo);

      *zeroQ = false;

      CPU_cmplx4_factors_qrbs
         (dim,dim,Qrehihi,Qrelohi,Qrehilo,Qrelolo,
                  Qimhihi,Qimlohi,Qimhilo,Qimlolo,
          Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
          rhsrehihi[0],rhsrelohi[0],rhsrehilo[0],rhsrelolo[0],
          rhsimhihi[0],rhsimlohi[0],rhsimhilo[0],rhsimlolo[0],
          solrehihi[0],solrelohi[0],solrehilo[0],solrelolo[0],
          solimhihi[0],solimlohi[0],solimhilo[0],solimlolo[0],
          wrkvecrehihi,wrkvecrelohi,wrkvecrehilo,wrkvecrelolo,
          wrkvecimhihi,wrkvecimlohi,wrkvecimhilo,wrkvecimlolo);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_cmplx_onenorm(dim,wrkvecrehihi,wrkvecimhihi,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_cmplx_onenorm(dim,solrehihi[0],solimhihi[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double acchihi,acclohi,acchilo,acclolo;

         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solrehihi[0][i] << "  " << solrelohi[0][i] << endl << "  "
                 << solrehilo[0][i] << "  " << solrelolo[0][i] << endl << "  "
                 << solimhihi[0][i] << "  " << solimlohi[0][i] << endl << "  "
                 << solimhilo[0][i] << "  " << solimlolo[0][i] << endl;

         for(int i=0; i<dim; i++)
         {
            wrkvecrehihi[i] = rhsrehihi[0][i];
            wrkvecrelohi[i] = rhsrelohi[0][i];
            wrkvecrehilo[i] = rhsrehilo[0][i];
            wrkvecrelolo[i] = rhsrelolo[0][i];
            wrkvecimhihi[i] = rhsimhihi[0][i];
            wrkvecimlohi[i] = rhsimlohi[0][i];
            wrkvecimhilo[i] = rhsimhilo[0][i];
            wrkvecimlolo[i] = rhsimlolo[0][i];

            for(int j=0; j<dim; j++)
            {
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
               // zre = matre[0][i][j]*solre[0][j] - matim[0][i][j]*solim[0][j];
               // wrkvecre[i] = wrkvecre[i] + zre;
               qdf_mul(matrehihi[0][i][j],matrelohi[0][i][j],
                       matrehilo[0][i][j],matrelolo[0][i][j],
                       solrehihi[0][j],solrelohi[0][j],
                       solrehilo[0][j],solrelolo[0][j],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_inc(&wrkvecrehihi[i],&wrkvecrelohi[i],
                       &wrkvecrehilo[i],&wrkvecrelolo[i],
                       acchihi,acclohi,acchilo,acclolo);
               qdf_mul(matimhihi[0][i][j],matimlohi[0][i][j],
                       matimhilo[0][i][j],matimlolo[0][i][j],
                       solimhihi[0][j],solimlohi[0][j],
                       solimhilo[0][j],solimlolo[0][j],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&wrkvecrehihi[i],&wrkvecrelohi[i],
                       &wrkvecrehilo[i],&wrkvecrelolo[i],
                       acchihi,acclohi,acchilo,acclolo);
               // zim = matre[0][i][j]*solim[0][j] + matim[0][i][j]*solre[0][j];
               // wrkvecim[i] = wrkvecim[i] + zim;
               qdf_mul(matrehihi[0][i][j],matrelohi[0][i][j],
                       matrehilo[0][i][j],matrelolo[0][i][j],
                       solimhihi[0][j],solimlohi[0][j],
                       solimhilo[0][j],solimlolo[0][j],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_inc(&wrkvecimhihi[i],&wrkvecimlohi[i],
                       &wrkvecimhilo[i],&wrkvecimlolo[i],
                       acchihi,acclohi,acchilo,acclolo);
               qdf_mul(matimhihi[0][i][j],matimlohi[0][i][j],
                       matimhilo[0][i][j],matimlolo[0][i][j],
                       solrehihi[0][j],solrelohi[0][j],
                       solrehilo[0][j],solrelolo[0][j],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_inc(&wrkvecimhihi[i],&wrkvecimlohi[i],
                       &wrkvecimhilo[i],&wrkvecimlolo[i],
                       acchihi,acclohi,acchilo,acclolo);
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
            cout << wrkvecrehihi[i] << "  " << wrkvecrelohi[i] << endl << "  "
                 << wrkvecrehilo[i] << "  " << wrkvecrelolo[i] << endl << "  "
                 << wrkvecimhihi[i] << "  " << wrkvecimlohi[i] << endl << "  "
                 << wrkvecimhilo[i] << "  " << wrkvecimlolo[i] << endl;
      }
   }
}

void CPU_dbl4_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo,
   int *upidx, int *bsidx, int *newtail, int vrblvl )
{
   double nrm,acchihi,acclohi,acchilo,acclolo;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      double *xihihi = solhihi[i-1];    // solution to do the update with
      double *xilohi = sollohi[i-1];
      double *xihilo = solhilo[i-1];
      double *xilolo = sollolo[i-1];

      CPU_dbl_onenorm(dim,xihihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-56)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "skip update with x[" << i-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0) cout << "updating with x[" << i-1 << "] ..." << endl;

         for(int j=i; j<degp1; j++)
         {
            double **Ajhihi = mathihi[j-i+1]; // always start with A[1]
            double **Ajlohi = matlohi[j-i+1]; 
            double **Ajhilo = mathilo[j-i+1];
            double **Ajlolo = matlolo[j-i+1]; 
            double *wjhihi = rhshihi[j];    // current right hand side vector
            double *wjlohi = rhslohi[j];
            double *wjhilo = rhshilo[j];
            double *wjlolo = rhslolo[j];

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  qdf_mul(Ajhihi[k][L],Ajlohi[k][L],Ajhilo[k][L],Ajlolo[k][L],
                          xihihi[L],xilohi[L],xihilo[L],xilolo[L],
                          &acchihi,&acclohi,&acchilo,&acclolo);
                  qdf_dec(&wjhihi[k],&wjlohi[k],&wjhilo[k],&wjlolo[k],
                          acchihi,acclohi,acchilo,acclolo);
               }
         }
      }
      // compute sol[i] with back substitution
      double *xhihi = solhihi[i];
      double *xlohi = sollohi[i];
      double *xhilo = solhilo[i];
      double *xlolo = sollolo[i];
      double *bhihi = rhshihi[i];
      double *blohi = rhslohi[i];
      double *bhilo = rhshilo[i];
      double *blolo = rhslolo[i];

      CPU_dbl_onenorm(dim,bhihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-56) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++)
         {
            xhihi[j] = 0.0; xlohi[j] = 0.0;
            xhilo[j] = 0.0; xlolo[j] = 0.0;
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
         CPU_dbl4_factors_qrbs
            (dim,dim,Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
             bhihi,blohi,bhilo,blolo,xhihi,xlohi,xhilo,xlolo,
             wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "Qtb[" << i << "] : "
                    << wrkvechihi[i] << "  " << wrkveclohi[i] << "  "
                    << wrkvechilo[i] << "  " << wrkveclolo[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xhihi[j] << "  " << xlohi[j] << "  "
                    << xhilo[j] << "  " << xlolo[j] << endl;
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

void CPU_cmplx4_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo, 
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo,
   int *upidx, int *bsidx, int *newtail, int vrblvl )
{
   double nrm,acchihi,acclohi,acchilo,acclolo;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1

      double *xirehihi = solrehihi[i-1]; // solution to do the update with
      double *xirelohi = solrelohi[i-1]; 
      double *xirehilo = solrehilo[i-1];
      double *xirelolo = solrelolo[i-1]; 
      double *xiimhihi = solimhihi[i-1]; 
      double *xiimlohi = solimlohi[i-1]; 
      double *xiimhilo = solimhilo[i-1]; 
      double *xiimlolo = solimlolo[i-1]; 

      CPU_cmplx_onenorm(dim,xirehihi,xiimhihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-56)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "skip update with x[" << i-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0) cout << "updating with x[" << i-1 << "] ..." << endl;

         for(int j=i; j<degp1; j++)
         {
            double **Ajrehihi = matrehihi[j-i+1]; // always start with A[1]
            double **Ajrelohi = matrelohi[j-i+1];
            double **Ajrehilo = matrehilo[j-i+1];
            double **Ajrelolo = matrelolo[j-i+1];
            double **Ajimhihi = matimhihi[j-i+1];
            double **Ajimlohi = matimlohi[j-i+1];
            double **Ajimhilo = matimhilo[j-i+1];
            double **Ajimlolo = matimlolo[j-i+1];
            double *wjrehihi = rhsrehihi[j]; // current right hand side vector
            double *wjrelohi = rhsrelohi[j];
            double *wjrehilo = rhsrehilo[j];
            double *wjrelolo = rhsrelolo[j];
            double *wjimhihi = rhsimhihi[j]; 
            double *wjimlohi = rhsimlohi[j]; 
            double *wjimhilo = rhsimhilo[j]; 
            double *wjimlolo = rhsimlolo[j]; 

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  // zre = Ajre[k][L]*xire[L] - Ajim[k][L]*xiim[L];
                  // wjre[k] = wjre[k] - zre;
                  qdf_mul(Ajrehihi[k][L],Ajrelohi[k][L],
                          Ajrehilo[k][L],Ajrelolo[k][L],
                          xirehihi[L],xirelohi[L], xirehilo[L],xirelolo[L],
                          &acchihi,&acclohi,&acchilo,&acclolo);
                  qdf_dec(&wjrehihi[k],&wjrelohi[k],
                          &wjrehilo[k],&wjrelolo[k],
                          acchihi,acclohi,acchilo,acclolo);
                  qdf_mul(Ajimhihi[k][L],Ajimlohi[k][L],
                          Ajimhilo[k][L],Ajimlolo[k][L],
                          xiimhihi[L],xiimlohi[L],xiimhilo[L],xiimlolo[L],
                          &acchihi,&acclohi,&acchilo,&acclolo);
                  qdf_inc(&wjrehihi[k],&wjrelohi[k],
                          &wjrehilo[k],&wjrelolo[k],
                          acchihi,acclohi,acchilo,acclolo);
                  // zim = Ajre[k][L]*xiim[L] + Ajim[k][L]*xire[L];
                  // wjim[k] = wjim[k] - zim;
                  qdf_mul(Ajrehihi[k][L],Ajrelohi[k][L],
                          Ajrehilo[k][L],Ajrelolo[k][L],
                          xiimhihi[L],xiimlohi[L],xiimhilo[L],xiimlolo[L],
                          &acchihi,&acclohi,&acchilo,&acclolo);
                  qdf_dec(&wjimhihi[k],&wjimlohi[k],
                          &wjimhilo[k],&wjimlolo[k],
                          acchihi,acclohi,acchilo,acclolo);
                  qdf_mul(Ajimhihi[k][L],Ajimlohi[k][L],
                          Ajimhilo[k][L],Ajimlolo[k][L],
                          xirehihi[L],xirelohi[L],xirehilo[L],xirelolo[L],
                          &acchihi,&acclohi,&acchilo,&acclolo);
                  qdf_dec(&wjimhihi[k],&wjimlohi[k],
                          &wjimhilo[k],&wjimlolo[k],
                          acchihi,acclohi,acchilo,acclolo);
               }
         }
      }
      // compute sol[i] with back substitution
      double *xrehihi = solrehihi[i];
      double *xrelohi = solrelohi[i];
      double *xrehilo = solrehilo[i];
      double *xrelolo = solrelolo[i];
      double *ximhihi = solimhihi[i];
      double *ximlohi = solimlohi[i];
      double *ximhilo = solimhilo[i];
      double *ximlolo = solimlolo[i];
      double *brehihi = rhsrehihi[i];
      double *brelohi = rhsrelohi[i];
      double *brehilo = rhsrehilo[i];
      double *brelolo = rhsrelolo[i];
      double *bimhihi = rhsimhihi[i];
      double *bimlohi = rhsimlohi[i];
      double *bimhilo = rhsimhilo[i];
      double *bimlolo = rhsimlolo[i];

      CPU_cmplx_onenorm(dim,brehihi,bimhihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-56) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++)
         {
            xrehihi[j] = 0.0; xrelohi[j] = 0.0;
            xrehilo[j] = 0.0; xrelolo[j] = 0.0;
            ximhihi[j] = 0.0; ximlohi[j] = 0.0;
            ximhilo[j] = 0.0; ximlolo[j] = 0.0;
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
         CPU_cmplx4_factors_qrbs
            (dim,dim,Qrehihi,Qrelohi,Qrehilo,Qrelolo,
                     Qimhihi,Qimlohi,Qimhilo,Qimlolo,
             Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
             brehihi,brelohi,brehilo,brelolo,bimhihi,bimlohi,bimhilo,bimlolo,
             xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
             wrkvecrehihi,wrkvecrelohi,wrkvecrehilo,wrkvecrelolo,
             wrkvecimhihi,wrkvecimlohi,wrkvecimhilo,wrkvecimlolo);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "QHb[" << i << "] : "
                    << wrkvecrehihi[i] << "  " << wrkvecrehihi[i] << endl
                    << "  "
                    << wrkvecrehilo[i] << "  " << wrkvecrehilo[i] << endl
                    << "  "
                    << wrkvecimhihi[i] << "  " << wrkvecimlohi[i] << endl
                    << "  "
                    << wrkvecimhilo[i] << "  " << wrkvecimlolo[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xrehihi[j] << "  " << xrelohi[j] << endl << "  "
                    << xrehilo[j] << "  " << xrelolo[j] << endl << "  "
                    << ximhihi[j] << "  " << ximlohi[j] << endl << "  "
                    << ximhilo[j] << "  " << ximlolo[j] << endl;
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

void CPU_dbl4_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_dbl4_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_dbl4_qrbs_head ..." << endl;

      CPU_dbl4_qrbs_head
         (dim,degp1,mathihi,matlohi,mathilo,matlolo,
          rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
          Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
          wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo,zeroQ,noqr,vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl4_qrbs_tail ..." << endl;

      CPU_dbl4_qrbs_tail
         (dim,degp1,tailidx,mathihi,matlohi,mathilo,matlolo,
          rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
          Qhihi,Qlohi,Qhilo,Qlolo,Rhihi,Rlohi,Rhilo,Rlolo,
          wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo,
          upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_cmplx4_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi, 
   double ***matimhilo, double ***matimlolo, 
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_cmplx4_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx4_qrbs_head ..." << endl;

      CPU_cmplx4_qrbs_head
         (dim,degp1,matrehihi,matrelohi,matrehilo,matrelolo,
                    matimhihi,matimlohi,matimhilo,matimlolo,
          rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
          rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
          solrehihi,solrelohi,solrehilo,solrelolo,
          solimhihi,solimlohi,solimhilo,solimlolo,
          Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
          Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
          wrkvecrehihi,wrkvecrelohi,wrkvecrehilo,wrkvecrelolo,
          wrkvecimhihi,wrkvecimlohi,wrkvecimhilo,wrkvecimlolo,zeroQ,noqr,
          vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx4_qrbs_tail ..." << endl;

      CPU_cmplx4_qrbs_tail
         (dim,degp1,tailidx,
          matrehihi,matrelohi,matrehilo,matrelolo,
          matimhihi,matimlohi,matimhilo,matimlolo,
          rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
          rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
          solrehihi,solrelohi,solrehilo,solrelolo,
          solimhihi,solimlohi,solimhilo,solimlolo,
          Qrehihi,Qrelohi,Qrehilo,Qrelolo,Qimhihi,Qimlohi,Qimhilo,Qimlolo,
          Rrehihi,Rrelohi,Rrehilo,Rrelolo,Rimhihi,Rimlohi,Rimhilo,Rimlolo,
          wrkvecrehihi,wrkvecrelohi,wrkvecrehilo,wrkvecrelolo,
          wrkvecimhihi,wrkvecimlohi,wrkvecimhilo,wrkvecimlolo,
          upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_dbl4_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl )
{
   *resmaxhihi = 0.0;
   *resmaxlohi = 0.0;
   *resmaxhilo = 0.0;
   *resmaxlolo = 0.0;
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rihihi = resvechihi[i];
      double *rilohi = resveclohi[i];
      double *rihilo = resvechilo[i];
      double *rilolo = resveclolo[i];

      for(int j=0; j<dim; j++)
      {
         rihihi[j] = rhshihi[i][j];
         rilohi[j] = rhslohi[i][j];
         rihilo[j] = rhshilo[i][j];
         rilolo[j] = rhslolo[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajhihi = mathihi[j];
         double **Ajlohi = matlohi[j];
         double **Ajhilo = mathilo[j];
         double **Ajlolo = matlolo[j];
         double *xhihi = solhihi[i-j];
         double *xlohi = sollohi[i-j];
         double *xhilo = solhilo[i-j];
         double *xlolo = sollolo[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               qdf_mul(Ajhihi[L][k],Ajlohi[L][k],Ajhilo[L][k],Ajlolo[L][k],
                       xhihi[k],xlohi[k],xhilo[k],xlolo[k],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&rihihi[L],&rilohi[L],&rihilo[L],&rilolo[L],
                       acchihi,acclohi,acchilo,acclolo);
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << solhihi[i][j] << "  " << sollohi[i][j] << endl;
            cout << solhilo[i][j] << "  " << sollolo[i][j] << endl;
         }
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << rihihi[j] << "  " << rilohi[j] << endl;
            cout << rihilo[j] << "  " << rilolo[j] << endl;
         }
      }
      for(int j=0; j<dim; j++)
         if(abs(rihihi[j]) > *resmaxhihi)
         {
            *resmaxhihi = abs(rihihi[j]);
            *resmaxlohi = abs(rilohi[j]);
            *resmaxhilo = abs(rihilo[j]);
            *resmaxlolo = abs(rilolo[j]);
         }
   }
}

void CPU_cmplx4_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi, 
   double **rhsimhilo, double **rhsimlolo, 
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl )
{
   *resmaxhihi = 0.0;
   *resmaxlohi = 0.0;
   *resmaxhilo = 0.0;
   *resmaxlolo = 0.0;
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rirehihi = resvecrehihi[i];
      double *rirelohi = resvecrelohi[i];
      double *rirehilo = resvecrehilo[i];
      double *rirelolo = resvecrelolo[i];
      double *riimhihi = resvecimhihi[i];
      double *riimlohi = resvecimlohi[i];
      double *riimhilo = resvecimhilo[i];
      double *riimlolo = resvecimlolo[i];

      for(int j=0; j<dim; j++)
      {
         rirehihi[j] = rhsrehihi[i][j]; rirelohi[j] = rhsrelohi[i][j];
         rirehilo[j] = rhsrehilo[i][j]; rirelolo[j] = rhsrelolo[i][j];
         riimhihi[j] = rhsimhihi[i][j]; riimlohi[j] = rhsimlohi[i][j];
         riimhilo[j] = rhsimhilo[i][j]; riimlolo[j] = rhsimlolo[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajrehihi = matrehihi[j];
         double **Ajrelohi = matrelohi[j];
         double **Ajrehilo = matrehilo[j];
         double **Ajrelolo = matrelolo[j];
         double **Ajimhihi = matimhihi[j];
         double **Ajimlohi = matimlohi[j];
         double **Ajimhilo = matimhilo[j];
         double **Ajimlolo = matimlolo[j];
         double *xrehihi = solrehihi[i-j];
         double *xrelohi = solrelohi[i-j];
         double *xrehilo = solrehilo[i-j];
         double *xrelolo = solrelolo[i-j];
         double *ximhihi = solimhihi[i-j];
         double *ximlohi = solimlohi[i-j];
         double *ximhilo = solimhilo[i-j];
         double *ximlolo = solimlolo[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               // zre = Ajre[L][k]*xre[k] - Ajim[L][k]*xim[k];
               // rire[L] = rire[L] - zre;
               qdf_mul(Ajrehihi[L][k],Ajrelohi[L][k],
                       Ajrehilo[L][k],Ajrelolo[L][k],
                       xrehihi[k],xrelohi[k],xrehilo[k],xrelolo[k],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&rirehihi[L],&rirelohi[L],
                       &rirehilo[L],&rirelolo[L],
                       acchihi,acclohi,acchilo,acclolo);
               qdf_mul(Ajimhihi[L][k],Ajimlohi[L][k],
                       Ajimhilo[L][k],Ajimlolo[L][k],
                       ximhihi[k],ximlohi[k],ximhilo[k],ximlolo[k],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_inc(&rirehihi[L],&rirelohi[L],
                       &rirehilo[L],&rirelolo[L],
                       acchihi,acclohi,acchilo,acclolo);
               // zim = Ajre[L][k]*xim[k] + Ajim[L][k]*xre[k];
               // riim[L] = riim[L] - zim;
               qdf_mul(Ajrehihi[L][k],Ajrelohi[L][k],
                       Ajrehilo[L][k],Ajrelolo[L][k],
                       ximhihi[k],ximlohi[k],ximhilo[k],ximlolo[k],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&riimhihi[L],&riimlohi[L],
                       &riimhilo[L],&riimlolo[L],
                       acchihi,acclohi,acchilo,acclolo);
               qdf_mul(Ajimhihi[L][k],Ajimlohi[L][k],
                       Ajimhilo[L][k],Ajimlolo[L][k],
                       xrehihi[k],xrelohi[k],xrehilo[k],xrelolo[k],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&riimhihi[L],&riimlohi[L],
                       &riimhilo[L],&riimlolo[L],
                       acchihi,acclohi,acchilo,acclolo);
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solrehihi[i][j] << "  " << solrelohi[i][j] << endl << "  "
                 << solrehilo[i][j] << "  " << solrelolo[i][j] << endl << "  "
                 << solimhihi[i][j] << "  " << solimlohi[i][j] << endl << "  "
                 << solimhilo[i][j] << "  " << solimlolo[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rirehihi[j] << "  " << rirelohi[j] << endl << "  "
                 << rirehilo[j] << "  " << rirelolo[j] << endl << "  "
                 << riimhihi[j] << "  " << riimlohi[j] << endl << "  "
                 << riimhilo[j] << "  " << riimlolo[j] << endl;
      }
      for(int j=0; j<dim; j++)
         if(abs(rirehihi[j]) + abs(riimhihi[j]) > *resmaxhihi)
         {
            *resmaxhihi = abs(rirehihi[j]) + abs(riimhihi[j]);
            *resmaxlohi = abs(rirelohi[j]) + abs(riimlohi[j]);
            *resmaxhilo = abs(rirehilo[j]) + abs(riimhilo[j]);
            *resmaxlolo = abs(rirelolo[j]) + abs(riimlolo[j]);
         }
   }
}
