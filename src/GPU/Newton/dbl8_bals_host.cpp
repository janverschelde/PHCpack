// The file dbl8_bals_host.cpp defines functions with prototypes in
// the file dbl8_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"
#include "dbl_onenorms_host.h"
#include "dbl8_bals_host.h"

using namespace std;

void CPU_dbl8_qrbs_head
 ( int dim, int degp1,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi, 
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo, 
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo,
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
            cout << "  " << mathihihi[0][i][j]
                 << "  " << matlohihi[0][i][j]
                 << "  " << mathilohi[0][i][j]
                 << "  " << matlolohi[0][i][j]
                 << "  " << mathihilo[0][i][j]
                 << "  " << matlohilo[0][i][j]
                 << "  " << mathilolo[0][i][j]
                 << "  " << matlololo[0][i][j];
         }
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhshihihi[0][i] << "  " << rhslohihi[0][i] << "  "
              << rhshilohi[0][i] << "  " << rhslolohi[0][i] << endl << "  "
              << rhshihilo[0][i] << "  " << rhslohilo[0][i] << "  "
              << rhshilolo[0][i] << "  " << rhslololo[0][i] << endl;
   }
   double nrm;
   CPU_dbl_onenorm(dim,rhshihihi[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm < 1.0e-120)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping CPU_dbl8_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_dbl8_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl8_factors_houseqr ..." << endl;

      CPU_dbl8_factors_houseqr
         (dim,dim,
          mathihihi[0],matlohihi[0],mathilohi[0],matlolohi[0],
          mathihilo[0],matlohilo[0],mathilolo[0],matlololo[0],
          Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
          Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo);

      *zeroQ = false;

      CPU_dbl8_factors_qrbs
         (dim,dim,
          Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
          Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
          rhshihihi[0],rhslohihi[0],rhshilohi[0],rhslolohi[0],
          rhshihilo[0],rhslohilo[0],rhshilolo[0],rhslololo[0],
          solhihihi[0],sollohihi[0],solhilohi[0],sollolohi[0],
          solhihilo[0],sollohilo[0],solhilolo[0],sollololo[0],
          wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
          wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_dbl_onenorm(dim,wrkvechihihi,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_dbl_onenorm(dim,solhihihi[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double acchihihi,acclohihi,acchilohi,acclolohi;
         double acchihilo,acclohilo,acchilolo,acclololo;

         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solhihihi[0][i] << "  " << sollohihi[0][i] << "  "
                 << solhilohi[0][i] << "  " << sollolohi[0][i] << endl << "  "
                 << solhihilo[0][i] << "  " << sollohilo[0][i] << "  "
                 << solhilolo[0][i] << "  " << sollololo[0][i] << endl;

         for(int i=0; i<dim; i++)
         {
            wrkvechihihi[i] = rhshihihi[0][i];
            wrkveclohihi[i] = rhslohihi[0][i];
            wrkvechilohi[i] = rhshilohi[0][i];
            wrkveclolohi[i] = rhslolohi[0][i];
            wrkvechihilo[i] = rhshihilo[0][i];
            wrkveclohilo[i] = rhslohilo[0][i];
            wrkvechilolo[i] = rhshilolo[0][i];
            wrkveclololo[i] = rhslololo[0][i];

            for(int j=0; j<dim; j++)
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
            {
               odf_mul(mathihihi[0][i][j],matlohihi[0][i][j],
                       mathilohi[0][i][j],matlolohi[0][i][j],
                       mathihilo[0][i][j],matlohilo[0][i][j],
                       mathilolo[0][i][j],matlololo[0][i][j],
                       solhihihi[0][j],sollohihi[0][j],
                       solhilohi[0][j],sollolohi[0][j],
                       solhihilo[0][j],sollohilo[0][j],
                       solhilolo[0][j],sollololo[0][j],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_dec(&wrkvechihihi[i],&wrkveclohihi[i],
                       &wrkvechilohi[i],&wrkveclolohi[i],
                       &wrkvechihilo[i],&wrkveclohilo[i],
                       &wrkvechilolo[i],&wrkveclololo[i],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
           cout << wrkvechihihi[i] << "  " << wrkveclohihi[i] << "  "
                << wrkvechilohi[i] << "  " << wrkveclolohi[i] << endl << "  "
                << wrkvechihilo[i] << "  " << wrkveclohilo[i] << "  "
                << wrkvechilolo[i] << "  " << wrkveclololo[i] << endl;
      }
   }
}

void CPU_cmplx8_qrbs_head
 ( int dim, int degp1,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo,
   bool *zeroQ, bool *noqr, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix : " << endl;
      // cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            cout << "  " << matrehihihi[0][i][j]
                 << "  " << matrelohihi[0][i][j] << endl
                 << "  " << matrehilohi[0][i][j]
                 << "  " << matrelolohi[0][i][j] << endl
                 << "  " << matrehihilo[0][i][j]
                 << "  " << matrelohilo[0][i][j] << endl
                 << "  " << matrehilolo[0][i][j]
                 << "  " << matrelololo[0][i][j] << endl
                 << "  " << matimhihihi[0][i][j]
                 << "  " << matimlohihi[0][i][j] << endl
                 << "  " << matimhilohi[0][i][j]
                 << "  " << matimlolohi[0][i][j] << endl
                 << "  " << matimhihilo[0][i][j]
                 << "  " << matimlohilo[0][i][j] << endl
                 << "  " << matimhilolo[0][i][j]
                 << "  " << matimlololo[0][i][j] << endl;
         cout << endl;
      }
      // cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhsrehihihi[0][i] << "  " << rhsrelohihi[0][i] << endl
              << "  " 
              << rhsrehilohi[0][i] << "  " << rhsrelolohi[0][i] << endl
              << "  " 
              << rhsrehihilo[0][i] << "  " << rhsrelohilo[0][i] << endl
              << "  " 
              << rhsrehilolo[0][i] << "  " << rhsrelololo[0][i] << endl
              << "  " 
              << rhsimhihihi[0][i] << "  " << rhsimlohihi[0][i] << endl
              << "  "
              << rhsimhilohi[0][i] << "  " << rhsimlolohi[0][i] << endl
              << "  "
              << rhsimhihilo[0][i] << "  " << rhsimlohilo[0][i] << endl
              << "  "
              << rhsimhilolo[0][i] << "  " << rhsimlololo[0][i] << endl;
   }
   double nrm;
   CPU_cmplx_onenorm(dim,rhsrehihihi[0],rhsimhihihi[0],&nrm);
   if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

   if(nrm < 1.0e-120)
   {
      if(*zeroQ)
      {
         if(vrblvl > 0)
            cout << "no skipping CPU_cmplx8_factors_houseqr because zeroQ"
                 << endl;

         *noqr = false;
      }
      else
      {
         if(vrblvl > 0)
            cout << "skip call to CPU_cmplx8_factors_houseqr ..." << endl;

         *noqr = true;
      }
   }
   if(!*noqr)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx8_factors_houseqr ..." << endl;

      CPU_cmplx8_factors_houseqr
         (dim,dim,
          matrehihihi[0],matrelohihi[0],matrehilohi[0],matrelolohi[0],
          matrehihilo[0],matrelohilo[0],matrehilolo[0],matrelololo[0],
          matimhihihi[0],matimlohihi[0],matimhilohi[0],matimlolohi[0],
          matimhihilo[0],matimlohilo[0],matimhilolo[0],matimlololo[0],
          Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
          Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
          Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
          Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
          Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
          Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
          Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
          Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo);

      *zeroQ = false;

      CPU_cmplx8_factors_qrbs
         (dim,dim,
          Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
          Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
          Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
          Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
          Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
          Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
          Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
          Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
          rhsrehihihi[0],rhsrelohihi[0],rhsrehilohi[0],rhsrelolohi[0],
          rhsrehihilo[0],rhsrelohilo[0],rhsrehilolo[0],rhsrelololo[0],
          rhsimhihihi[0],rhsimlohihi[0],rhsimhilohi[0],rhsimlolohi[0],
          rhsimhihilo[0],rhsimlohilo[0],rhsimhilolo[0],rhsimlololo[0],
          solrehihihi[0],solrelohihi[0],solrehilohi[0],solrelolohi[0],
          solrehihilo[0],solrelohilo[0],solrehilolo[0],solrelololo[0],
          solimhihihi[0],solimlohihi[0],solimhilohi[0],solimlolohi[0],
          solimhihilo[0],solimlohilo[0],solimhilolo[0],solimlololo[0],
          wrkvecrehihihi,wrkvecrelohihi,wrkvecrehilohi,wrkvecrelolohi,
          wrkvecrehihilo,wrkvecrelohilo,wrkvecrehilolo,wrkvecrelololo,
          wrkvecimhihihi,wrkvecimlohihi,wrkvecimhilohi,wrkvecimlolohi,
          wrkvecimhihilo,wrkvecimlohilo,wrkvecimhilolo,wrkvecimlololo);

      if(vrblvl > 0)
      {
         double nrm;

         CPU_cmplx_onenorm(dim,wrkvecrehihihi,wrkvecimhihihi,&nrm);
         cout << "1-norm of Q^T*b : " << nrm << endl;
         CPU_cmplx_onenorm(dim,solrehihihi[0],solimhihihi[0],&nrm);
         cout << "1-norm of x : " << nrm << endl;
      }
      if(vrblvl > 1)
      {
         double acchihihi,acclohihi,acchilohi,acclolohi;
         double acchihilo,acclohilo,acchilolo,acclololo;

         cout << "The leading coefficients of the solution :" << endl;
         for(int i=0; i<dim; i++)
            cout << solrehihihi[0][i] << "  " << solrelohihi[0][i] << endl
                 << "  "
                 << solrehilohi[0][i] << "  " << solrelolohi[0][i] << endl
                 << "  "
                 << solrehihilo[0][i] << "  " << solrelohilo[0][i] << endl
                 << "  "
                 << solrehilolo[0][i] << "  " << solrelololo[0][i] << endl
                 << "  "
                 << solimhihihi[0][i] << "  " << solimlohihi[0][i] << endl
                 << "  "
                 << solimhilohi[0][i] << "  " << solimlolohi[0][i] << endl
                 << "  "
                 << solimhihilo[0][i] << "  " << solimlohilo[0][i] << endl
                 << "  "
                 << solimhilolo[0][i] << "  " << solimlololo[0][i] << endl;

         for(int i=0; i<dim; i++)
         {
            wrkvecrehihihi[i] = rhsrehihihi[0][i];
            wrkvecrelohihi[i] = rhsrelohihi[0][i];
            wrkvecrehilohi[i] = rhsrehilohi[0][i];
            wrkvecrelolohi[i] = rhsrelolohi[0][i];
            wrkvecrehihilo[i] = rhsrehihilo[0][i];
            wrkvecrelohilo[i] = rhsrelohilo[0][i];
            wrkvecrehilolo[i] = rhsrehilolo[0][i];
            wrkvecrelololo[i] = rhsrelololo[0][i];
            wrkvecimhihihi[i] = rhsimhihihi[0][i];
            wrkvecimlohihi[i] = rhsimlohihi[0][i];
            wrkvecimhilohi[i] = rhsimhilohi[0][i];
            wrkvecimlolohi[i] = rhsimlolohi[0][i];
            wrkvecimhihilo[i] = rhsimhihilo[0][i];
            wrkvecimlohilo[i] = rhsimlohilo[0][i];
            wrkvecimhilolo[i] = rhsimhilolo[0][i];
            wrkvecimlololo[i] = rhsimlololo[0][i];
   
            for(int j=0; j<dim; j++)
            {
               // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
               // zre = matre[0][i][j]*solre[0][j] - matim[0][i][j]*solim[0][j];
               // wrkvecre[i] = wrkvecre[i] + zre;
               odf_mul(matrehihihi[0][i][j],matrelohihi[0][i][j],
                       matrehilohi[0][i][j],matrelolohi[0][i][j],
                       matrehihilo[0][i][j],matrelohilo[0][i][j],
                       matrehilolo[0][i][j],matrelololo[0][i][j],
                       solrehihihi[0][j],solrelohihi[0][j],
                       solrehilohi[0][j],solrelolohi[0][j],
                       solrehihilo[0][j],solrelohilo[0][j],
                       solrehilolo[0][j],solrelololo[0][j],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_inc(&wrkvecrehihihi[i],&wrkvecrelohihi[i],
                       &wrkvecrehilohi[i],&wrkvecrelolohi[i],
                       &wrkvecrehihilo[i],&wrkvecrelohilo[i],
                       &wrkvecrehilolo[i],&wrkvecrelololo[i],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
               odf_mul(matimhihihi[0][i][j],matimlohihi[0][i][j],
                       matimhilohi[0][i][j],matimlolohi[0][i][j],
                       matimhihilo[0][i][j],matimlohilo[0][i][j],
                       matimhilolo[0][i][j],matimlololo[0][i][j],
                       solimhihihi[0][j],solimlohihi[0][j],
                       solimhilohi[0][j],solimlolohi[0][j],
                       solimhihilo[0][j],solimlohilo[0][j],
                       solimhilolo[0][j],solimlololo[0][j],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_dec(&wrkvecrehihihi[i],&wrkvecrelohihi[i],
                       &wrkvecrehilohi[i],&wrkvecrelolohi[i],
                       &wrkvecrehihilo[i],&wrkvecrelohilo[i],
                       &wrkvecrehilolo[i],&wrkvecrelololo[i],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
               // zim = matre[0][i][j]*solim[0][j] + matim[0][i][j]*solre[0][j];
               // wrkvecim[i] = wrkvecim[i] + zim;
               odf_mul(matrehihihi[0][i][j],matrelohihi[0][i][j],
                       matrehilohi[0][i][j],matrelolohi[0][i][j],
                       matrehihilo[0][i][j],matrelohilo[0][i][j],
                       matrehilolo[0][i][j],matrelololo[0][i][j],
                       solimhihihi[0][j],solimlohihi[0][j],
                       solimhilohi[0][j],solimlolohi[0][j],
                       solimhihilo[0][j],solimlohilo[0][j],
                       solimhilolo[0][j],solimlololo[0][j],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_inc(&wrkvecimhihihi[i],&wrkvecimlohihi[i],
                       &wrkvecimhilohi[i],&wrkvecimlolohi[i],
                       &wrkvecimhihilo[i],&wrkvecimlohilo[i],
                       &wrkvecimhilolo[i],&wrkvecimlololo[i],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
               odf_mul(matimhihihi[0][i][j],matimlohihi[0][i][j],
                       matimhilohi[0][i][j],matimlolohi[0][i][j],
                       matimhihilo[0][i][j],matimlohilo[0][i][j],
                       matimhilolo[0][i][j],matimlololo[0][i][j],
                       solrehihihi[0][j],solrelohihi[0][j],
                       solrehilohi[0][j],solrelolohi[0][j],
                       solrehihilo[0][j],solrelohilo[0][j],
                       solrehilolo[0][j],solrelololo[0][j],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_inc(&wrkvecimhihihi[i],&wrkvecimlohihi[i],
                       &wrkvecimhilohi[i],&wrkvecimlolohi[i],
                       &wrkvecimhihilo[i],&wrkvecimlohilo[i],
                       &wrkvecimhilolo[i],&wrkvecimlololo[i],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
            }
         }
         cout << "The residual vector :" << endl;
         for(int i=0; i<dim; i++)
            cout << wrkvecrehihihi[i] << "  " << wrkvecrelohihi[i] << endl
                 << "  "
                 << wrkvecrehilohi[i] << "  " << wrkvecrelolohi[i] << endl
                 << "  "
                 << wrkvecrehihilo[i] << "  " << wrkvecrelohilo[i] << endl
                 << "  "
                 << wrkvecrehilolo[i] << "  " << wrkvecrelololo[i] << endl
                 << "  "
                 << wrkvecimhihihi[i] << "  " << wrkvecimlohihi[i] << endl
                 << "  "
                 << wrkvecimhilohi[i] << "  " << wrkvecimlolohi[i] << endl
                 << "  "
                 << wrkvecimhihilo[i] << "  " << wrkvecimlohilo[i] << endl
                 << "  "
                 << wrkvecimhilolo[i] << "  " << wrkvecimlololo[i] << endl;
      }
   }
}

void CPU_dbl8_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo,
   int *upidx, int *bsidx, int *newtail, int vrblvl )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double nrm;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1

      double *xihihihi = solhihihi[i-1]; // solution to do the update with
      double *xilohihi = sollohihi[i-1];
      double *xihilohi = solhilohi[i-1];
      double *xilolohi = sollolohi[i-1];
      double *xihihilo = solhihilo[i-1];
      double *xilohilo = sollohilo[i-1];
      double *xihilolo = solhilolo[i-1];
      double *xilololo = sollololo[i-1];

      CPU_dbl_onenorm(dim,xihihihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-120)
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
            double **Ajhihihi = mathihihi[j-i+1]; // always start with A[1]
            double **Ajlohihi = matlohihi[j-i+1]; 
            double **Ajhilohi = mathilohi[j-i+1];
            double **Ajlolohi = matlolohi[j-i+1]; 
            double **Ajhihilo = mathihilo[j-i+1];
            double **Ajlohilo = matlohilo[j-i+1]; 
            double **Ajhilolo = mathilolo[j-i+1];
            double **Ajlololo = matlololo[j-i+1]; 
            double *wjhihihi = rhshihihi[j]; // current right hand side vector
            double *wjlohihi = rhslohihi[j];
            double *wjhilohi = rhshilohi[j];
            double *wjlolohi = rhslolohi[j];
            double *wjhihilo = rhshihilo[j];
            double *wjlohilo = rhslohilo[j];
            double *wjhilolo = rhshilolo[j];
            double *wjlololo = rhslololo[j];

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  odf_mul(Ajhihihi[k][L],Ajlohihi[k][L],
                          Ajhilohi[k][L],Ajlolohi[k][L],
                          Ajhihilo[k][L],Ajlohilo[k][L],
                          Ajhilolo[k][L],Ajlololo[k][L],
                          xihihihi[L],xilohihi[L],xihilohi[L],xilolohi[L],
                          xihihilo[L],xilohilo[L],xihilolo[L],xilololo[L],
                          &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                          &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  odf_dec(&wjhihihi[k],&wjlohihi[k],&wjhilohi[k],&wjlolohi[k],
                          &wjhihilo[k],&wjlohilo[k],&wjhilolo[k],&wjlololo[k],
                          acchihihi,acclohihi,acchilohi,acclolohi,
                          acchihilo,acclohilo,acchilolo,acclololo);
               }
         }
      }
      // compute sol[i] with back substitution
      double *xhihihi = solhihihi[i];
      double *xlohihi = sollohihi[i];
      double *xhilohi = solhilohi[i];
      double *xlolohi = sollolohi[i];
      double *xhihilo = solhihilo[i];
      double *xlohilo = sollohilo[i];
      double *xhilolo = solhilolo[i];
      double *xlololo = sollololo[i];
      double *bhihihi = rhshihihi[i];
      double *blohihi = rhslohihi[i];
      double *bhilohi = rhshilohi[i];
      double *blolohi = rhslolohi[i];
      double *bhihilo = rhshihilo[i];
      double *blohilo = rhslohilo[i];
      double *bhilolo = rhshilolo[i];
      double *blololo = rhslololo[i];

      CPU_dbl_onenorm(dim,bhihihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-120) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++)
         {
            xhihihi[j] = 0.0; xlohihi[j] = 0.0;
            xhilohi[j] = 0.0; xlolohi[j] = 0.0;
            xhihilo[j] = 0.0; xlohilo[j] = 0.0;
            xhilolo[j] = 0.0; xlololo[j] = 0.0;
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
         CPU_dbl8_factors_qrbs
            (dim,dim,
             Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
             Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
             bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
             xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
             wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
             wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "Qtb[" << i << "] : "
                    << wrkvechihihi[i] << "  " << wrkveclohihi[i] << "  "
                    << wrkvechilohi[i] << "  " << wrkveclolohi[i] << endl
                    << "  "
                    << wrkvechihilo[i] << "  " << wrkveclohilo[i] << "  "
                    << wrkvechilolo[i] << "  " << wrkveclololo[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xhihihi[j] << "  " << xlohihi[j] << "  "
                    << xhilohi[j] << "  " << xlolohi[j] << endl << "  "
                    << xhihilo[j] << "  " << xlohilo[j] << "  "
                    << xhilolo[j] << "  " << xlololo[j] << endl;
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

void CPU_cmplx8_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi, 
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo, 
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo,
   int *upidx, int *bsidx, int *newtail, int vrblvl )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double nrm;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   for(int i=tailidx; i<degp1; i++)
   {
      if(vrblvl > 0) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1

      double *xirehihihi = solrehihihi[i-1]; // solution in the update
      double *xirelohihi = solrelohihi[i-1]; 
      double *xirehilohi = solrehilohi[i-1];
      double *xirelolohi = solrelolohi[i-1]; 
      double *xirehihilo = solrehihilo[i-1];
      double *xirelohilo = solrelohilo[i-1]; 
      double *xirehilolo = solrehilolo[i-1];
      double *xirelololo = solrelololo[i-1]; 
      double *xiimhihihi = solimhihihi[i-1]; 
      double *xiimlohihi = solimlohihi[i-1]; 
      double *xiimhilohi = solimhilohi[i-1]; 
      double *xiimlolohi = solimlolohi[i-1]; 
      double *xiimhihilo = solimhihilo[i-1]; 
      double *xiimlohilo = solimlohilo[i-1]; 
      double *xiimhilolo = solimhilolo[i-1]; 
      double *xiimlololo = solimlololo[i-1]; 

      CPU_cmplx_onenorm(dim,xirehihihi,xiimhihihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of x[" << i-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-120)
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
            double **Ajrehihihi = matrehihihi[j-i+1]; // always start with A[1]
            double **Ajrelohihi = matrelohihi[j-i+1];
            double **Ajrehilohi = matrehilohi[j-i+1];
            double **Ajrelolohi = matrelolohi[j-i+1];
            double **Ajrehihilo = matrehihilo[j-i+1];
            double **Ajrelohilo = matrelohilo[j-i+1];
            double **Ajrehilolo = matrehilolo[j-i+1];
            double **Ajrelololo = matrelololo[j-i+1];
            double **Ajimhihihi = matimhihihi[j-i+1];
            double **Ajimlohihi = matimlohihi[j-i+1];
            double **Ajimhilohi = matimhilohi[j-i+1];
            double **Ajimlolohi = matimlolohi[j-i+1];
            double **Ajimhihilo = matimhihilo[j-i+1];
            double **Ajimlohilo = matimlohilo[j-i+1];
            double **Ajimhilolo = matimhilolo[j-i+1];
            double **Ajimlololo = matimlololo[j-i+1];
            double *wjrehihihi = rhsrehihihi[j]; // current right hand side
            double *wjrelohihi = rhsrelohihi[j];
            double *wjrehilohi = rhsrehilohi[j];
            double *wjrelolohi = rhsrelolohi[j];
            double *wjrehihilo = rhsrehihilo[j];
            double *wjrelohilo = rhsrelohilo[j];
            double *wjrehilolo = rhsrehilolo[j];
            double *wjrelololo = rhsrelololo[j];
            double *wjimhihihi = rhsimhihihi[j]; 
            double *wjimlohihi = rhsimlohihi[j]; 
            double *wjimhilohi = rhsimhilohi[j]; 
            double *wjimlolohi = rhsimlolohi[j]; 
            double *wjimhihilo = rhsimhihilo[j]; 
            double *wjimlohilo = rhsimlohilo[j]; 
            double *wjimhilolo = rhsimhilolo[j]; 
            double *wjimlololo = rhsimlololo[j]; 

            for(int k=0; k<dim; k++)
               for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
               {
                  // zre = Ajre[k][L]*xire[L] - Ajim[k][L]*xiim[L];
                  // wjre[k] = wjre[k] - zre;
                  odf_mul(Ajrehihihi[k][L],Ajrelohihi[k][L],
                          Ajrehilohi[k][L],Ajrelolohi[k][L],
                          Ajrehihilo[k][L],Ajrelohilo[k][L],
                          Ajrehilolo[k][L],Ajrelololo[k][L],
                          xirehihihi[L],xirelohihi[L],
                          xirehilohi[L],xirelolohi[L],
                          xirehihilo[L],xirelohilo[L],
                          xirehilolo[L],xirelololo[L],
                          &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                          &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  odf_dec(&wjrehihihi[k],&wjrelohihi[k],
                          &wjrehilohi[k],&wjrelolohi[k],
                          &wjrehihilo[k],&wjrelohilo[k],
                          &wjrehilolo[k],&wjrelololo[k],
                          acchihihi,acclohihi,acchilohi,acclolohi,
                          acchihilo,acclohilo,acchilolo,acclololo);
                  odf_mul(Ajimhihihi[k][L],Ajimlohihi[k][L],
                          Ajimhilohi[k][L],Ajimlolohi[k][L],
                          Ajimhihilo[k][L],Ajimlohilo[k][L],
                          Ajimhilolo[k][L],Ajimlololo[k][L],
                          xiimhihihi[L],xiimlohihi[L],
                          xiimhilohi[L],xiimlolohi[L],
                          xiimhihilo[L],xiimlohilo[L],
                          xiimhilolo[L],xiimlololo[L],
                          &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                          &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  odf_inc(&wjrehihihi[k],&wjrelohihi[k],
                          &wjrehilohi[k],&wjrelolohi[k],
                          &wjrehihilo[k],&wjrelohilo[k],
                          &wjrehilolo[k],&wjrelololo[k],
                          acchihihi,acclohihi,acchilohi,acclolohi,
                          acchihilo,acclohilo,acchilolo,acclololo);
                  // zim = Ajre[k][L]*xiim[L] + Ajim[k][L]*xire[L];
                  // wjim[k] = wjim[k] - zim;
                  odf_mul(Ajrehihihi[k][L],Ajrelohihi[k][L],
                          Ajrehilohi[k][L],Ajrelolohi[k][L],
                          Ajrehihilo[k][L],Ajrelohilo[k][L],
                          Ajrehilolo[k][L],Ajrelololo[k][L],
                          xiimhihihi[L],xiimlohihi[L],
                          xiimhilohi[L],xiimlolohi[L],
                          xiimhihilo[L],xiimlohilo[L],
                          xiimhilolo[L],xiimlololo[L],
                          &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                          &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  odf_dec(&wjimhihihi[k],&wjimlohihi[k],
                          &wjimhilohi[k],&wjimlolohi[k],
                          &wjimhihilo[k],&wjimlohilo[k],
                          &wjimhilolo[k],&wjimlololo[k],
                          acchihihi,acclohihi,acchilohi,acclolohi,
                          acchihilo,acclohilo,acchilolo,acclololo);
                  odf_mul(Ajimhihihi[k][L],Ajimlohihi[k][L],
                          Ajimhilohi[k][L],Ajimlolohi[k][L],
                          Ajimhihilo[k][L],Ajimlohilo[k][L],
                          Ajimhilolo[k][L],Ajimlololo[k][L],
                          xirehihihi[L],xirelohihi[L],
                          xirehilohi[L],xirelolohi[L],
                          xirehihilo[L],xirelohilo[L],
                          xirehilolo[L],xirelololo[L],
                          &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                          &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  odf_dec(&wjimhihihi[k],&wjimlohihi[k],
                          &wjimhilohi[k],&wjimlolohi[k],
                          &wjimhihilo[k],&wjimlohilo[k],
                          &wjimhilolo[k],&wjimlololo[k],
                          acchihihi,acclohihi,acchilohi,acclolohi,
                          acchihilo,acclohilo,acchilolo,acclololo);
               }
         }
      }
      // compute sol[i] with back substitution
      double *xrehihihi = solrehihihi[i];
      double *xrelohihi = solrelohihi[i];
      double *xrehilohi = solrehilohi[i];
      double *xrelolohi = solrelolohi[i];
      double *xrehihilo = solrehihilo[i];
      double *xrelohilo = solrelohilo[i];
      double *xrehilolo = solrehilolo[i];
      double *xrelololo = solrelololo[i];
      double *ximhihihi = solimhihihi[i];
      double *ximlohihi = solimlohihi[i];
      double *ximhilohi = solimhilohi[i];
      double *ximlolohi = solimlolohi[i];
      double *ximhihilo = solimhihilo[i];
      double *ximlohilo = solimlohilo[i];
      double *ximhilolo = solimhilolo[i];
      double *ximlololo = solimlololo[i];
      double *brehihihi = rhsrehihihi[i];
      double *brelohihi = rhsrelohihi[i];
      double *brehilohi = rhsrehilohi[i];
      double *brelolohi = rhsrelolohi[i];
      double *brehihilo = rhsrehihilo[i];
      double *brelohilo = rhsrelohilo[i];
      double *brehilolo = rhsrehilolo[i];
      double *brelololo = rhsrelololo[i];
      double *bimhihihi = rhsimhihihi[i];
      double *bimlohihi = rhsimlohihi[i];
      double *bimhilohi = rhsimhilohi[i];
      double *bimlolohi = rhsimlolohi[i];
      double *bimhihilo = rhsimhihilo[i];
      double *bimlohilo = rhsimlohilo[i];
      double *bimhilolo = rhsimhilolo[i];
      double *bimlololo = rhsimlololo[i];

      CPU_cmplx_onenorm(dim,brehihihi,bimhihihi,&nrm);
      if(vrblvl > 0) cout << "1-norm of b : " << nrm << endl;

      if((nrm < 1.0e-120) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << i << "] ..." << endl;

         for(int j=0; j<dim; j++)
         {
            xrehihihi[j] = 0.0; xrelohihi[j] = 0.0;
            xrehilohi[j] = 0.0; xrelolohi[j] = 0.0;
            xrehihilo[j] = 0.0; xrelohilo[j] = 0.0;
            xrehilolo[j] = 0.0; xrelololo[j] = 0.0;
            ximhihihi[j] = 0.0; ximlohihi[j] = 0.0;
            ximhilohi[j] = 0.0; ximlolohi[j] = 0.0;
            ximhihilo[j] = 0.0; ximlohilo[j] = 0.0;
            ximhilolo[j] = 0.0; ximlololo[j] = 0.0;
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
         CPU_cmplx8_factors_qrbs
            (dim,dim,Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
                     Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
                     Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
                     Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
             Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
             Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
             Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
             Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
             brehihihi,brelohihi,brehilohi,brelolohi,
             brehihilo,brelohilo,brehilolo,brelololo,
             bimhihihi,bimlohihi,bimhilohi,bimlolohi,
             bimhihilo,bimlohilo,bimhilolo,bimlololo,
             xrehihihi,xrelohihi,xrehilohi,xrelolohi,
             xrehihilo,xrelohilo,xrehilolo,xrelololo,
             ximhihihi,ximlohihi,ximhilohi,ximlolohi,
             ximhihilo,ximlohilo,ximhilolo,ximlololo,
             wrkvecrehihihi,wrkvecrelohihi,wrkvecrehilohi,wrkvecrelolohi,
             wrkvecrehihilo,wrkvecrelohilo,wrkvecrehilolo,wrkvecrelololo,
             wrkvecimhihihi,wrkvecimlohihi,wrkvecimhilohi,wrkvecimlolohi,
             wrkvecimhihilo,wrkvecimlohilo,wrkvecimhilolo,wrkvecimlololo);

         if(vrblvl > 1)
         {
            for(int i=0; i<dim; i++)
               cout << "QHb[" << i << "] : "
                    << wrkvecrehihihi[i] << "  " << wrkvecrehihihi[i] << endl
                    << "  "
                    << wrkvecrehilohi[i] << "  " << wrkvecrehilohi[i] << endl
                    << "  "
                    << wrkvecrehihilo[i] << "  " << wrkvecrehihilo[i] << endl
                    << "  "
                    << wrkvecrehilolo[i] << "  " << wrkvecrehilolo[i] << endl
                    << "  "
                    << wrkvecimhihihi[i] << "  " << wrkvecimlohihi[i] << endl
                    << "  "
                    << wrkvecimhilohi[i] << "  " << wrkvecimlolohi[i] << endl
                    << "  "
                    << wrkvecimhihilo[i] << "  " << wrkvecimlohilo[i] << endl
                    << "  "
                    << wrkvecimhilolo[i] << "  " << wrkvecimlololo[i] << endl;

            cout << "the solution : " << endl;
            for(int j=0; j<dim; j++)
               cout << xrehihihi[j] << "  " << xrelohihi[j] << endl << "  "
                    << xrehilohi[j] << "  " << xrelolohi[j] << endl << "  "
                    << xrehihilo[j] << "  " << xrelohilo[j] << endl << "  "
                    << xrehilolo[j] << "  " << xrelololo[j] << endl << "  "
                    << ximhihihi[j] << "  " << ximlohihi[j] << endl << "  "
                    << ximhilohi[j] << "  " << ximlolohi[j] << endl << "  "
                    << ximhihilo[j] << "  " << ximlohilo[j] << endl << "  "
                    << ximhilolo[j] << "  " << ximlololo[j] << endl;
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

void CPU_dbl8_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_dbl8_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_dbl8_qrbs_head ..." << endl;

      CPU_dbl8_qrbs_head
         (dim,degp1,
          mathihihi,matlohihi,mathilohi,matlolohi,
          mathihilo,matlohilo,mathilolo,matlololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
          Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
          wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
          wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo,zeroQ,noqr,
          vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl8_qrbs_tail ..." << endl;

      CPU_dbl8_qrbs_tail
         (dim,degp1,tailidx,
          mathihihi,matlohihi,mathilohi,matlolohi,
          mathihilo,matlohilo,mathilolo,matlololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          Qhihihi,Qlohihi,Qhilohi,Qlolohi,Qhihilo,Qlohilo,Qhilolo,Qlololo,
          Rhihihi,Rlohihi,Rhilohi,Rlolohi,Rhihilo,Rlohilo,Rhilolo,Rlololo,
          wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
          wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo,
          upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_cmplx8_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi, 
   double ***matimhilohi, double ***matimlolohi, 
   double ***matimhihilo, double ***matimlohilo, 
   double ***matimhilolo, double ***matimlololo, 
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl )
{
   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping CPU_cmplx8_qrbs_head ..." << endl;
   }
   else
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx8_qrbs_head ..." << endl;

      CPU_cmplx8_qrbs_head
         (dim,degp1,matrehihihi,matrelohihi,matrehilohi,matrelolohi,
                    matrehihilo,matrelohilo,matrehilolo,matrelololo,
                    matimhihihi,matimlohihi,matimhilohi,matimlolohi,
                    matimhihilo,matimlohilo,matimhilolo,matimlololo,
          rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
          rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
          rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
          rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
          solrehihihi,solrelohihi,solrehilohi,solrelolohi,
          solrehihilo,solrelohilo,solrehilolo,solrelololo,
          solimhihihi,solimlohihi,solimhilohi,solimlolohi,
          solimhihilo,solimlohilo,solimhilolo,solimlololo,
          Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
          Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
          Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
          Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
          Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
          Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
          Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
          Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
          wrkvecrehihihi,wrkvecrelohihi,wrkvecrehilohi,wrkvecrelolohi,
          wrkvecrehihilo,wrkvecrelohilo,wrkvecrehilolo,wrkvecrelololo,
          wrkvecimhihihi,wrkvecimlohihi,wrkvecimhilohi,wrkvecimlolohi,
          wrkvecimhihilo,wrkvecimlohilo,wrkvecimhilolo,wrkvecimlololo,
          zeroQ,noqr,vrblvl);
   }
   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_cmplx8_qrbs_tail ..." << endl;

      CPU_cmplx8_qrbs_tail
         (dim,degp1,tailidx,
          matrehihihi,matrelohihi,matrehilohi,matrelolohi,
          matrehihilo,matrelohilo,matrehilolo,matrelololo,
          matimhihihi,matimlohihi,matimhilohi,matimlolohi,
          matimhihilo,matimlohilo,matimhilolo,matimlololo,
          rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
          rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
          rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
          rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
          solrehihihi,solrelohihi,solrehilohi,solrelolohi,
          solrehihilo,solrelohilo,solrehilolo,solrelololo,
          solimhihihi,solimlohihi,solimhilohi,solimlolohi,
          solimhihilo,solimlohilo,solimhilolo,solimlololo,
          Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
          Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
          Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
          Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
          Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
          Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
          Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
          Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
          wrkvecrehihihi,wrkvecrelohihi,wrkvecrehilohi,wrkvecrelolohi,
          wrkvecrehihilo,wrkvecrelohilo,wrkvecrehilolo,wrkvecrelololo,
          wrkvecimhihihi,wrkvecimlohihi,wrkvecimhilohi,wrkvecimlolohi,
          wrkvecimhihilo,wrkvecimlohilo,wrkvecimhilolo,wrkvecimlololo,
          upidx,bsidx,newtail,vrblvl);
   }
}

void CPU_dbl8_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl )
{
   *resmaxhihihi = 0.0;
   *resmaxlohihi = 0.0;
   *resmaxhilohi = 0.0;
   *resmaxlolohi = 0.0;
   *resmaxhihilo = 0.0;
   *resmaxlohilo = 0.0;
   *resmaxhilolo = 0.0;
   *resmaxlololo = 0.0;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rihihihi = resvechihihi[i];
      double *rilohihi = resveclohihi[i];
      double *rihilohi = resvechilohi[i];
      double *rilolohi = resveclolohi[i];
      double *rihihilo = resvechihilo[i];
      double *rilohilo = resveclohilo[i];
      double *rihilolo = resvechilolo[i];
      double *rilololo = resveclololo[i];

      for(int j=0; j<dim; j++)
      {
         rihihihi[j] = rhshihihi[i][j]; rilohihi[j] = rhslohihi[i][j];
         rihilohi[j] = rhshilohi[i][j]; rilolohi[j] = rhslolohi[i][j];
         rihihilo[j] = rhshihilo[i][j]; rilohilo[j] = rhslohilo[i][j];
         rihilolo[j] = rhshilolo[i][j]; rilololo[j] = rhslololo[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajhihihi = mathihihi[j];
         double **Ajlohihi = matlohihi[j];
         double **Ajhilohi = mathilohi[j];
         double **Ajlolohi = matlolohi[j];
         double **Ajhihilo = mathihilo[j];
         double **Ajlohilo = matlohilo[j];
         double **Ajhilolo = mathilolo[j];
         double **Ajlololo = matlololo[j];
         double *xhihihi = solhihihi[i-j];
         double *xlohihi = sollohihi[i-j];
         double *xhilohi = solhilohi[i-j];
         double *xlolohi = sollolohi[i-j];
         double *xhihilo = solhihilo[i-j];
         double *xlohilo = sollohilo[i-j];
         double *xhilolo = solhilolo[i-j];
         double *xlololo = sollololo[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               odf_mul(Ajhihihi[L][k],Ajlohihi[L][k],
                       Ajhilohi[L][k],Ajlolohi[L][k],
                       Ajhihilo[L][k],Ajlohilo[L][k],
                       Ajhilolo[L][k],Ajlololo[L][k],
                       xhihihi[k],xlohihi[k],xhilohi[k],xlolohi[k],
                       xhihilo[k],xlohilo[k],xhilolo[k],xlololo[k],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_dec(&rihihihi[L],&rilohihi[L],&rihilohi[L],&rilolohi[L],
                       &rihihilo[L],&rilohilo[L],&rihilolo[L],&rilololo[L],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << solhihihi[i][j] << "  " << sollohihi[i][j] << endl;
            cout << solhilohi[i][j] << "  " << sollolohi[i][j] << endl;
            cout << solhihilo[i][j] << "  " << sollohilo[i][j] << endl;
            cout << solhilolo[i][j] << "  " << sollololo[i][j] << endl;
         }
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << rihihihi[j] << "  " << rilohihi[j] << endl;
            cout << rihilohi[j] << "  " << rilolohi[j] << endl;
            cout << rihihilo[j] << "  " << rilohilo[j] << endl;
            cout << rihilolo[j] << "  " << rilololo[j] << endl;
         }
      }
      for(int j=0; j<dim; j++)
         if(abs(rihihihi[j]) > *resmaxhihihi)
         {
            *resmaxhihihi = abs(rihihihi[j]);
            *resmaxlohihi = abs(rilohihi[j]);
            *resmaxhilohi = abs(rihilohi[j]);
            *resmaxlolohi = abs(rilolohi[j]);
            *resmaxhihilo = abs(rihihilo[j]);
            *resmaxlohilo = abs(rilohilo[j]);
            *resmaxhilolo = abs(rihilolo[j]);
            *resmaxlololo = abs(rilololo[j]);
         }
   }
}

void CPU_cmplx8_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi, 
   double **rhsimhilohi, double **rhsimlolohi, 
   double **rhsimhihilo, double **rhsimlohilo, 
   double **rhsimhilolo, double **rhsimlololo, 
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl )
{
   *resmaxhihihi = 0.0;
   *resmaxlohihi = 0.0;
   *resmaxhilohi = 0.0;
   *resmaxlolohi = 0.0;
   *resmaxhihilo = 0.0;
   *resmaxlohilo = 0.0;
   *resmaxhilolo = 0.0;
   *resmaxlololo = 0.0;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=tailidx; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rirehihihi = resvecrehihihi[i];
      double *rirelohihi = resvecrelohihi[i];
      double *rirehilohi = resvecrehilohi[i];
      double *rirelolohi = resvecrelolohi[i];
      double *rirehihilo = resvecrehihilo[i];
      double *rirelohilo = resvecrelohilo[i];
      double *rirehilolo = resvecrehilolo[i];
      double *rirelololo = resvecrelololo[i];
      double *riimhihihi = resvecimhihihi[i];
      double *riimlohihi = resvecimlohihi[i];
      double *riimhilohi = resvecimhilohi[i];
      double *riimlolohi = resvecimlolohi[i];
      double *riimhihilo = resvecimhihilo[i];
      double *riimlohilo = resvecimlohilo[i];
      double *riimhilolo = resvecimhilolo[i];
      double *riimlololo = resvecimlololo[i];

      for(int j=0; j<dim; j++)
      {
         rirehihihi[j] = rhsrehihihi[i][j]; rirelohihi[j] = rhsrelohihi[i][j];
         rirehilohi[j] = rhsrehilohi[i][j]; rirelolohi[j] = rhsrelolohi[i][j];
         rirehihilo[j] = rhsrehihilo[i][j]; rirelohilo[j] = rhsrelohilo[i][j];
         rirehilolo[j] = rhsrehilolo[i][j]; rirelololo[j] = rhsrelololo[i][j];
         riimhihihi[j] = rhsimhihihi[i][j]; riimlohihi[j] = rhsimlohihi[i][j];
         riimhilohi[j] = rhsimhilohi[i][j]; riimlolohi[j] = rhsimlolohi[i][j];
         riimhihilo[j] = rhsimhihilo[i][j]; riimlohilo[j] = rhsimlohilo[i][j];
         riimhilolo[j] = rhsimhilolo[i][j]; riimlololo[j] = rhsimlololo[i][j];
      }
      for(int j=0; j<=(i-tailidx); j++)
      {
         double **Ajrehihihi = matrehihihi[j];
         double **Ajrelohihi = matrelohihi[j];
         double **Ajrehilohi = matrehilohi[j];
         double **Ajrelolohi = matrelolohi[j];
         double **Ajrehihilo = matrehihilo[j];
         double **Ajrelohilo = matrelohilo[j];
         double **Ajrehilolo = matrehilolo[j];
         double **Ajrelololo = matrelololo[j];
         double **Ajimhihihi = matimhihihi[j];
         double **Ajimlohihi = matimlohihi[j];
         double **Ajimhilohi = matimhilohi[j];
         double **Ajimlolohi = matimlolohi[j];
         double **Ajimhihilo = matimhihilo[j];
         double **Ajimlohilo = matimlohilo[j];
         double **Ajimhilolo = matimhilolo[j];
         double **Ajimlololo = matimlololo[j];
         double *xrehihihi = solrehihihi[i-j];
         double *xrelohihi = solrelohihi[i-j];
         double *xrehilohi = solrehilohi[i-j];
         double *xrelolohi = solrelolohi[i-j];
         double *xrehihilo = solrehihilo[i-j];
         double *xrelohilo = solrelohilo[i-j];
         double *xrehilolo = solrehilolo[i-j];
         double *xrelololo = solrelololo[i-j];
         double *ximhihihi = solimhihihi[i-j];
         double *ximlohihi = solimlohihi[i-j];
         double *ximhilohi = solimhilohi[i-j];
         double *ximlolohi = solimlolohi[i-j];
         double *ximhihilo = solimhihilo[i-j];
         double *ximlohilo = solimlohilo[i-j];
         double *ximhilolo = solimhilolo[i-j];
         double *ximlololo = solimlololo[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // ri[L] = ri[L] - Aj[L][k]*x[k];
            {
               // zre = Ajre[L][k]*xre[k] - Ajim[L][k]*xim[k];
               // rire[L] = rire[L] - zre;
               odf_mul(Ajrehihihi[L][k],Ajrelohihi[L][k],
                       Ajrehilohi[L][k],Ajrelolohi[L][k],
                       Ajrehihilo[L][k],Ajrelohilo[L][k],
                       Ajrehilolo[L][k],Ajrelololo[L][k],
                       xrehihihi[k],xrelohihi[k],xrehilohi[k],xrelolohi[k],
                       xrehihilo[k],xrelohilo[k],xrehilolo[k],xrelololo[k],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_dec(&rirehihihi[L],&rirelohihi[L],
                       &rirehilohi[L],&rirelolohi[L],
                       &rirehihilo[L],&rirelohilo[L],
                       &rirehilolo[L],&rirelololo[L],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
               odf_mul(Ajimhihihi[L][k],Ajimlohihi[L][k],
                       Ajimhilohi[L][k],Ajimlolohi[L][k],
                       Ajimhihilo[L][k],Ajimlohilo[L][k],
                       Ajimhilolo[L][k],Ajimlololo[L][k],
                       ximhihihi[k],ximlohihi[k],ximhilohi[k],ximlolohi[k],
                       ximhihilo[k],ximlohilo[k],ximhilolo[k],ximlololo[k],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_inc(&rirehihihi[L],&rirelohihi[L],
                       &rirehilohi[L],&rirelolohi[L],
                       &rirehihilo[L],&rirelohilo[L],
                       &rirehilolo[L],&rirelololo[L],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
               // zim = Ajre[L][k]*xim[k] + Ajim[L][k]*xre[k];
               // riim[L] = riim[L] - zim;
               odf_mul(Ajrehihihi[L][k],Ajrelohihi[L][k],
                       Ajrehilohi[L][k],Ajrelolohi[L][k],
                       Ajrehihilo[L][k],Ajrelohilo[L][k],
                       Ajrehilolo[L][k],Ajrelololo[L][k],
                       ximhihihi[k],ximlohihi[k],ximhilohi[k],ximlolohi[k],
                       ximhihilo[k],ximlohilo[k],ximhilolo[k],ximlololo[k],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_dec(&riimhihihi[L],&riimlohihi[L],
                       &riimhilohi[L],&riimlolohi[L],
                       &riimhihilo[L],&riimlohilo[L],
                       &riimhilolo[L],&riimlololo[L],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
               odf_mul(Ajimhihihi[L][k],Ajimlohihi[L][k],
                       Ajimhilohi[L][k],Ajimlolohi[L][k],
                       Ajimhihilo[L][k],Ajimlohilo[L][k],
                       Ajimhilolo[L][k],Ajimlololo[L][k],
                       xrehihihi[k],xrelohihi[k],xrehilohi[k],xrelolohi[k],
                       xrehihilo[k],xrelohilo[k],xrehilolo[k],xrelololo[k],
                       &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                       &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_dec(&riimhihihi[L],&riimlohihi[L],
                       &riimhilohi[L],&riimlolohi[L],
                       &riimhihilo[L],&riimlohilo[L],
                       &riimhilolo[L],&riimlololo[L],
                       acchihihi,acclohihi,acchilohi,acclolohi,
                       acchihilo,acclohilo,acchilolo,acclololo);
            }
      }
      if(vrblvl > 1)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solrehihihi[i][j] << "  " << solrelohihi[i][j] << endl
                 << "  "
                 << solrehilohi[i][j] << "  " << solrelolohi[i][j] << endl
                 << "  "
                 << solrehihilo[i][j] << "  " << solrelohilo[i][j] << endl
                 << "  "
                 << solrehilolo[i][j] << "  " << solrelololo[i][j] << endl
                 << "  "
                 << solimhihihi[i][j] << "  " << solimlohihi[i][j] << endl
                 << "  "
                 << solimhilohi[i][j] << "  " << solimlolohi[i][j] << endl
                 << "  "
                 << solimhihilo[i][j] << "  " << solimlohilo[i][j] << endl
                 << "  "
                 << solimhilolo[i][j] << "  " << solimlololo[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rirehihihi[j] << "  " << rirelohihi[j] << endl << "  "
                 << rirehilohi[j] << "  " << rirelolohi[j] << endl << "  "
                 << rirehihilo[j] << "  " << rirelohilo[j] << endl << "  "
                 << rirehilolo[j] << "  " << rirelololo[j] << endl << "  "
                 << riimhihihi[j] << "  " << riimlohihi[j] << endl << "  "
                 << riimhilohi[j] << "  " << riimlolohi[j] << endl << "  "
                 << riimhihilo[j] << "  " << riimlohilo[j] << endl << "  "
                 << riimhilolo[j] << "  " << riimlololo[j] << endl;
      }
      for(int j=0; j<dim; j++)
         if(abs(rirehihihi[j]) + abs(riimhihihi[j]) > *resmaxhihihi)
         {
            *resmaxhihihi = abs(rirehihihi[j]) + abs(riimhihihi[j]);
            *resmaxlohihi = abs(rirelohihi[j]) + abs(riimlohihi[j]);
            *resmaxhilohi = abs(rirehilohi[j]) + abs(riimhilohi[j]);
            *resmaxlolohi = abs(rirelolohi[j]) + abs(riimlolohi[j]);
            *resmaxhihilo = abs(rirehihilo[j]) + abs(riimhihilo[j]);
            *resmaxlohilo = abs(rirelohilo[j]) + abs(riimlohilo[j]);
            *resmaxhilolo = abs(rirehilolo[j]) + abs(riimhilolo[j]);
            *resmaxlololo = abs(rirelololo[j]) + abs(riimlololo[j]);
         }
   }
}
