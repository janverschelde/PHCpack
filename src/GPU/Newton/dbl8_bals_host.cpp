// The file dbl8_bals_host.cpp defines functions with prototypes in
// the file dbl8_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"
#include "dbl8_bals_host.h"

using namespace std;

void CPU_dbl8_linear_head
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         wrkmathihihi[i][j] = mathihihi[0][i][j];
         wrkmatlohihi[i][j] = matlohihi[0][i][j];
         wrkmathilohi[i][j] = mathilohi[0][i][j];
         wrkmatlolohi[i][j] = matlolohi[0][i][j];
         wrkmathihilo[i][j] = mathihilo[0][i][j];
         wrkmatlohilo[i][j] = matlohilo[0][i][j];
         wrkmathilolo[i][j] = mathilolo[0][i][j];
         wrkmatlololo[i][j] = matlololo[0][i][j];
      }
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
   }
   if(verbose)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << " " << wrkmathihihi[i][j] << " " << wrkmatlohihi[i][j];
            cout << " " << wrkmathilohi[i][j] << " " << wrkmatlolohi[i][j];
            cout << " " << wrkmathihilo[i][j] << " " << wrkmatlohilo[i][j];
            cout << " " << wrkmathilolo[i][j] << " " << wrkmatlololo[i][j];
         }
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << wrkvechihihi[i] << "  " << wrkveclohihi[i] << endl;
         cout << wrkvechilohi[i] << "  " << wrkveclolohi[i] << endl;
         cout << wrkvechihilo[i] << "  " << wrkveclohilo[i] << endl;
         cout << wrkvechilolo[i] << "  " << wrkveclololo[i] << endl;
      }
   }
   if(verbose) cout << "calling CPU_dbl4_factors_lusolve ..." << endl;

   CPU_dbl8_factors_lusolve
      (dim,wrkmathihihi,wrkmatlohihi,wrkmathilohi,wrkmatlolohi,
           wrkmathihilo,wrkmatlohilo,wrkmathilolo,wrkmatlololo,pivots,
       wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
       wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo,
       solhihihi[0],sollohihi[0],solhilohi[0],sollolohi[0],
       solhihilo[0],sollohilo[0],solhilolo[0],sollololo[0]);

   if(verbose)
   {
      double acchihihi,acclohihi,acchilohi,acclolohi;
      double acchihilo,acclohilo,acchilolo,acclololo;

      cout << "The pivots :";
      for(int i=0; i<dim; i++) cout << " " << pivots[i];
      cout << endl;
      cout << "The leading coefficients of the solution :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << solhihihi[0][i] << "  " << sollohihi[0][i] << endl;
         cout << solhilohi[0][i] << "  " << sollolohi[0][i] << endl;
         cout << solhihilo[0][i] << "  " << sollohilo[0][i] << endl;
         cout << solhilolo[0][i] << "  " << sollololo[0][i] << endl;
      }
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
      cout << "The residual vector :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << wrkvechihihi[i] << "  " << wrkveclohihi[i] << endl;
         cout << wrkvechilohi[i] << "  " << wrkveclolohi[i] << endl;
         cout << wrkvechihilo[i] << "  " << wrkveclohilo[i] << endl;
         cout << wrkvechilolo[i] << "  " << wrkveclololo[i] << endl;
      }
   }
}

void CPU_dbl8_linear_tail
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=1; i<degp1; i++)
   {
      if(verbose) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
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
         double *xihihihi = solhihihi[i-1]; // solution to do the update with
         double *xilohihi = sollohihi[i-1];
         double *xihilohi = solhilohi[i-1];
         double *xilolohi = sollolohi[i-1];
         double *xihihilo = solhihilo[i-1];
         double *xilohilo = sollohilo[i-1];
         double *xihilolo = solhilolo[i-1];
         double *xilololo = sollololo[i-1];
         double *wjhihihi = rhshihihi[j];  // current right hand side vector
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
      // the rhs[i] is used as work space
      for(int j=0; j<dim; j++)
      {
         xhihihi[j] = bhihihi[pivots[j]];
         xlohihi[j] = blohihi[pivots[j]];
         xhilohi[j] = bhilohi[pivots[j]];
         xlolohi[j] = blolohi[pivots[j]];
         xhihilo[j] = bhihilo[pivots[j]];
         xlohilo[j] = blohilo[pivots[j]];
         xhilolo[j] = bhilolo[pivots[j]];
         xlololo[j] = blololo[pivots[j]];
      }
      CPU_dbl8_factors_forward
         (dim,wrkmathihihi,wrkmatlohihi,wrkmathilohi,wrkmatlolohi,
              wrkmathihilo,wrkmatlohilo,wrkmathilolo,wrkmatlololo,
          xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
          bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo);
      CPU_dbl8_factors_backward
         (dim,wrkmathihihi,wrkmatlohihi,wrkmathilohi,wrkmatlolohi,
              wrkmathihilo,wrkmatlohilo,wrkmathilolo,wrkmatlololo,
          bhihihi,blohihi,bhilohi,blolohi,bhihilo,blohilo,bhilolo,blololo,
          xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo);

      if(verbose)
      {
         cout << "the solution : " << endl;
         for(int j=0; j<dim; j++)
         {
            cout << xhihihi[j] << "  " << xlohihi[j] << endl;
            cout << xhilohi[j] << "  " << xlolohi[j] << endl;
            cout << xhihilo[j] << "  " << xlohilo[j] << endl;
            cout << xhilolo[j] << "  " << xlololo[j] << endl;
         }
      }
   }
}

void CPU_dbl8_linear_solve
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo, int *pivots, int vrblvl )
{
   if(vrblvl > 0) cout << "calling CPU_dbl8_linear_head ..." << endl;

   CPU_dbl8_linear_head
      (dim,degp1,mathihihi,matlohihi,mathilohi,matlolohi,
                 mathihilo,matlohilo,mathilolo,matlololo,
       rhshihihi,rhslohihi,rhshilohi,rhslolohi,
       rhshihilo,rhslohilo,rhshilolo,rhslololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,
       wrkmathihihi,wrkmatlohihi,wrkmathilohi,wrkmatlolohi,
       wrkmathihilo,wrkmatlohilo,wrkmathilolo,wrkmatlololo,
       wrkvechihihi,wrkveclohihi,wrkvechilohi,wrkveclolohi,
       wrkvechihilo,wrkveclohilo,wrkvechilolo,wrkveclololo,pivots,vrblvl);

   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl8_linear_tail ..." << endl;

      CPU_dbl8_linear_tail
         (dim,degp1,mathihihi,matlohihi,mathilohi,matlolohi,
                    mathihilo,matlohilo,mathilolo,matlololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          wrkmathihihi,wrkmatlohihi,wrkmathilohi,wrkmatlolohi,
          wrkmathihilo,wrkmatlohilo,wrkmathilolo,wrkmatlololo,pivots,vrblvl);
   }
}

void CPU_dbl8_linear_residue
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

   for(int i=0; i<degp1; i++)  // compute the i-th residual vector
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
      for(int j=0; j<=i; j++)
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
      if(vrblvl > 0)
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
