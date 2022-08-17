// The file dbl4_bals_host.cpp defines functions with prototypes in
// the file dbl4_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"
#include "dbl4_bals_host.h"

using namespace std;

void CPU_dbl4_linear_head
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **wrkmathihi, double **wrkmatlohi,
   double **wrkmathilo, double **wrkmatlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         wrkmathihi[i][j] = mathihi[0][i][j];
         wrkmatlohi[i][j] = matlohi[0][i][j];
         wrkmathilo[i][j] = mathilo[0][i][j];
         wrkmatlolo[i][j] = matlolo[0][i][j];
      }
   for(int i=0; i<dim; i++)
   {
      wrkvechihi[i] = rhshihi[0][i];
      wrkveclohi[i] = rhslohi[0][i];
      wrkvechilo[i] = rhshilo[0][i];
      wrkveclolo[i] = rhslolo[0][i];
   }
   if(verbose)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << " " << wrkmathihi[i][j] << " " << wrkmatlohi[i][j];
            cout << " " << wrkmathilo[i][j] << " " << wrkmatlolo[i][j];
         }
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << wrkvechihi[i] << "  " << wrkveclohi[i] << endl;
         cout << wrkvechilo[i] << "  " << wrkveclolo[i] << endl;
      }
   }
   if(verbose) cout << "calling CPU_dbl4_factors_lusolve ..." << endl;

   CPU_dbl4_factors_lusolve
      (dim,wrkmathihi,wrkmatlohi,wrkmathilo,wrkmatlolo,pivots,
       wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo,
       solhihi[0],sollohi[0],solhilo[0],sollolo[0]);

   if(verbose)
   {
      double acchihi,acclohi,acchilo,acclolo;

      cout << "The pivots :";
      for(int i=0; i<dim; i++) cout << " " << pivots[i];
      cout << endl;
      cout << "The leading coefficients of the solution :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << solhihi[0][i] << "  " << sollohi[0][i] << endl;
         cout << solhilo[0][i] << "  " << sollolo[0][i] << endl;
      }
      for(int i=0; i<dim; i++)
      {
         wrkvechihi[i] = rhshihi[0][i];
         wrkveclohi[i] = rhslohi[0][i];
         wrkvechilo[i] = rhshilo[0][i];
         wrkveclolo[i] = rhslolo[0][i];

         for(int j=0; j<dim; j++)
            // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
            qdf_mul(mathihi[0][i][j],matlohi[0][i][j],
                    mathilo[0][i][j],matlolo[0][i][j],
                    solhihi[0][j],sollohi[0][j],solhilo[0][j],sollolo[0][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_dec(&wrkvechihi[i],&wrkveclohi[i],
                    &wrkvechilo[i],&wrkveclolo[i],
                    acchihi,acclohi,acchilo,acclolo);
      }
      cout << "The residual vector :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << wrkvechihi[i] << "  " << wrkveclohi[i] << endl;
         cout << wrkvechilo[i] << "  " << wrkveclolo[i] << endl;
      }
   }
}

void CPU_dbl4_linear_tail
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **wrkmathihi, double **wrkmatlohi,
   double **wrkmathilo, double **wrkmatlolo, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=1; i<degp1; i++)
   {
      if(verbose) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      for(int j=i; j<degp1; j++)
      {
         double **Ajhihi = mathihi[j-i+1]; // always start with A[1]
         double **Ajlohi = matlohi[j-i+1]; 
         double **Ajhilo = mathilo[j-i+1];
         double **Ajlolo = matlolo[j-i+1]; 
         double *xihihi = solhihi[i-1];    // solution to do the update with
         double *xilohi = sollohi[i-1];
         double *xihilo = solhilo[i-1];
         double *xilolo = sollolo[i-1];
         double *wjhihi = rhshihi[j];      // current right hand side vector
         double *wjlohi = rhslohi[j];
         double *wjhilo = rhshilo[j];
         double *wjlolo = rhslolo[j];

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
            {
               qdf_mul(Ajhihi[k][L],Ajlohi[k][L],
                       Ajhilo[k][L],Ajlolo[k][L],
                       xihihi[L],xilohi[L],xihilo[L],xilolo[L],
                       &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_dec(&wjhihi[k],&wjlohi[k],&wjhilo[k],&wjlolo[k],
                       acchihi,acclohi,acchilo,acclolo);
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
      // the rhs[i] is used as work space
      for(int j=0; j<dim; j++)
      {
         xhihi[j] = bhihi[pivots[j]];
         xlohi[j] = blohi[pivots[j]];
         xhilo[j] = bhilo[pivots[j]];
         xlolo[j] = blolo[pivots[j]];
      }
      CPU_dbl4_factors_forward
         (dim,wrkmathihi,wrkmatlohi,wrkmathilo,wrkmatlolo,
          xhihi,xlohi,xhilo,xlolo,bhihi,blohi,bhilo,blolo);
      CPU_dbl4_factors_backward
         (dim,wrkmathihi,wrkmatlohi,wrkmathilo,wrkmatlolo,
          bhihi,blohi,bhilo,blolo,xhihi,xlohi,xhilo,xlolo);
      if(verbose)
      {
         cout << "the solution : " << endl;
         for(int j=0; j<dim; j++)
         {
            cout << xhihi[j] << "  " << xlohi[j] << endl;
            cout << xhilo[j] << "  " << xlolo[j] << endl;
         }
      }
   }
}

void CPU_dbl4_linear_solve
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **wrkmathihi, double **wrkmatlohi,
   double **wrkmathilo, double **wrkmatlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo, int *pivots, int vrblvl )
{
   if(vrblvl > 0) cout << "calling CPU_dbl4_linear_head ..." << endl;

   CPU_dbl4_linear_head
      (dim,degp1,mathihi,matlohi,mathilo,matlolo,
       rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
       wrkmathihi,wrkmatlohi,wrkmathilo,wrkmatlolo,
       wrkvechihi,wrkveclohi,wrkvechilo,wrkveclolo,pivots,vrblvl);

   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl4_linear_tail ..." << endl;

      CPU_dbl4_linear_tail
         (dim,degp1,mathihi,matlohi,mathilo,matlolo,
          rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
          wrkmathihi,wrkmatlohi,wrkmathilo,wrkmatlolo,pivots,vrblvl);
   }
}

void CPU_dbl4_linear_residue
 ( int dim, int degp1,
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

   for(int i=0; i<degp1; i++)  // compute the i-th residual vector
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
      for(int j=0; j<=i; j++)
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
      if(vrblvl > 0)
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
