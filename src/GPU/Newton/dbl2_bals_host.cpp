// The file dbl2_bals_host.cpp defines functions with prototypes in
// the file dbl2_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "dbl2_factorizations.h"
#include "dbl2_bals_host.h"

using namespace std;

void CPU_dbl2_lusb_head
 ( int dim, int degp1, double ***mathi, double ***matlo, 
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, double *wrkvechi, double *wrkveclo,
   int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         wrkmathi[i][j] = mathi[0][i][j];
         wrkmatlo[i][j] = matlo[0][i][j];
      }
   for(int i=0; i<dim; i++)
   {
      wrkvechi[i] = rhshi[0][i];
      wrkveclo[i] = rhslo[0][i];
   }
   if(verbose)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            cout << " " << wrkmathi[i][j] << " " << wrkmatlo[i][j];
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << wrkvechi[i] << "  " << wrkveclo[i] << endl;
   }
   if(verbose) cout << "calling CPU_dbl2_factors_lusolve ..." << endl;

   CPU_dbl2_factors_lusolve
      (dim,wrkmathi,wrkmatlo,pivots,wrkvechi,wrkveclo,solhi[0],sollo[0]);

   if(verbose)
   {
      double acchi,acclo;

      cout << "The pivots :";
      for(int i=0; i<dim; i++) cout << " " << pivots[i];
      cout << endl;
      cout << "The leading coefficients of the solution :" << endl;
      for(int i=0; i<dim; i++)
         cout << solhi[0][i] << "  " << sollo[0][i] << endl;

      for(int i=0; i<dim; i++)
      {
         wrkvechi[i] = rhshi[0][i];
         wrkveclo[i] = rhslo[0][i];

         for(int j=0; j<dim; j++)
            // wrkvec[i] = wrkvec[i] - mat[0][i][j]*sol[0][j];
            ddf_mul(mathi[0][i][j],matlo[0][i][j],
                    solhi[0][j],sollo[0][j],&acchi,&acclo);
            ddf_dec(&wrkvechi[i],&wrkveclo[i],acchi,acclo);
      }
      cout << "The residual vector :" << endl;
      for(int i=0; i<dim; i++)
         cout << wrkvechi[i] << "  " << wrkveclo[i] << endl;
   }
}

void CPU_dbl2_qrbs_head
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double *wrkvechi, double *wrkveclo,
   int vrblvl )
{
   bool verbose = (vrblvl > 0);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         wrkmathi[i][j] = mathi[0][i][j];
         wrkmatlo[i][j] = matlo[0][i][j];
      }

   if(verbose)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "  " << wrkmathi[i][j]
                 << "  " << wrkmatlo[i][j];
         }
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++)
         cout << rhshi[0][i] << "  " << rhslo[0][i] << endl;
   }
   if(verbose) cout << "calling CPU_dbl2_factors_houseqr ..." << endl;

   CPU_dbl2_factors_houseqr(dim,dim,wrkmathi,wrkmatlo,Qhi,Qlo,Rhi,Rlo);
   CPU_dbl2_factors_qrbs
      (dim,dim,Qhi,Qlo,Rhi,Rlo,rhshi[0],rhslo[0],solhi[0],sollo[0],
       wrkvechi,wrkveclo);

   if(verbose)
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

void CPU_dbl2_lusb_tail
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   double acchi,acclo;

   for(int i=1; i<degp1; i++)
   {
      if(verbose) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      for(int j=i; j<degp1; j++)
      {
         double **Ajhi = mathi[j-i+1]; // always start with A[1]
         double **Ajlo = matlo[j-i+1]; 
         double *xihi = solhi[i-1];    // solution to do the update with
         double *xilo = sollo[i-1];
         double *wjhi = rhshi[j];      // current right hand side vector
         double *wjlo = rhslo[j];

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
            {
               ddf_mul(Ajhi[k][L],Ajlo[k][L],xihi[L],xilo[L],&acchi,&acclo);
               ddf_dec(&wjhi[k],&wjlo[k],acchi,acclo);
            }
      }
      // compute sol[i] with back substitution
      double *xhi = solhi[i];
      double *xlo = sollo[i];
      double *bhi = rhshi[i];
      double *blo = rhslo[i];
      // the rhs[i] is used as work space
      for(int j=0; j<dim; j++)
      {
         xhi[j] = bhi[pivots[j]];
         xlo[j] = blo[pivots[j]];
      }
      CPU_dbl2_factors_forward(dim,wrkmathi,wrkmatlo,xhi,xlo,bhi,blo);
      CPU_dbl2_factors_backward(dim,wrkmathi,wrkmatlo,bhi,blo,xhi,xlo);
      if(verbose)
      {
         cout << "the solution : " << endl;
         for(int j=0; j<dim; j++)
            cout << xhi[j] << "  " << xlo[j] << endl;
      }
   }
}

void CPU_dbl2_qrbs_tail
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   double acchi,acclo;

   for(int i=1; i<degp1; i++)
   {
      if(verbose) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      for(int j=i; j<degp1; j++)
      {
         double **Ajhi = mathi[j-i+1]; // always start with A[1]
         double **Ajlo = matlo[j-i+1]; 
         double *xihi = solhi[i-1];    // solution to do the update with
         double *xilo = sollo[i-1];
         double *wjhi = rhshi[j];      // current right hand side vector
         double *wjlo = rhslo[j];

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) // wj[k] = wj[k] - Aj[k][L]*xi[L];
            {
               ddf_mul(Ajhi[k][L],Ajlo[k][L],xihi[L],xilo[L],&acchi,&acclo);
               ddf_dec(&wjhi[k],&wjlo[k],acchi,acclo);
            }
      }
      // compute sol[i] with back substitution
      double *xhi = solhi[i];
      double *xlo = sollo[i];
      double *bhi = rhshi[i];
      double *blo = rhslo[i];

      CPU_dbl2_factors_qrbs
         (dim,dim,Qhi,Qlo,Rhi,Rlo,bhi,blo,xhi,xlo,wrkvechi,wrkveclo);

      if(verbose)
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

void CPU_dbl2_lusb_solve
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, double *wrkvechi, double *wrkveclo,
   int *pivots, int vrblvl )
{
   if(vrblvl > 0) cout << "calling dbl2_linear_solve_head ..." << endl;

   CPU_dbl2_lusb_head
      (dim,degp1,mathi,matlo,rhshi,rhslo,solhi,sollo,
       wrkmathi,wrkmatlo,wrkvechi,wrkveclo,pivots,vrblvl);

   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling dbl2_linear_solve_tail ..." << endl;

      CPU_dbl2_lusb_tail
         (dim,degp1,mathi,matlo,rhshi,rhslo,solhi,sollo,
          wrkmathi,wrkmatlo,pivots,vrblvl);
   }
}

void CPU_dbl2_qrbs_solve
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo, int vrblvl )
{
   if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_head ..." << endl;

   CPU_dbl2_qrbs_head
      (dim,degp1,mathi,matlo,rhshi,rhslo,solhi,sollo,wrkmathi,wrkmatlo,
       Qhi,Qlo,Rhi,Rlo,wrkvechi,wrkveclo,vrblvl);

   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_tail ..." << endl;

      CPU_dbl2_qrbs_tail
         (dim,degp1,mathi,matlo,rhshi,rhslo,solhi,sollo,Qhi,Qlo,Rhi,Rlo,
          wrkvechi,wrkveclo,vrblvl);
   }
}

void CPU_dbl2_linear_residue
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   int vrblvl )
{
   *resmaxhi = 0.0;
   *resmaxlo = 0.0;
   double acchi,acclo;

   for(int i=0; i<degp1; i++)  // compute the i-th residual vector
   {
      double *rihi = resvechi[i];
      double *rilo = resveclo[i];
      for(int j=0; j<dim; j++)
      {
         rihi[j] = rhshi[i][j];
         rilo[j] = rhslo[i][j];
      }
      for(int j=0; j<=i; j++)
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
      if(vrblvl > 0)
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
