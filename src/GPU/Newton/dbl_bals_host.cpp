// The file dbl_bals_host.cpp defines functions with prototypes in
// the file dbl_bals_host.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "dbl_factorizations.h"
#include "dbl_bals_host.h"

using namespace std;

void CPU_dbl_lusb_head
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **wrkmat, double *wrkvec, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) wrkmat[i][j] = mat[0][i][j];

   for(int i=0; i<dim; i++) wrkvec[i] = rhs[0][i];

   if(verbose)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++) cout << " " << wrkmat[i][j];
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++) cout << wrkvec[i] << endl;
   }
   if(verbose) cout << "calling CPU_dbl_factors_lusolve ..." << endl;

   CPU_dbl_factors_lusolve(dim,wrkmat,pivots,wrkvec,sol[0]);

   if(verbose)
   {
      cout << "The pivots :";
      for(int i=0; i<dim; i++) cout << " " << pivots[i];
      cout << endl;
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

void CPU_dbl_qrbs_head
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **wrkmat, double **Q, double **R, double *wrkvec, int vrblvl )
{
   bool verbose = (vrblvl > 0);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) wrkmat[i][j] = mat[0][i][j];

   if(verbose)
   {
      cout << "The matrix : " << endl;
      cout << setprecision(2);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++) cout << " " << wrkmat[i][j];
         cout << endl;
      }
      cout << setprecision(16);
      cout << "The right hand side vector : " << endl;
      for(int i=0; i<dim; i++) cout << rhs[0][i] << endl;
   }
   if(verbose) cout << "calling CPU_dbl_factors_houseqr ..." << endl;

   CPU_dbl_factors_houseqr(dim,dim,wrkmat,Q,R);
   CPU_dbl_factors_qrbs(dim,dim,Q,R,rhs[0],sol[0],wrkvec);

   if(verbose)
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

void CPU_dbl_lusb_tail
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **wrkmat, int *pivots, int vrblvl )
{
   bool verbose = (vrblvl > 0);

   for(int i=1; i<degp1; i++)
   {
      if(verbose) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      for(int j=i; j<degp1; j++)
      {
         double **Aj = mat[j-i+1]; // always start with A[1]
         double *xi = sol[i-1];    // solution to do the update with
         double *wj = rhs[j];      // current right hand side vector

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) wj[k] = wj[k] - Aj[k][L]*xi[L];
      }
      // compute sol[i] with back substitution
      double *x = sol[i];
      double *b = rhs[i];
      // the rhs[i] is used as work space
      for(int j=0; j<dim; j++) x[j] = b[pivots[j]];
      CPU_dbl_factors_forward(dim,wrkmat,x,b);
      CPU_dbl_factors_backward(dim,wrkmat,b,x);
      if(verbose)
      {
         cout << "the solution : " << endl;
         for(int j=0; j<dim; j++) cout << x[j] << endl;
      }
   }
}

void CPU_dbl_qrbs_tail
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec, int vrblvl )
{
   bool verbose = (vrblvl > 0);

   for(int i=1; i<degp1; i++)
   {
      if(verbose) cout << "stage " << i << " in solve tail ..." << endl;
      // use sol[i-1] to update rhs[j] for j in i to degp1
      for(int j=i; j<degp1; j++)
      {
         double **Aj = mat[j-i+1]; // always start with A[1]
         double *xi = sol[i-1];    // solution to do the update with
         double *wj = rhs[j];      // current right hand side vector

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) wj[k] = wj[k] - Aj[k][L]*xi[L];
      }
      // compute sol[i] with back substitution
      double *x = sol[i];
      double *b = rhs[i];

      CPU_dbl_factors_qrbs(dim,dim,Q,R,b,x,wrkvec);

      if(verbose)
      {
         for(int i=0; i<dim; i++)
            cout << "Qtb[" << i << "] : " << b[i] << endl;

         cout << "the solution : " << endl;
         for(int j=0; j<dim; j++) cout << x[j] << endl;
      }
   }
}

void CPU_dbl_lusb_solve
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **wrkmat, double *wrkvec, int *pivots, int vrblvl )
{
   if(vrblvl > 0) cout << "calling CPU_dbl_lusb_head ..." << endl;

   CPU_dbl_lusb_head
      (dim,degp1,mat,rhs,sol,wrkmat,wrkvec,pivots,vrblvl);

   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl_lusb_tail ..." << endl;

      CPU_dbl_lusb_tail(dim,degp1,mat,rhs,sol,wrkmat,pivots,vrblvl);
   }
}

void CPU_dbl_qrbs_solve
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **wrkmat, double **Q, double **R, double *wrkvec, int vrblvl )
{
   if(vrblvl > 0) cout << "calling CPU_dbl_qrbs_head ..." << endl;

   CPU_dbl_qrbs_head
      (dim,degp1,mat,rhs,sol,wrkmat,Q,R,wrkvec,vrblvl);

   if(degp1 > 1)
   {
      if(vrblvl > 0) cout << "calling CPU_dbl_qrbs_tail ..." << endl;

      CPU_dbl_qrbs_tail(dim,degp1,mat,rhs,sol,Q,R,wrkvec,vrblvl);
   }
}

void CPU_dbl_linear_residue
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **resvec, double *resmax, int vrblvl )
{
   *resmax = 0.0;

   for(int i=0; i<degp1; i++)  // compute the i-th residual vector
   {
      double *ri = resvec[i];
      for(int j=0; j<dim; j++) ri[j] = rhs[i][j];
      for(int j=0; j<=i; j++)
      {
         double **Aj = mat[j];
         double *x = sol[i-j];

         // if(vrblvl > 0)
         //    cout << "A[" << j << "] and x[" << i-j << "] ..." << endl;

         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) ri[L] = ri[L] - Aj[L][k]*x[k];
      }
      if(vrblvl > 0)
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
