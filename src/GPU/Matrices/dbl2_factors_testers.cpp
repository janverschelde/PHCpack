/* The file dbl2_factors_testers.cpp define the functions specified in
   the file dbl2_factors_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "random2_matrices.h"
#include "dbl2_factorizations.h"

using namespace std;

void test_factors_real2_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Ahi = new double*[dim];
   double **Alo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }
   random_dbl2_matrix(dim,dim,Ahi,Alo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   double *solhi = new double[dim];
   double *sollo = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solhi[i] = 1.0;
      sollo[i] = 0.0;
   }
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];
   double acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshi[i] << "  " << rhslo[i] << endl;
   }
   double *xhi = new double[dim];
   double *xlo = new double[dim];
   int *pivots = new int[dim];

   CPU_dbl2_factors_lusolve(dim,Ahi,Alo,pivots,rhshi,rhslo,xhi,xlo);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhi[i] << "  " << xlo[i] << endl;
   }
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xhi[i] - 1.0) + abs(xlo[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;
}

void test_factors_cmplx2_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Arehi = new double*[dim];
   double **Arelo = new double*[dim];
   double **Aimhi = new double*[dim];
   double **Aimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Arehi[i] = new double[dim];
      Arelo[i] = new double[dim];
      Aimhi[i] = new double[dim];
      Aimlo[i] = new double[dim];
   }
   random_cmplx2_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi[i][j] << "  " << Arelo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
         }
   }
   double *solrehi = new double[dim];
   double *solrelo = new double[dim];
   double *solimhi = new double[dim];
   double *solimlo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solrehi[i] = 1.0; solrelo[i] = 0.0;
      solimhi[i] = 0.0; solimlo[i] = 0.0;
   }
   double *rhsrehi = new double[dim];
   double *rhsrelo = new double[dim];
   double *rhsimhi = new double[dim];
   double *rhsimlo = new double[dim];
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;

   for(int i=0; i<dim; i++)
   {
      rhsrehi[i] = 0.0; rhsrelo[i] = 0.0;
      rhsimhi[i] = 0.0; rhsimlo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Arehi[i][j],Arelo[i][j],solrehi[j],solrelo[j],
                 &acc1hi,&acc1lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solimhi[j],solimlo[j],
                 &acc2hi,&acc2lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solrehi[j],solrelo[j],
                 &acc3hi,&acc3lo);
         ddf_mul(Arehi[i][j],Arelo[i][j],solimhi[j],solimlo[j],
                 &acc4hi,&acc4lo);
         ddf_inc(&rhsrehi[i],&rhsrelo[i],acc1hi,acc1lo);
         ddf_dec(&rhsrehi[i],&rhsrelo[i],acc2hi,acc2lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc3hi,acc3lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc4hi,acc4lo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehi[i] << "  " << rhsrelo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhi[i] << "  " << rhsimlo[i] << endl;
      }
   }
   double *xrehi = new double[dim];
   double *xrelo = new double[dim];
   double *ximhi = new double[dim];
   double *ximlo = new double[dim];
   int *pivots = new int[dim];

   CPU_cmplx2_factors_lusolve
      (dim,Arehi,Arelo,Aimhi,Aimlo,pivots,
       rhsrehi,rhsrelo,rhsimhi,rhsimlo,xrehi,xrelo,ximhi,ximlo);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehi[i] << "  " << xrelo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhi[i] << "  " << ximlo[i] << endl;
      }
   }
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xrehi[i] - 1.0) + abs(xrelo[i])
                    + abs(ximhi[i]) + abs(ximlo[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;
}

int test_real2_qr_factors_probe
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumhi,Qsumlo,Rsumhi,Rsumlo,acchi,acclo;
   double errorQ = 0.0;
   double errorR = 0.0;

   for(int p=0; p<nbprobes; p++)
   {
      rowidx = rand() % nrows;
      colidx = rand() % ncols;

      if(verbose > 0)
      {
         cout << "Probing row index : " << rowidx
              << ", column index : " << colidx << "." << endl;
      }
      Qsumhi = 0.0; Qsumlo = 0.0;
      Rsumhi = 0.0; Rsumlo = 0.0;

      for(int i=0; i<nrows; i++)
      {
         ddf_mul(Qhi[i][rowidx],Qlo[i][rowidx],
                 Qhi[i][colidx],Qlo[i][colidx],&acchi,&acclo);
         ddf_inc(&Qsumhi,&Qsumlo,acchi,acclo);
         ddf_mul(Qhi[i][rowidx],Qlo[i][rowidx],
                 Ahi[i][colidx],Alo[i][colidx],&acchi,&acclo);
         ddf_inc(&Rsumhi,&Rsumlo,acchi,acclo);
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "] : "
              << Qsumhi << "  " << Qsumlo << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "] : "
              << Rsumhi << "  " << Rsumlo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "] : "
              << Rhi[rowidx][colidx] << "  "
              << Rlo[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
         errorQ = errorQ + fabs(Qsumhi - 1.0) + fabs(Qsumlo);
      else
         errorQ = errorQ + fabs(Qsumhi) + fabs(Qsumlo);
      errorR = errorR + fabs(Rsumhi - Rhi[rowidx][colidx])
                      + fabs(Rsumlo - Rlo[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

int test_real2_qr_factors
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double tol, int verbose )
{
   double **QThi = new double*[nrows];
   double **QTlo = new double*[nrows];
   double **QTQhi = new double*[nrows];
   double **QTQlo = new double*[nrows];
   double **QTAhi = new double*[nrows];
   double **QTAlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QThi[i] = new double[nrows];
      QTlo[i] = new double[nrows];
      QTQhi[i] = new double[nrows];
      QTQlo[i] = new double[nrows];
      QTAhi[i] = new double[ncols];
      QTAlo[i] = new double[ncols];
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix Q :" << endl;
   }
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhi[i][j] << "  " << Qlo[i][j] << endl;
         QThi[j][i] = Qhi[i][j];
         QTlo[j][i] = Qlo[i][j];
      }
   CPU_dbl2_factors_matmatmul
      (nrows,nrows,nrows,QThi,QTlo,Qhi,Qlo,QTQhi,QTQlo);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q'*Q[" << i << "][" << j << "] : "
                 << QTQhi[i][j] << "  " << QTQlo[i][j] << endl;
         if(i == j)
            errorQ = errorQ + fabs(QTQhi[i][j] - 1.0) + fabs(QTQlo[i][j]);
         else
            errorQ = errorQ + fabs(QTQhi[i][j]) + fabs(QTQlo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhi[i][j] << "  " << Rlo[i][j] << endl;
   }
   CPU_dbl2_factors_matmatmul
      (nrows,nrows,ncols,QThi,QTlo,Ahi,Alo,QTAhi,QTAlo);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q'*A[" << i << "][" << j << "] : "
                 << QTAhi[i][j] << "  " << QTAlo[i][j] << endl;
         errorR = errorR + fabs(Rhi[i][j] - QTAhi[i][j])
                         + fabs(Rlo[i][j] - QTAlo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QThi[i]); free(QTQhi[i]); free(QTAhi[i]);
      free(QTlo[i]); free(QTQlo[i]); free(QTAlo[i]);
   }
   free(QThi); free(QTQhi); free(QTAhi);
   free(QTlo); free(QTQlo); free(QTAlo);

   return int(errorQ + errorR > tol);
}

int test_cmplx2_qr_factors_probe
 ( int nrows, int ncols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumrehi,Qsumrelo,Rsumrehi,Rsumrelo;
   double Qsumimhi,Qsumimlo,Rsumimhi,Rsumimlo;
   double acchi,acclo;
   double errorQ = 0.0;
   double errorR = 0.0;

   for(int p=0; p<nbprobes; p++)
   {
      rowidx = rand() % nrows;
      colidx = rand() % ncols;

      if(verbose > 0)
      {
         cout << "Probing row index : " << rowidx
              << ", column index : " << colidx << "." << endl;
      }
      Qsumrehi = 0.0; Qsumrelo = 0.0;
      Qsumimhi = 0.0; Qsumimlo = 0.0;
      Rsumrehi = 0.0; Rsumrelo = 0.0;
      Rsumimhi = 0.0; Rsumimlo = 0.0;

      for(int i=0; i<nrows; i++)
      {
         // multiply Q^H with Q
         ddf_mul(Qrehi[i][rowidx],Qrelo[i][rowidx],
                 Qrehi[i][colidx],Qrelo[i][colidx],&acchi,&acclo);
         ddf_inc(&Qsumrehi,&Qsumrelo,acchi,acclo);
         ddf_mul(Qimhi[i][rowidx],Qimlo[i][rowidx],
                 Qimhi[i][colidx],Qimlo[i][colidx],&acchi,&acclo);
         ddf_inc(&Qsumrehi,&Qsumrelo,acchi,acclo);
         ddf_mul(Qimhi[i][rowidx],Qimlo[i][rowidx],
                 Qrehi[i][colidx],Qrelo[i][colidx],&acchi,&acclo);
         ddf_dec(&Qsumimhi,&Qsumimlo,acchi,acclo);
         ddf_mul(Qrehi[i][rowidx],Qrelo[i][rowidx],
                 Qimhi[i][colidx],Qimlo[i][colidx],&acchi,&acclo);
         ddf_inc(&Qsumimhi,&Qsumimlo,acchi,acclo);
         // multiply Q^H with A
         ddf_mul(Qrehi[i][rowidx],Qrelo[i][rowidx],
                 Arehi[i][colidx],Arelo[i][colidx],&acchi,&acclo);
         ddf_inc(&Rsumrehi,&Rsumrelo,acchi,acclo);
         ddf_mul(Qimhi[i][rowidx],Qimlo[i][rowidx],
                 Aimhi[i][colidx],Aimlo[i][colidx],&acchi,&acclo);
         ddf_inc(&Rsumrehi,&Rsumrelo,acchi,acclo);
         ddf_mul(Qimhi[i][rowidx],Qimlo[i][rowidx],
                 Arehi[i][colidx],Arelo[i][colidx],&acchi,&acclo);
         ddf_dec(&Rsumimhi,&Rsumimlo,acchi,acclo);
         ddf_mul(Qrehi[i][rowidx],Qrelo[i][rowidx],
                 Aimhi[i][colidx],Aimlo[i][colidx],&acchi,&acclo);
         ddf_inc(&Rsumimhi,&Rsumimlo,acchi,acclo);
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "]re : "
              << Qsumrehi << "  " << Qsumrelo << endl;
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "]im : "
              << Qsumimhi << "  " << Qsumimlo << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "]re : "
              << Rsumrehi << "  " << Rsumrelo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "]re : "
              << Rrehi[rowidx][colidx] << "  "
              << Rrelo[rowidx][colidx] << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "]im : "
              << Rsumimhi << "  " << Rsumimlo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "]im : "
              << Rimhi[rowidx][colidx] << "  "
              << Rimlo[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
      {
         errorQ = errorQ + fabs(Qsumrehi - 1.0) + fabs(Qsumrelo)
                         + fabs(Qsumimhi)       + fabs(Qsumimlo);
      }
      else
      {
         errorQ = errorQ + fabs(Qsumrehi) + fabs(Qsumrelo)
                         + fabs(Qsumimhi) + fabs(Qsumimlo);
      }
      errorR = errorR + fabs(Rsumrehi - Rrehi[rowidx][colidx])
                      + fabs(Rsumrelo - Rrelo[rowidx][colidx])
                      + fabs(Rsumimhi - Rimhi[rowidx][colidx])
                      + fabs(Rsumimlo - Rimlo[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

int test_cmplx2_qr_factors
 ( int nrows, int ncols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double tol, int verbose )
{
   double **QHrehi = new double*[nrows];
   double **QHrelo = new double*[nrows];
   double **QHimhi = new double*[nrows];
   double **QHimlo = new double*[nrows];
   double **QHQrehi = new double*[nrows];
   double **QHQrelo = new double*[nrows];
   double **QHQimhi = new double*[nrows];
   double **QHQimlo = new double*[nrows];
   double **QHArehi = new double*[nrows];
   double **QHArelo = new double*[nrows];
   double **QHAimhi = new double*[nrows];
   double **QHAimlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QHrehi[i] = new double[nrows];
      QHrelo[i] = new double[nrows];
      QHimhi[i] = new double[nrows];
      QHimlo[i] = new double[nrows];
      QHQrehi[i] = new double[nrows];
      QHQrelo[i] = new double[nrows];
      QHQimhi[i] = new double[nrows];
      QHQimlo[i] = new double[nrows];
      QHArehi[i] = new double[ncols];
      QHArelo[i] = new double[ncols];
      QHAimhi[i] = new double[ncols];
      QHAimlo[i] = new double[ncols];
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix Q :" << endl;
   }
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehi[i][j] << "  " << Qrelo[i][j] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhi[i][j] << "  " << Qimlo[i][j] << endl;
         }
         QHrehi[j][i] = Qrehi[i][j]; QHimhi[j][i] = Qimhi[i][j];
         QHrelo[j][i] = Qrelo[i][j]; QHimlo[j][i] = Qimlo[i][j];
         ddf_minus(&QHimhi[j][i],&QHimlo[j][i]); // Hermitian transpose
      }

   CPU_cmplx2_factors_matmatmul
      (nrows,nrows,nrows,QHrehi, QHrelo, QHimhi, QHimlo,
                          Qrehi,  Qrelo,  Qimhi,  Qimlo,
                        QHQrehi,QHQrelo,QHQimhi,QHQimlo);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
         {
            cout << "Q^H*Q[" << i << "][" << j << "]re : "
                 << QHQrehi[i][j] << "  " << QHQrelo[i][j] << endl;
            cout << "Q^H*Q[" << i << "][" << j << "]im : "
                 << QHQimhi[i][j] << "  " << QHQimlo[i][j] << endl;
         }
         if(i == j)
            errorQ = errorQ + abs(QHQrehi[i][j] - 1.0) + abs(QHQrelo[i][j])
                            + abs(QHQimhi[i][j]) + abs(QHQimlo[i][j]);
         else
            errorQ = errorQ + abs(QHQrehi[i][j]) + abs(QHQrelo[i][j])
                            + abs(QHQimhi[i][j]) + abs(QHQimlo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*Q - I| : " << errorQ << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "R[" << i << "][" << j << "]re : "
                 << Rrehi[i][j] << "  " << Rrelo[i][j] << endl;
            cout << "R[" << i << "][" << j << "]im : "
                 << Rimhi[i][j] << "  " << Rimlo[i][j] << endl;
         }
   }
   CPU_cmplx2_factors_matmatmul
      (nrows,nrows,ncols,QHrehi, QHrelo, QHimhi, QHimlo,
                          Arehi,  Arelo,  Aimhi,  Aimlo,
                        QHArehi,QHArelo,QHAimhi,QHAimlo);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
         {
            cout << "Q^H*A[" << i << "][" << j << "]re : "
                 << QHArehi[i][j] << "  " << QHArelo[i][j] << endl;
            cout << "Q^H*A[" << i << "][" << j << "]im : "
                 << QHAimhi[i][j] << "  " << QHAimlo[i][j] << endl;
         }
         errorR = errorR + abs(Rrehi[i][j] - QHArehi[i][j])
                         + abs(Rrelo[i][j] - QHArelo[i][j])
                         + abs(Rimhi[i][j] - QHAimhi[i][j])
                         + abs(Rimlo[i][j] - QHAimlo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QHrehi[i]); free(QHimhi[i]);
      free(QHrelo[i]); free(QHimlo[i]);
      free(QHQrehi[i]); free(QHQimhi[i]);
      free(QHQrelo[i]); free(QHQimlo[i]);
      free(QHArehi[i]); free(QHAimhi[i]);
      free(QHArelo[i]); free(QHAimlo[i]);
   }
   free(QHrehi); free(QHQrehi); free(QHArehi);
   free(QHrelo); free(QHQrelo); free(QHArelo);
   free(QHimhi); free(QHQimhi); free(QHAimhi);
   free(QHimlo); free(QHQimlo); free(QHAimlo);

   return int(errorQ + errorR > tol);
}

void test_factors_real2_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Ahi = new double*[nrows];
   double **Alo = new double*[nrows];
   double **Qhi = new double*[nrows];
   double **Qlo = new double*[nrows];
   double **Rhi = new double*[nrows];
   double **Rlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahi[i] = new double[ncols];
      Alo[i] = new double[ncols];
      Qhi[i] = new double[nrows];
      Qlo[i] = new double[nrows];
      Rhi[i] = new double[ncols];
      Rlo[i] = new double[ncols];
   }
   random_dbl2_matrix(nrows,ncols,Ahi,Alo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   CPU_dbl2_factors_houseqr(nrows,ncols,Ahi,Alo,Qhi,Qlo,Rhi,Rlo);

   const double tol = 1.0E-26;
   const int fail = test_real2_qr_factors
      (nrows,ncols,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Ahi[i]); free(Qhi[i]); free(Rhi[i]);
      free(Alo[i]); free(Qlo[i]); free(Rlo[i]);
   }
   free(Ahi); free(Qhi); free(Rhi);
   free(Alo); free(Qlo); free(Rlo);
}

void test_factors_cmplx2_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Arehi = new double*[nrows];
   double **Arelo = new double*[nrows];
   double **Aimhi = new double*[nrows];
   double **Aimlo = new double*[nrows];
   double **Qrehi = new double*[nrows];
   double **Qrelo = new double*[nrows];
   double **Qimhi = new double*[nrows];
   double **Qimlo = new double*[nrows];
   double **Rrehi = new double*[nrows];
   double **Rrelo = new double*[nrows];
   double **Rimhi = new double*[nrows];
   double **Rimlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehi[i] = new double[ncols];
      Arelo[i] = new double[ncols];
      Aimhi[i] = new double[ncols];
      Aimlo[i] = new double[ncols];
      Qrehi[i] = new double[nrows];
      Qrelo[i] = new double[nrows];
      Qimhi[i] = new double[nrows];
      Qimlo[i] = new double[nrows];
      Rrehi[i] = new double[ncols];
      Rrelo[i] = new double[ncols];
      Rimhi[i] = new double[ncols];
      Rimlo[i] = new double[ncols];
   }
   random_cmplx2_matrix(nrows,ncols,Arehi,Arelo,Aimhi,Aimlo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi[i][j] << "  " << Arelo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
         }
   }
   CPU_cmplx2_factors_houseqr
      (nrows,ncols,Arehi,Arelo,Aimhi,Aimlo,
                   Qrehi,Qrelo,Qimhi,Qimlo,
                   Rrehi,Rrelo,Rimhi,Rimlo);

   const double tol = 1.0e-12;
   const int fail = test_cmplx2_qr_factors
      (nrows,ncols,Arehi,Arelo,Aimhi,Aimlo,
                   Qrehi,Qrelo,Qimhi,Qimlo,
                   Rrehi,Rrelo,Rimhi,Rimlo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Arehi[i]); free(Qrehi[i]); free(Rrehi[i]);
      free(Arelo[i]); free(Qrelo[i]); free(Rrelo[i]);
      free(Aimhi[i]); free(Qimhi[i]); free(Rimhi[i]);
      free(Aimlo[i]); free(Qimlo[i]); free(Rimlo[i]);
   }
   free(Arehi); free(Qrehi); free(Rrehi);
   free(Arelo); free(Qrelo); free(Rrelo);
   free(Aimhi); free(Qimhi); free(Rimhi);
   free(Aimlo); free(Qimlo); free(Rimlo);
}
