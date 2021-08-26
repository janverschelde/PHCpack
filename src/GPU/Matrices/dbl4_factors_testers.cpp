/* The file dbl4_factors_testers.cpp define the functions specified in
   the file dbl4_factors_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "quad_double_functions.h"
#include "random4_matrices.h"
#include "dbl4_factorizations.h"

using namespace std;

void test_factors_real4_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Ahihi = new double*[dim];
   double **Alohi = new double*[dim];
   double **Ahilo = new double*[dim];
   double **Alolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Ahihi[i] = new double[dim];
      Alohi[i] = new double[dim];
      Ahilo[i] = new double[dim];
      Alolo[i] = new double[dim];
   }
   random_dbl4_matrix(dim,dim,Ahihi,Alohi,Ahilo,Alolo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi[i][j] << "  " << Alohi[i][j] << endl
                 << "          "
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }
   double *solhihi = new double[dim];
   double *sollohi = new double[dim];
   double *solhilo = new double[dim];
   double *sollolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solhihi[i] = 1.0; sollohi[i] = 0.0;
      solhilo[i] = 0.0; sollolo[i] = 0.0;
   }
   double *rhshihi = new double[dim];
   double *rhslohi = new double[dim];
   double *rhshilo = new double[dim];
   double *rhslolo = new double[dim];
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<dim; i++)
   {
      rhshihi[i] = 0.0; rhslohi[i] = 0.0;
      rhshilo[i] = 0.0; rhslolo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Ahihi[i][j],Alohi[i][j],Ahilo[i][j],Alolo[i][j],
               solhihi[j], sollohi[j], solhilo[j], sollolo[j],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_inc(&rhshihi[i],&rhslohi[i],&rhshilo[i],&rhslolo[i],
                  acchihi,    acclohi,    acchilo,    acclolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshihi[i] << "  " << rhslohi[i] << endl
              << "       "
              << rhshilo[i] << "  " << rhslolo[i] << endl;
   }
   double *xhihi = new double[dim];
   double *xlohi = new double[dim];
   double *xhilo = new double[dim];
   double *xlolo = new double[dim];
   int *pivots = new int[dim];

   CPU_dbl4_factors_lusolve
      (dim,Ahihi,  Alohi,  Ahilo,  Alolo,pivots,
         rhshihi,rhslohi,rhshilo,rhslolo,xhihi,xlohi,xhilo,xlolo);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihi[i] << "  " << xlohi[i] << endl
              << "       "
              << xhilo[i] << "  " << xlolo[i] << endl;
   }
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xhihi[i] - 1.0) + abs(xlohi[i])
                    + abs(xhilo[i]) + abs(xlolo[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;

   for(int i=0; i<dim; i++)
   {
      free(Ahihi[i]); free(Alohi[i]);
      free(Ahilo[i]); free(Alolo[i]);
   }
   free(Ahihi); free(Alohi); free(Ahilo); free(Alolo);
   free(solhihi); free(sollohi); free(solhilo); free(sollolo);
   free(rhshihi); free(rhslohi); free(rhshilo); free(rhslolo);
   free(xhihi); free(xlohi); free(xhilo); free(xlolo);
}

void test_factors_cmplx4_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Arehihi = new double*[dim];
   double **Arelohi = new double*[dim];
   double **Arehilo = new double*[dim];
   double **Arelolo = new double*[dim];
   double **Aimhihi = new double*[dim];
   double **Aimlohi = new double*[dim];
   double **Aimhilo = new double*[dim];
   double **Aimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Arehihi[i] = new double[dim];
      Arelohi[i] = new double[dim];
      Arehilo[i] = new double[dim];
      Arelolo[i] = new double[dim];
      Aimhihi[i] = new double[dim];
      Aimlohi[i] = new double[dim];
      Aimhilo[i] = new double[dim];
      Aimlolo[i] = new double[dim];
   }
   random_cmplx4_matrix
      (dim,dim,Arehihi,Arelohi,Arehilo,Arelolo,
               Aimhihi,Aimlohi,Aimhilo,Aimlolo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << "            "
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << "            "
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   double *solrehihi = new double[dim];
   double *solrehilo = new double[dim];
   double *solrelohi = new double[dim];
   double *solrelolo = new double[dim];
   double *solimhihi = new double[dim];
   double *solimhilo = new double[dim];
   double *solimlohi = new double[dim];
   double *solimlolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solrehihi[i] = 1.0; solrelohi[i] = 0.0;
      solrehilo[i] = 0.0; solrelolo[i] = 0.0;
      solimhihi[i] = 0.0; solimlohi[i] = 0.0;
      solimhilo[i] = 0.0; solimlolo[i] = 0.0;
   }
   double *rhsrehihi = new double[dim];
   double *rhsrehilo = new double[dim];
   double *rhsrelohi = new double[dim];
   double *rhsrelolo = new double[dim];
   double *rhsimhihi = new double[dim];
   double *rhsimhilo = new double[dim];
   double *rhsimlohi = new double[dim];
   double *rhsimlolo = new double[dim];
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;

   for(int i=0; i<dim; i++)
   {
      rhsrehihi[i] = 0.0; rhsrelohi[i] = 0.0;
      rhsrehilo[i] = 0.0; rhsrelolo[i] = 0.0;
      rhsimhihi[i] = 0.0; rhsimlohi[i] = 0.0;
      rhsimhilo[i] = 0.0; rhsimlolo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc2hihi,    &acc2lohi,    &acc2hilo,    &acc2lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc3hihi,    &acc3lohi,    &acc3hilo,    &acc3lolo);
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc4hihi,    &acc4lohi,    &acc4hilo,    &acc4lolo);

         qdf_inc(&rhsrehihi[i],&rhsrelohi[i],&rhsrehilo[i],&rhsrelolo[i],
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
         qdf_dec(&rhsrehihi[i],&rhsrelohi[i],&rhsrehilo[i],&rhsrelolo[i],
                   acc2hihi,     acc2lohi,     acc2hilo,     acc2lolo);
         qdf_inc(&rhsimhihi[i],&rhsimlohi[i],&rhsimhilo[i],&rhsimlolo[i],
                   acc3hihi,     acc3lohi,     acc3hilo,     acc3lolo);
         qdf_inc(&rhsimhihi[i],&rhsimlohi[i],&rhsimhilo[i],&rhsimlolo[i],
                   acc4hihi,     acc4lohi,     acc4hilo,     acc4lolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihi[i] << "  " << rhsrelohi[i] << endl
              << "         "
              << rhsrehilo[i] << "  " << rhsrelolo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihi[i] << "  " << rhsimlohi[i] << endl
              << "         "
              << rhsimhilo[i] << "  " << rhsimlolo[i] << endl;
      }
   }
   double *xrehihi = new double[dim];
   double *xrelohi = new double[dim];
   double *xrehilo = new double[dim];
   double *xrelolo = new double[dim];
   double *ximhihi = new double[dim];
   double *ximlohi = new double[dim];
   double *ximhilo = new double[dim];
   double *ximlolo = new double[dim];
   int *pivots = new int[dim];

   CPU_cmplx4_factors_lusolve
      (dim,Arehihi,  Arelohi,  Arehilo,  Arelolo,
           Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,pivots,
         rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
         rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
           xrehihi,  xrelohi,  xrehilo,  xrelolo,
           ximhihi,  ximlohi,  ximhilo,  ximlolo);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihi[i] << "  " << xrelohi[i] << endl
              << "         "
              << xrehilo[i] << "  " << xrelolo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihi[i] << "  " << ximlohi[i] << endl
              << "         "
              << ximhilo[i] << "  " << ximlolo[i] << endl;
      }
   }
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xrehihi[i] - 1.0) + abs(xrelohi[i])
                    + abs(xrehilo[i]) + abs(xrelolo[i])
                    + abs(ximhihi[i]) + abs(ximlohi[i])
                    + abs(ximhilo[i]) + abs(ximlolo[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;

   for(int i=0; i<dim; i++)
   {
      free(Arehihi[i]); free(Arelohi[i]);
      free(Arehilo[i]); free(Arelolo[i]);
      free(Aimhihi[i]); free(Aimlohi[i]);
      free(Aimhilo[i]); free(Aimlolo[i]);
   }
   free(Arehihi); free(Arelohi); free(Arehilo); free(Arelolo);
   free(Aimhihi); free(Aimlohi); free(Aimhilo); free(Aimlolo);
   free(solrehihi); free(solrelohi); free(solrehilo); free(solrelolo);
   free(solimhihi); free(solimlohi); free(solimhilo); free(solimlolo);
   free(rhsrehihi); free(rhsrelohi); free(rhsrehilo); free(rhsrelolo);
   free(rhsimhihi); free(rhsimlohi); free(rhsimhilo); free(rhsimlolo);
   free(xrehihi); free(xrelohi); free(xrehilo); free(xrelolo);
   free(ximhihi); free(ximlohi); free(ximhilo); free(ximlolo);
}

int test_real4_qr_factors_probe
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumhihi,Qsumlohi,Qsumhilo,Qsumlolo;
   double Rsumhihi,Rsumlohi,Rsumhilo,Rsumlolo;
   double acchihi,acclohi,acchilo,acclolo;
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
      Qsumhihi = 0.0; Qsumlohi = 0.0;
      Qsumhilo = 0.0; Qsumlolo = 0.0;
      Rsumhihi = 0.0; Rsumlohi = 0.0;
      Rsumhilo = 0.0; Rsumlolo = 0.0;

      for(int i=0; i<nrows; i++)
      {
         qdf_mul(Qhihi[i][rowidx],Qlohi[i][rowidx],
                 Qhilo[i][rowidx],Qlolo[i][rowidx],
                 Qhihi[i][colidx],Qlohi[i][colidx],
                 Qhilo[i][colidx],Qlolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Qsumhihi,&Qsumlohi,&Qsumhilo,&Qsumlolo,
                   acchihi,  acclohi,  acchilo,  acclolo);
         qdf_mul(Qhihi[i][rowidx],Qlohi[i][rowidx],
                 Qhilo[i][rowidx],Qlolo[i][rowidx],
                 Ahihi[i][colidx],Alohi[i][colidx],
                 Ahilo[i][colidx],Alolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Rsumhihi,&Rsumlohi,&Rsumhilo,&Rsumlolo,
                   acchihi,  acclohi,  acchilo,  acclolo);
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "] : "
              << Qsumhihi << "  " << Qsumlohi << endl
              << "              "
              << Qsumhilo << "  " << Qsumlolo << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "] : "
              << Rsumhihi << "  " << Rsumlohi << endl
              << "              "
              << Rsumhilo << "  " << Rsumlolo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "] : "
              << Rhihi[rowidx][colidx] << "  "
              << Rlohi[rowidx][colidx] << endl
              << "              "
              << Rhilo[rowidx][colidx] << "  "
              << Rlolo[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
         errorQ = errorQ + fabs(Qsumhihi - 1.0) + fabs(Qsumlohi)
                         + fabs(Qsumhilo) + fabs(Qsumlolo);
      else
         errorQ = errorQ + fabs(Qsumhihi) + fabs(Qsumlohi)
                         + fabs(Qsumhilo) + fabs(Qsumlolo);

      errorR = errorR + fabs(Rsumhihi - Rhihi[rowidx][colidx])
                      + fabs(Rsumlohi - Rlohi[rowidx][colidx])
                      + fabs(Rsumhilo - Rhilo[rowidx][colidx])
                      + fabs(Rsumlolo - Rlolo[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

int test_real4_qr_factors
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double tol, int verbose )
{
   double **QThihi = new double*[nrows];
   double **QTlohi = new double*[nrows];
   double **QThilo = new double*[nrows];
   double **QTlolo = new double*[nrows];
   double **QTQhihi = new double*[nrows];
   double **QTQlohi = new double*[nrows];
   double **QTQhilo = new double*[nrows];
   double **QTQlolo = new double*[nrows];
   double **QTAhihi = new double*[nrows];
   double **QTAlohi = new double*[nrows];
   double **QTAhilo = new double*[nrows];
   double **QTAlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QThihi[i] = new double[nrows];
      QTlohi[i] = new double[nrows];
      QThilo[i] = new double[nrows];
      QTlolo[i] = new double[nrows];
      QTQhihi[i] = new double[nrows];
      QTQlohi[i] = new double[nrows];
      QTQhilo[i] = new double[nrows];
      QTQlolo[i] = new double[nrows];
      QTAhihi[i] = new double[ncols];
      QTAlohi[i] = new double[ncols];
      QTAhilo[i] = new double[ncols];
      QTAlolo[i] = new double[ncols];
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
                 << Qhihi[i][j] << "  " << Qlohi[i][j] << endl
                 << "          "
                 << Qhilo[i][j] << "  " << Qlolo[i][j] << endl;
         QThihi[j][i] = Qhihi[i][j];
         QTlohi[j][i] = Qlohi[i][j];
         QThilo[j][i] = Qhilo[i][j];
         QTlolo[j][i] = Qlolo[i][j];
      }
   CPU_dbl4_factors_matmatmul
      (nrows,nrows,nrows,QThihi, QTlohi, QThilo, QTlolo,
                          Qhihi,  Qlohi,  Qhilo,  Qlolo,
                        QTQhihi,QTQlohi,QTQhilo,QTQlolo);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q^T*Q[" << i << "][" << j << "] : "
                 << QTQhihi[i][j] << "  " << QTQlohi[i][j] << endl
                 << "              "
                 << QTQhilo[i][j] << "  " << QTQlolo[i][j] << endl;

         if(i == j)
            errorQ = errorQ + fabs(QTQhihi[i][j] - 1.0) + fabs(QTQlohi[i][j])
                            + fabs(QTQhilo[i][j]) + fabs(QTQlolo[i][j]);
         else
            errorQ = errorQ + fabs(QTQhihi[i][j]) + fabs(QTQlohi[i][j])
                            + fabs(QTQhilo[i][j]) + fabs(QTQlolo[i][j]);
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
                 << Rhihi[i][j] << "  " << Rlohi[i][j] << endl
                 << "          "
                 << Rhilo[i][j] << "  " << Rlolo[i][j] << endl;
   }
   CPU_dbl4_factors_matmatmul
      (nrows,nrows,ncols,QThihi, QTlohi, QThilo, QTlolo,
                          Ahihi,  Alohi,  Ahilo,  Alolo,
                        QTAhihi,QTAlohi,QTAhilo,QTAlolo);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q^T*A[" << i << "][" << j << "] : "
                 << QTAhihi[i][j] << "  " << QTAlohi[i][j] << endl
                 << "              "
                 << QTAhilo[i][j] << "  " << QTAlolo[i][j] << endl;

         errorR = errorR + fabs(Rhihi[i][j] - QTAhihi[i][j])
                         + fabs(Rhilo[i][j] - QTAhilo[i][j])
                         + fabs(Rlohi[i][j] - QTAlohi[i][j])
                         + fabs(Rlolo[i][j] - QTAlolo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QThihi[i]); free(QTQhihi[i]); free(QTAhihi[i]);
      free(QTlohi[i]); free(QTQlohi[i]); free(QTAlohi[i]);
      free(QThilo[i]); free(QTQhilo[i]); free(QTAhilo[i]);
      free(QTlolo[i]); free(QTQlolo[i]); free(QTAlolo[i]);
   }
   free(QThihi); free(QTQhihi); free(QTAhihi);
   free(QTlohi); free(QTQlohi); free(QTAlohi);
   free(QThilo); free(QTQhilo); free(QTAhilo);
   free(QTlolo); free(QTQlolo); free(QTAlolo);

   return int(errorQ + errorR > tol);
}

int test_cmplx4_qr_factors_probe
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumrehihi,Qsumrelohi,Qsumrehilo,Qsumrelolo;
   double Rsumrehihi,Rsumrelohi,Rsumrehilo,Rsumrelolo;
   double Qsumimhihi,Qsumimlohi,Qsumimhilo,Qsumimlolo;
   double Rsumimhihi,Rsumimlohi,Rsumimhilo,Rsumimlolo;
   double acchihi,acclohi,acchilo,acclolo;
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
      Qsumrehihi = 0.0; Qsumrelohi = 0.0;
      Qsumrehilo = 0.0; Qsumrelolo = 0.0;
      Qsumimhihi = 0.0; Qsumimlohi = 0.0;
      Qsumimhilo = 0.0; Qsumimlolo = 0.0;
      Rsumrehihi = 0.0; Rsumrelohi = 0.0;
      Rsumrehilo = 0.0; Rsumrelolo = 0.0;
      Rsumimhihi = 0.0; Rsumimlohi = 0.0;
      Rsumimhilo = 0.0; Rsumimlolo = 0.0;

      for(int i=0; i<nrows; i++)
      {
         // multiply Q^H with Q
         qdf_mul(Qrehihi[i][rowidx],Qrelohi[i][rowidx],
                 Qrehilo[i][rowidx],Qrelolo[i][rowidx],
                 Qrehihi[i][colidx],Qrelohi[i][colidx],
                 Qrehilo[i][colidx],Qrelolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Qsumrehihi,&Qsumrelohi,&Qsumrehilo,&Qsumrelolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         qdf_mul(Qimhihi[i][rowidx],Qimlohi[i][rowidx],
                 Qimhilo[i][rowidx],Qimlolo[i][rowidx],
                 Qimhihi[i][colidx],Qimlohi[i][colidx],
                 Qimhilo[i][colidx],Qimlolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Qsumrehihi,&Qsumrelohi,&Qsumrehilo,&Qsumrelolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         qdf_mul(Qimhihi[i][rowidx],Qimlohi[i][rowidx],
                 Qimhilo[i][rowidx],Qimlolo[i][rowidx],
                 Qrehihi[i][colidx],Qrelohi[i][colidx],
                 Qrehilo[i][colidx],Qrelolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_dec(&Qsumimhihi,&Qsumimlohi,&Qsumimhilo,&Qsumimlolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         qdf_mul(Qrehihi[i][rowidx],Qrelohi[i][rowidx],
                 Qrehilo[i][rowidx],Qrelolo[i][rowidx],
                 Qimhihi[i][colidx],Qimlohi[i][colidx],
                 Qimhilo[i][colidx],Qimlolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Qsumimhihi,&Qsumimlohi,&Qsumimhilo,&Qsumimlolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         // multiply Q^H with A
         qdf_mul(Qrehihi[i][rowidx],Qrelohi[i][rowidx],
                 Qrehilo[i][rowidx],Qrelolo[i][rowidx],
                 Arehihi[i][colidx],Arelohi[i][colidx],
                 Arehilo[i][colidx],Arelolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Rsumrehihi,&Rsumrelohi,&Rsumrehilo,&Rsumrelolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         qdf_mul(Qimhihi[i][rowidx],Qimlohi[i][rowidx],
                 Qimhilo[i][rowidx],Qimlolo[i][rowidx],
                 Aimhihi[i][colidx],Aimlohi[i][colidx],
                 Aimhilo[i][colidx],Aimlolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Rsumrehihi,&Rsumrelohi,&Rsumrehilo,&Rsumrelolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         qdf_mul(Qimhihi[i][rowidx],Qimlohi[i][rowidx],
                 Qimhilo[i][rowidx],Qimlolo[i][rowidx],
                 Arehihi[i][colidx],Arelohi[i][colidx],
                 Arehilo[i][colidx],Arelolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_dec(&Rsumimhihi,&Rsumimlohi,&Rsumimhilo,&Rsumimlolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
         qdf_mul(Qrehihi[i][rowidx],Qrelohi[i][rowidx],
                 Qrehilo[i][rowidx],Qrelolo[i][rowidx],
                 Aimhihi[i][colidx],Aimlohi[i][colidx],
                 Aimhilo[i][colidx],Aimlolo[i][colidx],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_inc(&Rsumimhihi,&Rsumimlohi,&Rsumimhilo,&Rsumimlolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "]re : "
              << Qsumrehihi << "  " << Qsumrelohi << endl
              << "                "
              << Qsumrehilo << "  " << Qsumrelolo << endl;
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "]im : "
              << Qsumimhihi << "  " << Qsumimlohi << endl
              << "                "
              << Qsumimhilo << "  " << Qsumimlolo << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "]re : "
              << Rsumrehihi << "  " << Rsumrelohi << endl
              << "                "
              << Rsumrehilo << "  " << Rsumrelolo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "]re : "
              << Rrehihi[rowidx][colidx] << "  "
              << Rrelohi[rowidx][colidx] << endl
              << "                "
              << Rrehilo[rowidx][colidx] << "  "
              << Rrelolo[rowidx][colidx] << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "]im : "
              << Rsumimhihi << "  " << Rsumimlohi << endl
              << "                "
              << Rsumimhilo << "  " << Rsumimlolo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "]im : "
              << Rimhihi[rowidx][colidx] << "  "
              << Rimlohi[rowidx][colidx] << endl
              << "                "
              << Rimhilo[rowidx][colidx] << "  "
              << Rimlolo[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
      {
         errorQ = errorQ + fabs(Qsumrehihi - 1.0) + fabs(Qsumrelohi)
                         + fabs(Qsumrehilo)       + fabs(Qsumrelolo)
                         + fabs(Qsumimhihi)       + fabs(Qsumimlohi)
                         + fabs(Qsumimhilo)       + fabs(Qsumimlolo);
      }
      else
      {
         errorQ = errorQ + fabs(Qsumrehihi) + fabs(Qsumrelohi)
                         + fabs(Qsumrehilo) + fabs(Qsumrelolo)
                         + fabs(Qsumimhihi) + fabs(Qsumimlohi)
                         + fabs(Qsumimhilo) + fabs(Qsumimlolo);
      }
      errorR = errorR + fabs(Rsumrehihi - Rrehihi[rowidx][colidx])
                      + fabs(Rsumrehilo - Rrehilo[rowidx][colidx])
                      + fabs(Rsumrelohi - Rrelohi[rowidx][colidx])
                      + fabs(Rsumrelolo - Rrelolo[rowidx][colidx])
                      + fabs(Rsumimhihi - Rimhihi[rowidx][colidx])
                      + fabs(Rsumimhilo - Rimhilo[rowidx][colidx])
                      + fabs(Rsumimlohi - Rimlohi[rowidx][colidx])
                      + fabs(Rsumimlolo - Rimlolo[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

int test_cmplx4_qr_factors
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double tol, int verbose )
{
   double **QHrehihi = new double*[nrows];
   double **QHrelohi = new double*[nrows];
   double **QHrehilo = new double*[nrows];
   double **QHrelolo = new double*[nrows];
   double **QHimhihi = new double*[nrows];
   double **QHimlohi = new double*[nrows];
   double **QHimhilo = new double*[nrows];
   double **QHimlolo = new double*[nrows];
   double **QHQrehihi = new double*[nrows];
   double **QHQrelohi = new double*[nrows];
   double **QHQrehilo = new double*[nrows];
   double **QHQrelolo = new double*[nrows];
   double **QHQimhihi = new double*[nrows];
   double **QHQimlohi = new double*[nrows];
   double **QHQimhilo = new double*[nrows];
   double **QHQimlolo = new double*[nrows];
   double **QHArehihi = new double*[nrows];
   double **QHArelohi = new double*[nrows];
   double **QHArehilo = new double*[nrows];
   double **QHArelolo = new double*[nrows];
   double **QHAimhihi = new double*[nrows];
   double **QHAimlohi = new double*[nrows];
   double **QHAimhilo = new double*[nrows];
   double **QHAimlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QHrehihi[i] = new double[nrows];
      QHrelohi[i] = new double[nrows];
      QHrehilo[i] = new double[nrows];
      QHrelolo[i] = new double[nrows];
      QHimhihi[i] = new double[nrows];
      QHimlohi[i] = new double[nrows];
      QHimhilo[i] = new double[nrows];
      QHimlolo[i] = new double[nrows];
      QHQrehihi[i] = new double[nrows];
      QHQrelohi[i] = new double[nrows];
      QHQrehilo[i] = new double[nrows];
      QHQrelolo[i] = new double[nrows];
      QHQimhihi[i] = new double[nrows];
      QHQimlohi[i] = new double[nrows];
      QHQimhilo[i] = new double[nrows];
      QHQimlolo[i] = new double[nrows];
      QHArehihi[i] = new double[ncols];
      QHArelohi[i] = new double[ncols];
      QHArehilo[i] = new double[ncols];
      QHArelolo[i] = new double[ncols];
      QHAimhihi[i] = new double[ncols];
      QHAimlohi[i] = new double[ncols];
      QHAimhilo[i] = new double[ncols];
      QHAimlolo[i] = new double[ncols];
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
                 << Qrehihi[i][j] << "  " << Qrelohi[i][j] << endl
                 << "            "
                 << Qrehilo[i][j] << "  " << Qrelolo[i][j] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhihi[i][j] << "  " << Qimlohi[i][j] << endl
                 << "            "
                 << Qimhilo[i][j] << "  " << Qimlolo[i][j] << endl;
         }
         QHrehihi[j][i] = Qrehihi[i][j]; QHrelohi[j][i] = Qrelohi[i][j];
         QHrehilo[j][i] = Qrehilo[i][j]; QHrelolo[j][i] = Qrelolo[i][j];
         QHimhihi[j][i] = Qimhihi[i][j]; QHimlohi[j][i] = Qimlohi[i][j];
         QHimhilo[j][i] = Qimhilo[i][j]; QHimlolo[j][i] = Qimlolo[i][j];

         qdf_minus(&QHimhihi[j][i],&QHimlohi[j][i],
                   &QHimhilo[j][i],&QHimlolo[j][i]); // Hermitian transpose
      }

   CPU_cmplx4_factors_matmatmul
      (nrows,nrows,nrows,QHrehihi, QHrelohi,  QHrehilo, QHrelolo,
                         QHimhihi, QHimlohi,  QHimhilo, QHimlolo,
                          Qrehihi,  Qrelohi,  Qrehilo,  Qrelolo, 
                          Qimhihi,  Qimlohi,  Qimhilo,  Qimlolo,
                        QHQrehihi,QHQrelohi,QHQrehilo,QHQrelolo,
                        QHQimhihi,QHQimlohi,QHQimhilo,QHQimlolo);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
         {
            cout << "Q^H*Q[" << i << "][" << j << "]re : "
                 << QHQrehihi[i][j] << "  " << QHQrelohi[i][j] << endl
                 << "                "
                 << QHQrehilo[i][j] << "  " << QHQrelolo[i][j] << endl;
            cout << "Q^H*Q[" << i << "][" << j << "]im : "
                 << QHQimhihi[i][j] << "  " << QHQimlohi[i][j] << endl
                 << "                "
                 << QHQimhilo[i][j] << "  " << QHQimlolo[i][j] << endl;
         }
         if(i == j)
            errorQ = errorQ + abs(QHQrehihi[i][j] - 1.0)
                            + abs(QHQrelohi[i][j])
                            + abs(QHQrehilo[i][j]) + abs(QHQrelolo[i][j])
                            + abs(QHQimhihi[i][j]) + abs(QHQimlohi[i][j])
                            + abs(QHQimhilo[i][j]) + abs(QHQimlolo[i][j]);
         else
            errorQ = errorQ + abs(QHQrehihi[i][j]) + abs(QHQrelohi[i][j])
                            + abs(QHQrehilo[i][j]) + abs(QHQrelolo[i][j])
                            + abs(QHQimhihi[i][j]) + abs(QHQimlohi[i][j])
                            + abs(QHQimhilo[i][j]) + abs(QHQimlolo[i][j]);
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
                 << Rrehihi[i][j] << "  " << Rrelohi[i][j] << endl
                 << "            "
                 << Rrehilo[i][j] << "  " << Rrelolo[i][j] << endl;
            cout << "R[" << i << "][" << j << "]im : "
                 << Rimhihi[i][j] << "  " << Rimlohi[i][j] << endl
                 << "            "
                 << Rimhilo[i][j] << "  " << Rimlolo[i][j] << endl;
         }
   }
   CPU_cmplx4_factors_matmatmul
      (nrows,nrows,ncols, QHrehihi, QHrelohi, QHrehilo, QHrelolo,
                          QHimhihi, QHimlohi, QHimhilo, QHimlolo,
                           Arehihi,  Arelohi,  Arehilo,  Arelolo, 
                           Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
                         QHArehihi,QHArelohi,QHArehilo,QHArelolo,
                         QHAimhihi,QHAimlohi,QHAimhilo,QHAimlolo);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
         {
            cout << "Q^H*A[" << i << "][" << j << "]re : "
                 << QHArehihi[i][j] << "  " << QHArelohi[i][j] << endl
                 << "                "
                 << QHArehilo[i][j] << "  " << QHArelolo[i][j] << endl;
            cout << "Q^H*A[" << i << "][" << j << "]im : "
                 << QHAimhihi[i][j] << "  " << QHAimlohi[i][j] << endl
                 << "                "
                 << QHAimhilo[i][j] << "  " << QHAimlolo[i][j] << endl;
         }
         errorR = errorR + abs(Rrehihi[i][j] - QHArehihi[i][j])
                         + abs(Rrelohi[i][j] - QHArelohi[i][j])
                         + abs(Rrehilo[i][j] - QHArehilo[i][j])
                         + abs(Rrelolo[i][j] - QHArelolo[i][j])
                         + abs(Rimhihi[i][j] - QHAimhihi[i][j])
                         + abs(Rimlohi[i][j] - QHAimlohi[i][j])
                         + abs(Rimhilo[i][j] - QHAimhilo[i][j])
                         + abs(Rimlolo[i][j] - QHAimlolo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QHrehihi[i]); free(QHimhihi[i]);
      free(QHrelohi[i]); free(QHimlohi[i]);
      free(QHrehilo[i]); free(QHimhilo[i]);
      free(QHrelolo[i]); free(QHimlolo[i]);
      free(QHQrehihi[i]); free(QHQimhihi[i]);
      free(QHQrelohi[i]); free(QHQimlohi[i]);
      free(QHQrehilo[i]); free(QHQimhilo[i]);
      free(QHQrelolo[i]); free(QHQimlolo[i]);
      free(QHArehihi[i]); free(QHAimhihi[i]);
      free(QHArelohi[i]); free(QHAimlohi[i]);
      free(QHArehilo[i]); free(QHAimhilo[i]);
      free(QHArelolo[i]); free(QHAimlolo[i]);
   }
   free(QHrehihi); free(QHQrehihi); free(QHArehihi);
   free(QHrelohi); free(QHQrelohi); free(QHArelohi);
   free(QHrehilo); free(QHQrehilo); free(QHArehilo);
   free(QHrelolo); free(QHQrelolo); free(QHArelolo);
   free(QHimhihi); free(QHQimhihi); free(QHAimhihi);
   free(QHimlohi); free(QHQimlohi); free(QHAimlohi);
   free(QHimhilo); free(QHQimhilo); free(QHAimhilo);
   free(QHimlolo); free(QHQimlolo); free(QHAimlolo);

   return int(errorQ + errorR > tol);
}

void test_factors_real4_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Ahihi = new double*[nrows];
   double **Alohi = new double*[nrows];
   double **Ahilo = new double*[nrows];
   double **Alolo = new double*[nrows];
   double **Qhihi = new double*[nrows];
   double **Qlohi = new double*[nrows];
   double **Qhilo = new double*[nrows];
   double **Qlolo = new double*[nrows];
   double **Rhihi = new double*[nrows];
   double **Rlohi = new double*[nrows];
   double **Rhilo = new double*[nrows];
   double **Rlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihi[i] = new double[ncols];
      Alohi[i] = new double[ncols];
      Ahilo[i] = new double[ncols];
      Alolo[i] = new double[ncols];
      Qhihi[i] = new double[nrows];
      Qlohi[i] = new double[nrows];
      Qhilo[i] = new double[nrows];
      Qlolo[i] = new double[nrows];
      Rhihi[i] = new double[ncols];
      Rlohi[i] = new double[ncols];
      Rhilo[i] = new double[ncols];
      Rlolo[i] = new double[ncols];
   }
   random_dbl4_matrix(nrows,ncols,Ahihi,Alohi,Ahilo,Alolo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi[i][j] << "  " << Alohi[i][j] << endl
                 << "          "
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }
   CPU_dbl4_factors_houseqr
      (nrows,ncols,Ahihi,Alohi,Ahilo,Alolo,
                   Qhihi,Qlohi,Qhilo,Qlolo,
                   Rhihi,Rlohi,Rhilo,Rlolo);

   const double tol = 1.0E-26;
   const int fail = test_real4_qr_factors
      (nrows,ncols,Ahihi,Alohi,Ahilo,Alolo,
                   Qhihi,Qlohi,Qhilo,Qlolo,
                   Rhihi,Rlohi,Rhilo,Rlolo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Ahihi[i]); free(Qhihi[i]); free(Rhihi[i]);
      free(Alohi[i]); free(Qlohi[i]); free(Rlohi[i]);
      free(Ahilo[i]); free(Qhilo[i]); free(Rhilo[i]);
      free(Alolo[i]); free(Qlolo[i]); free(Rlolo[i]);
   }
   free(Ahihi); free(Qhihi); free(Rhihi);
   free(Alohi); free(Qlohi); free(Rlohi);
   free(Ahilo); free(Qhilo); free(Rhilo);
   free(Alolo); free(Qlolo); free(Rlolo);
}

void test_factors_cmplx4_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Arehihi = new double*[nrows];
   double **Arelohi = new double*[nrows];
   double **Arehilo = new double*[nrows];
   double **Arelolo = new double*[nrows];
   double **Aimhihi = new double*[nrows];
   double **Aimlohi = new double*[nrows];
   double **Aimhilo = new double*[nrows];
   double **Aimlolo = new double*[nrows];
   double **Qrehihi = new double*[nrows];
   double **Qrelohi = new double*[nrows];
   double **Qrehilo = new double*[nrows];
   double **Qrelolo = new double*[nrows];
   double **Qimhihi = new double*[nrows];
   double **Qimlohi = new double*[nrows];
   double **Qimhilo = new double*[nrows];
   double **Qimlolo = new double*[nrows];
   double **Rrehihi = new double*[nrows];
   double **Rrelohi = new double*[nrows];
   double **Rrehilo = new double*[nrows];
   double **Rrelolo = new double*[nrows];
   double **Rimhihi = new double*[nrows];
   double **Rimlohi = new double*[nrows];
   double **Rimhilo = new double*[nrows];
   double **Rimlolo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehihi[i] = new double[ncols];
      Arelohi[i] = new double[ncols];
      Arehilo[i] = new double[ncols];
      Arelolo[i] = new double[ncols];
      Aimhihi[i] = new double[ncols];
      Aimlohi[i] = new double[ncols];
      Aimhilo[i] = new double[ncols];
      Aimlolo[i] = new double[ncols];
      Qrehihi[i] = new double[nrows];
      Qrelohi[i] = new double[nrows];
      Qrehilo[i] = new double[nrows];
      Qrelolo[i] = new double[nrows];
      Qimhihi[i] = new double[nrows];
      Qimlohi[i] = new double[nrows];
      Qimhilo[i] = new double[nrows];
      Qimlolo[i] = new double[nrows];
      Rrehihi[i] = new double[ncols];
      Rrelohi[i] = new double[ncols];
      Rrehilo[i] = new double[ncols];
      Rrelolo[i] = new double[ncols];
      Rimhihi[i] = new double[ncols];
      Rimlohi[i] = new double[ncols];
      Rimhilo[i] = new double[ncols];
      Rimlolo[i] = new double[ncols];
   }
   random_cmplx4_matrix
      (nrows,ncols,Arehihi,Arelohi,Arehilo,Arelolo,
                   Aimhihi,Aimlohi,Aimhilo,Aimlolo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << "            "
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << "            "
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   CPU_cmplx4_factors_houseqr
      (nrows,ncols,Arehihi,Arelohi,Arehilo,Arelolo,
                   Aimhihi,Aimlohi,Aimhilo,Aimlolo,
                   Qrehihi,Qrelohi,Qrehilo,Qrelolo,
                   Qimhihi,Qimlohi,Qimhilo,Qimlolo,
                   Rrehihi,Rrelohi,Rrehilo,Rrelolo,
                   Rimhihi,Rimlohi,Rimhilo,Rimlolo);

   const double tol = 1.0e-12;
   const int fail = test_cmplx4_qr_factors
      (nrows,ncols,Arehihi,Arelohi,Arehilo,Arelolo,
                   Aimhihi,Aimlohi,Aimhilo,Aimlolo,
                   Qrehihi,Qrelohi,Qrehilo,Qrelolo,
                   Qimhihi,Qimlohi,Qimhilo,Qimlolo,
                   Rrehihi,Rrelohi,Rrehilo,Rrelolo,
                   Rimhihi,Rimlohi,Rimhilo,Rimlolo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Arehihi[i]); free(Qrehihi[i]); free(Rrehihi[i]);
      free(Arelohi[i]); free(Qrelohi[i]); free(Rrelohi[i]);
      free(Arehilo[i]); free(Qrehilo[i]); free(Rrehilo[i]);
      free(Arelolo[i]); free(Qrelolo[i]); free(Rrelolo[i]);
      free(Aimhihi[i]); free(Qimhihi[i]); free(Rimhihi[i]);
      free(Aimlohi[i]); free(Qimlohi[i]); free(Rimlohi[i]);
      free(Aimhilo[i]); free(Qimhilo[i]); free(Rimhilo[i]);
      free(Aimlolo[i]); free(Qimlolo[i]); free(Rimlolo[i]);
   }
   free(Arehihi); free(Qrehihi); free(Rrehihi);
   free(Arelohi); free(Qrelohi); free(Rrelohi);
   free(Arehilo); free(Qrehilo); free(Rrehilo);
   free(Arelolo); free(Qrelolo); free(Rrelolo);
   free(Aimhihi); free(Qimhihi); free(Rimhihi);
   free(Aimlohi); free(Qimlohi); free(Rimlohi);
   free(Aimhilo); free(Qimhilo); free(Rimhilo);
   free(Aimlolo); free(Qimlolo); free(Rimlolo);
}
