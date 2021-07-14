// The file dbl2_baqr_testers.cpp defines the function with prototypes in
// the file dbl2_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random2_matrices.h"
#include "dbl2_factorizations.h"
#include "dbl2_factors_testers.h"
#include "dbl2_baqr_host.h"

using namespace std;

void test_real2_blocked_qr ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   const int ncols = sizetile*numtiles;

   cout << "Give the number of rows (" << " >= " << ncols << " ) : ";
   int nrows; cin >> nrows;

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
   bool vrb = (verbose > 0);
   CPU_dbl2_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,vrb);

   test_real2_qr_factors(nrows,ncols,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,verbose);

   for(int i=0; i<nrows; i++)
   {
      free(Ahi[i]); free(Qhi[i]); free(Rhi[i]);
      free(Alo[i]); free(Qlo[i]); free(Rlo[i]);
   }
   free(Ahi); free(Qhi); free(Rhi);
   free(Alo); free(Qlo); free(Rlo);
}

void test_cmplx2_blocked_qr ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   const int ncols = sizetile*numtiles;

   cout << "Give the number of rows (" << " >= " << ncols << " ) : ";
   int nrows; cin >> nrows;

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

   bool vrb = (verbose > 0);
   CPU_cmplx2_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,
       Arehi,Arelo,Aimhi,Aimlo,
       Qrehi,Qrelo,Qimhi,Qimlo,
       Rrehi,Rrelo,Rimhi,Rimlo,vrb);

   test_cmplx2_qr_factors
      (nrows,ncols,Arehi,Arelo,Aimhi,Aimlo,
                   Qrehi,Qrelo,Qimhi,Qimlo,
                   Rrehi,Rrelo,Rimhi,Rimlo,verbose);

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
