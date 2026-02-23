/* The file smDMMA_host.cpp contains the definitions of the functions with
 * prototypes in smDMMA_host.h. */

#include "smDMMA_dims.h"
#include "smDMMA_host.h"
#include <iostream>
#include "double_matrix_multiplications.h"
#include "random2_matrices.h"
#include "vectored_double_doubles.h"

void init_host_matrices ( double *a, double *b, double *c, int nbrange )
{
   for(int i=0; i<K_GLOBAL; i++) // #rows: M_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<i; j++) a[i*K_GLOBAL+j] = 0.0;
      for(int j=i; j<K_GLOBAL; j++) a[i*K_GLOBAL+j] = 1.0;
   }
   for(int i=K_GLOBAL; i<M_GLOBAL; i++) // row major with M_GLOBAL rows
      for(int j=0; j<K_GLOBAL; j++) a[i*K_GLOBAL+j] = 0.0;

   for(int i=0; i<K_GLOBAL; i++) // #columns: N_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<=i; j++) b[i*K_GLOBAL+j] = 1.0;
      for(int j=i+1; j<K_GLOBAL; j++) b[i*K_GLOBAL+j] = 0.0;
   }
   for(int i=K_GLOBAL; i<N_GLOBAL; i++) // #columns: N_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<K_GLOBAL; j++) 
      {
         if(j % nbrange == 0)
            b[i*K_GLOBAL+j] = 1.0;
         else
            b[i*K_GLOBAL+j] = 2.0*b[i*K_GLOBAL+j-1];
      }
   }

   for(int t=0; t<M_GLOBAL*N_GLOBAL; t++) c[t] = 0.0;
}

void random_double_double_matrices
 ( double *A, double *B, double *C,
   int numArows, int numAcols, int numBrows, int numBcols,
   int numCrows, int numCcols, int vrblvl )
{
   using namespace std;

   if(vrblvl > 0)
   {
      cout << "-> in smDMMA_host.random_double_double_matrices ..." << endl;
      cout << "#rows in A : " << numArows
           << ", #columns in A : " << numAcols << endl;
      cout << "#rows in B : " << numBrows
           << ", #columns in B : " << numBcols << endl;
      cout << "#rows in C : " << numCrows
           << ", #columns in C : " << numCcols << endl;
   }
   const int nddArows = numArows/12;
   const int nddAcols = numAcols/12;
   const int nddBrows = numBrows/12;
   const int nddBcols = numBcols;
   const int nddCrows = numCrows/12;
   const int nddCcols = numCcols;

   if(vrblvl > 0)
   {
      cout << "#rows in dd A : " << nddArows
           << ", #columns in dd A : " << nddAcols << endl;
      cout << "#rows in dd B : " << nddBrows
           << ", #columns in dd B : " << nddBcols << endl;
      cout << "#rows in dd C : " << nddCrows
           << ", #columns in dd C : " << nddCcols << endl;
   }

   double **Ahi = new double*[nddArows];
   double **Alo = new double*[nddArows];
   double **Bhi = new double*[nddBrows];
   double **Blo = new double*[nddBrows];

   if(vrblvl > 0)
      cout << "-> generating a random " << nddArows
           << "-by-" << nddAcols << " double double matrix A ..." << endl;

   for(int i=0; i<nddArows; i++)
   {
      Ahi[i] = new double[nddAcols];
      Alo[i] = new double[nddAcols];
   }
   random_dbl2_matrix(nddArows, nddAcols, Ahi, Alo);

   for(int i=0; i<nddArows; i++)
      for(int j=0; j<nddAcols; j++)
      {
         if(Ahi[i][j] < 0.0) Ahi[i][j] = -Ahi[i][j];
         if(Alo[i][j] < 0.0) Alo[i][j] = -Alo[i][j];
         make_dd_exponent_zero(&Ahi[i][j], &Alo[i][j]);
      }

   if(vrblvl > 0)
      cout << "-> generating a random " << nddBrows
           << "-by-" << nddBcols << " double double matrix B ..." << endl;

   for(int i=0; i<nddBrows; i++)
   {
      Bhi[i] = new double[nddBcols];
      Blo[i] = new double[nddBcols];
   }
   random_dbl2_matrix(nddBrows, nddBcols, Bhi, Blo);

   for(int i=0; i<nddBrows; i++)
      for(int j=0; j<nddBcols; j++)
      {
         if(Bhi[i][j] < 0.0) Bhi[i][j] = -Bhi[i][j];
         if(Blo[i][j] < 0.0) Blo[i][j] = -Blo[i][j];
         make_dd_exponent_zero(&Bhi[i][j], &Blo[i][j]);
      }

   double **Ahi0 = new double*[nddArows];
   double **Ahi1 = new double*[nddArows];
   double **Ahi2 = new double*[nddArows];
   double **Ahi3 = new double*[nddArows];
   double **Alo0 = new double*[nddArows];
   double **Alo1 = new double*[nddArows];
   double **Alo2 = new double*[nddArows];
   double **Alo3 = new double*[nddArows];
   double **Alo4 = new double*[nddArows];
   double **Alo5 = new double*[nddArows];
   double **Alo6 = new double*[nddArows];
   double **Alo7 = new double*[nddArows];

   for(int i=0; i<nddArows; i++)
   {
      Ahi0[i] = new double[nddAcols];
      Ahi1[i] = new double[nddAcols];
      Ahi2[i] = new double[nddAcols];
      Ahi3[i] = new double[nddAcols];
      Alo0[i] = new double[nddAcols];
      Alo1[i] = new double[nddAcols];
      Alo2[i] = new double[nddAcols];
      Alo3[i] = new double[nddAcols];
      Alo4[i] = new double[nddAcols];
      Alo5[i] = new double[nddAcols];
      Alo6[i] = new double[nddAcols];
      Alo7[i] = new double[nddAcols];
   }
   if(vrblvl > 0)
      cout << "-> 12-splitting a " << nddArows
           << "-by-" << nddAcols << " double double matrix A ..." << endl;

   split_dd_matrix
      (nddArows, nddAcols, Ahi, Alo, Ahi0, Ahi1, Ahi2, Ahi3,
       Alo0, Alo1, Alo2, Alo3, Alo4, Alo5, Alo6, Alo7);

   double **Bhi0 = new double*[nddBrows];
   double **Bhi1 = new double*[nddBrows];
   double **Bhi2 = new double*[nddBrows];
   double **Bhi3 = new double*[nddBrows];
   double **Blo0 = new double*[nddBrows];
   double **Blo1 = new double*[nddBrows];
   double **Blo2 = new double*[nddBrows];
   double **Blo3 = new double*[nddBrows];
   double **Blo4 = new double*[nddBrows];
   double **Blo5 = new double*[nddBrows];
   double **Blo6 = new double*[nddBrows];
   double **Blo7 = new double*[nddBrows];

   for(int i=0; i<nddBrows; i++)
   {
      Bhi0[i] = new double[nddBcols];
      Bhi1[i] = new double[nddBcols];
      Bhi2[i] = new double[nddBcols];
      Bhi3[i] = new double[nddBcols];
      Blo0[i] = new double[nddBcols];
      Blo1[i] = new double[nddBcols];
      Blo2[i] = new double[nddBcols];
      Blo3[i] = new double[nddBcols];
      Blo4[i] = new double[nddBcols];
      Blo5[i] = new double[nddBcols];
      Blo6[i] = new double[nddBcols];
      Blo7[i] = new double[nddBcols];
   }
   if(vrblvl > 0)
      cout << "-> 12-splitting a " << nddBrows
           << "-by-" << nddBcols << " double double matrix B ..." << endl;

   split_dd_matrix
      (nddBrows, nddBcols, Bhi, Blo, Bhi0, Bhi1, Bhi2, Bhi3,
       Blo0, Blo1, Blo2, Blo3, Blo4, Blo5, Blo6, Blo7);

   if(vrblvl > 0)
      cout << "-> convoluting the 12-splitted matrix A ..." << endl;

   double **cA = new double*[12*nddArows];
   for(int i=0; i<12*nddArows; i++) cA[i] = new double[12*nddAcols];

   dd_convolute_12splits
      (nddArows, nddAcols, Ahi0, Ahi1, Ahi2, Ahi3,
       Alo0, Alo1, Alo2, Alo3, Alo4, Alo5, Alo6, Alo7, cA);

   if(vrblvl > 0)
      cout << "-> stacking the 12-splitted matrix B ..." << endl;

   double **sB = new double*[12*nddBrows];
   for(int i=0; i<12*nddBrows; i++) sB[i] = new double[nddBcols];

   dd_stack_12splits
      (nddBrows, nddBcols, Bhi0, Bhi1, Bhi2, Bhi3,
       Blo0, Blo1, Blo2, Blo3, Blo4, Blo5, Blo6, Blo7, sB);

   const int padArows = numArows - 12*nddArows;
   const int padAcols = numAcols - 12*nddAcols;
   
   if(vrblvl > 0)
   {
      cout << "  #rows A : " << numArows
           << ",   12 x #rows dd A : " << 12*nddArows
           << " => pad rows A : " << padArows << endl;
      cout << "#colums A : " << numAcols
           << ", 12 x #colums dd A : " << 12*nddAcols
           << " => pad columns A : " << padAcols << endl;
      cout << "-> padding the convoluted 12-splitted matrix A ..." << endl;
   }
   double **pcA = new double*[numArows];
   for(int i=0; i<numArows; i++)
   {
      pcA[i] = new double[numAcols];
      for(int j=0; j<numAcols; j++) pcA[i][j] = 0.0;
   }
   for(int i=0; i<12*nddArows; i++)
      for(int j=0; j<12*nddAcols; j++) pcA[i][j] = cA[i][j];

   if(vrblvl > 0)
      cout << "-> single indexing the convoluted 12-splitted matrix A ..."
           << endl;

   double2single_row_major(numArows, numAcols, pcA, A);

   const int padBrows = numBrows - 12*nddBrows;
   const int padBcols = numBcols - nddBcols;

   if(vrblvl > 0)
   {
      cout << "  #rows B : " << numBrows
           << ",   12 x #rows dd B : " << 12*nddBrows
           << " => pad rows B : " << padBrows << endl;
      cout << "#colums B : " << numBcols
           << ", 12 x #colums dd B : " << 12*nddBcols
           << " => pad columns B : " << padBcols << endl;
      cout << "-> padding the stacked 12-splitted matrix B ..." << endl;
   }
   double **psB = new double*[numBrows];
   for(int i=0; i<numBrows; i++)
   {
      psB[i] = new double[numBcols];
      for(int j=0; j<numBcols; j++) psB[i][j] = 0.0;
   }
   for(int i=0; i<12*nddBrows; i++)
      for(int j=0; j<nddBcols; j++) psB[i][j] = sB[i][j];

   if(vrblvl > 0)
      cout << "-> transposing the padded stacked 12-splitted matrix"
           << " B ..." << endl; 

   double **sBT = new double*[numBcols];
   for(int i=0; i<numBcols; i++) sBT[i] = new double[numBrows];
   transpose_rows_columns(numBrows, numBcols, psB, sBT);

   if(vrblvl > 0)
      cout << "-> single indexing the transposed padded stacked 12-splitted"
           << " matrix B ..." << endl;
  
   double2single_column_major(numBrows, numBcols, sBT, B);

   if(vrblvl > 0)
      cout << "-> initializing the matrix C ..." << endl;

   for(int k=0; k<numCrows*numCcols; k++) C[k] = 0.0;
}


void matMultiplyOnHost
 ( double *A, double *B, double *C, float alpha, float beta,
   int numARows, int numAColumns, int numBRows, int numBColumns,
   int numCRows, int numCColumns )
{
   for(int i = 0; i < numCRows; i++)
   {
      for(int j = 0; j < numCColumns; j++)
      {
         double temp = 0.0;

         for(int k = 0; k < numAColumns; k++)
         {
            // B matrix is column major. A matrix is row major.
            temp += A[i*numAColumns+k] * B[j*numBRows+k];
         }
         C[i*numCColumns+j] = temp*alpha + beta*C[i*numCColumns+j];
        }
    }
}
