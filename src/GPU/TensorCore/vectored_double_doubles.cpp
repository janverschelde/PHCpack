/* Collection of functions for vectored double double arithmetic. */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "double_double.h"
#include "double_double_functions.h"
#include "splitting_doubles.h"

using namespace std;

bool is_dd_quarter_balanced
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3, int vrblvl=0 )
{
   if(vrblvl > 0)
      cout << "-> in vectored_double_doubles.is_dd_quarter_balanced ..."
           << endl;

   bool b01 = is_quarter_balanced(xhi0, xhi1, vrblvl-1);
   bool b12 = is_quarter_balanced(xhi1, xhi2, vrblvl-1);
   bool b23 = is_quarter_balanced(xhi2, xhi3, vrblvl-1);
   bool b34 = is_quarter_balanced(xhi3, xlo0, vrblvl-1);
   bool b45 = is_quarter_balanced(xlo0, xlo1, vrblvl-1);
   bool b56 = is_quarter_balanced(xlo1, xlo2, vrblvl-1);
   bool b67 = is_quarter_balanced(xlo2, xlo3, vrblvl-1);

   bool result = b01 and b12 and b23 and b34 and b45 and b56 and b67;
   if(vrblvl > 0)
   {
      if(result)
         cout << "All eight quarters of the double double are balanced."
              << endl;
      else
         cout << "Not all eight quarters of the double double are balanced."
              << endl;
   }
   return result;
}

void dd_balance_quarters
 ( double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3, int vrblvl=0 )
{
   if(vrblvl > 0)
      cout << "-> in vectored_double_doubles.dd_balance_quarters ..." << endl;

   if(not is_quarter_balanced(*xhi0, *xhi1, vrblvl-1))
      quarter_balance(xhi0, xhi1, vrblvl-1);
   if(not is_quarter_balanced(*xhi1, *xhi2, vrblvl-1))
      quarter_balance(xhi1, xhi2, vrblvl-1);
   if(not is_quarter_balanced(*xhi2, *xhi3, vrblvl-1))
      quarter_balance(xhi2, xhi3, vrblvl-1);
   if(not is_quarter_balanced(*xhi3, *xlo0, vrblvl-1))
      quarter_balance(xhi3, xlo0, vrblvl-1);
   if(not is_quarter_balanced(*xlo0, *xlo1, vrblvl-1))
      quarter_balance(xlo0, xlo1, vrblvl-1);
   if(not is_quarter_balanced(*xlo1, *xlo2, vrblvl-1))
      quarter_balance(xlo1, xlo2, vrblvl-1);
   if(not is_quarter_balanced(*xlo2, *xlo3, vrblvl-1))
      quarter_balance(xlo2, xlo3, vrblvl-1);
}

void make_dd_exponent_zero ( double *xhi, double *xlo, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in vectored_double_doubles.make_dd_exponent_zero ..."
           << endl;

   double factor;

   make_exponent_zero(xhi, &factor, vrblvl-1);
   *xlo = (*xlo)*factor;
}

void quarter_double_double
 ( double xhi, double xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in vectored_double_doubles.quarter_double_double ..."
           << endl;

   quarter_split(xhi, xhi0, xhi1, xhi2, xhi3, vrblvl-1);
   if(vrblvl > 0)
   {
      if(*xhi0 == 0.0) cout << "xhi0 is zero!" << endl;
      if(*xhi1 == 0.0) cout << "xhi1 is zero!" << endl;
      if(*xhi2 == 0.0) cout << "xhi2 is zero!" << endl;
      if(*xhi3 == 0.0) cout << "xhi3 is zero!" << endl;
   }
   quarter_split(xlo, xlo0, xlo1, xlo2, xlo3, vrblvl-1);
   if(vrblvl > 0)
   {
      if(*xlo0 == 0.0) cout << "xlo0 is zero!" << endl;
      if(*xlo1 == 0.0) cout << "xlo1 is zero!" << endl;
      if(*xlo2 == 0.0) cout << "xlo2 is zero!" << endl;
      if(*xlo3 == 0.0) cout << "xlo3 is zero!" << endl;
   }
   dd_balance_quarters
      (xhi0, xhi1, xhi2, xhi3, xlo0, xlo1, xlo2, xlo3, vrblvl-1);
}

void quarter_dd_vector
 ( int dim, double *xhi, double *xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in vectored_double_doubles.quarter_dd_vector ..."
           << endl;

   for(int i=0; i<dim; i++)
   {
      quarter_double_double
         (xhi[i], xlo[i],
          &xhi0[i], &xhi1[i], &xhi2[i], &xhi3[i],
          &xlo0[i], &xlo1[i], &xlo2[i], &xlo3[i], vrblvl-1);
   }
   if(vrblvl > 0)
   {
      bool check = true;
      bool fail = false;

      for(int i=0; i<dim; i++)
      {
         check = is_dd_quarter_balanced
                    (xhi0[i], xhi1[i], xhi2[i], xhi3[i],
                     xlo0[i], xlo1[i], xlo2[i], xlo3[i], vrblvl-1);
         if(not check)
         {
            cout << i << "-th quarters are not balanced!" << endl;
            fail = true;
            // redo the test for the output
            check = is_dd_quarter_balanced
                      (xhi0[i], xhi1[i], xhi2[i], xhi3[i],
                       xlo0[i], xlo1[i], xlo2[i], xlo3[i], 2);
         }
      }
      if(fail)
         cout << "Not all quarters in the vector are balanced." << endl;
      else
         cout << "All quarters in the vector are balanced." << endl;
   }
}

void quarter_dd_matrix
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in vectored_double_doubles.quarter_dd_matrix ..."
           << endl;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         quarter_double_double
            (Ahi[i][j], Alo[i][j],
             &Ahi0[i][j], &Ahi1[i][j], &Ahi2[i][j], &Ahi3[i][j],
             &Alo0[i][j], &Alo1[i][j], &Alo2[i][j], &Alo3[i][j], vrblvl-1);
      }
}

void to_double_double
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3,
   double *xhi, double *xlo )
{
/*
   *xhi = xhi0;
   *xlo = 0.0;

   ddf_inc_d(xhi, xlo, xhi1);
   ddf_inc_d(xhi, xlo, xhi2);
   ddf_inc_d(xhi, xlo, xhi3);
   ddf_inc_d(xhi, xlo, xlo0);
   ddf_inc_d(xhi, xlo, xlo1);
   ddf_inc_d(xhi, xlo, xlo2);
   ddf_inc_d(xhi, xlo, xlo3);
*/
   double z0,z1,z2,z3,e1,e2,e3,e4;

   z0 = ddf_two_sum(xhi0, xhi1, &e1);
   z1 = ddf_two_sum(xhi2, xhi3, &e2);
   z2 = ddf_two_sum(xlo0, xlo1, &e3);
   z3 = ddf_two_sum(xlo2, xlo3, &e4);

   *xhi = z3; *xlo = e4;
   ddf_inc(xhi, xlo, z2, e3);
   ddf_inc(xhi, xlo, z1, e2);
   ddf_inc(xhi, xlo, z0, e1);
}

void to_double_double_matrix
 ( int nrows, int ncols,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3,
   double **Ahi, double **Alo )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         to_double_double
            (Ahi0[i][j], Ahi1[i][j], Ahi2[i][j], Ahi3[i][j],
             Alo0[i][j], Alo1[i][j], Alo2[i][j], Alo3[i][j],
             &Ahi[i][j], &Alo[i][j]);
      }
}

void dd_write_vector ( int dim, double *xhi, double *xlo )
{
   double x[2];

   for(int i=0; i<dim; i++)
   {
      x[0] = xhi[i]; x[1] = xlo[i];
      dd_write(x, 32); cout << endl;
   }
}

void double_double_product
 ( int dim, double *xhi, double *xlo, double *yhi, double *ylo,
   double *prdhi, double *prdlo )
{
   double acchi,acclo;

   *prdhi = 0.0;
   *prdlo = 0.0;

   for(int i=0; i<dim; i++)
   {
      ddf_mul(xhi[i], xlo[i], yhi[i], ylo[i], &acchi, &acclo);
      ddf_inc(prdhi, prdlo, acchi, acclo);
   }
}

void double_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo )
{
   double acchi,acclo;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chi[i][j] = 0.0; Clo[i][j] = 0.0;

         for(int k=0; k<dim; k++)
         {
            ddf_mul(Ahi[i][k], Alo[i][k], Bhi[k][j], Blo[k][j],
                    &acchi, &acclo);
            ddf_inc(&Chi[i][j], &Clo[i][j], acchi, acclo);
         }
      }
}

void vectored_dd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7 )
{
   *s0 = 0.0; *s1 = 0.0; *s2 = 0.0; *s3 = 0.0;
   *s4 = 0.0; *s5 = 0.0; *s6 = 0.0; *s7 = 0.0;

   cout << scientific << setprecision(16);

   for(int i=0; i<dim; i++)
   {
      *s0 += x0[i]*y0[i];
      *s1 += x0[i]*y1[i] + x1[i]*y0[i];
      *s2 += x0[i]*y2[i] + x1[i]*y1[i] + x2[i]*y0[i];
      *s3 += x0[i]*y3[i] + x1[i]*y2[i] + x2[i]*y1[i] + x3[i]*y0[i];
      *s4 += x0[i]*y4[i] + x1[i]*y3[i] + x2[i]*y2[i] + x3[i]*y1[i]
           + x4[i]*y0[i];
      *s5 += x0[i]*y5[i] + x1[i]*y4[i] + x2[i]*y3[i] + x3[i]*y2[i]
           + x4[i]*y1[i] + x5[i]*y0[i];
      *s6 += x0[i]*y6[i] + x1[i]*y5[i] + x2[i]*y4[i] + x3[i]*y3[i]
           + x4[i]*y2[i] + x5[i]*y1[i] + x6[i]*y0[i];
      *s7 += x0[i]*y7[i] + x1[i]*y6[i] + x2[i]*y5[i] + x3[i]*y4[i]
           + x4[i]*y3[i] + x5[i]*y2[i] + x6[i]*y1[i] + x7[i]*y0[i];
   }
}

void transpose_dd_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         T0[j][i] = A0[i][j]; T1[j][i] = A1[i][j];
         T2[j][i] = A2[i][j]; T3[j][i] = A3[i][j];
         T4[j][i] = A4[i][j]; T5[j][i] = A5[i][j];
         T6[j][i] = A6[i][j]; T7[j][i] = A7[i][j];
      }
}

void vectored_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7 )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         vectored_dd_product
            (dim, A0[i], A1[i], A2[i], A3[i], A4[i], A5[i], A6[i], A7[i],
                  B0[j], B1[j], B2[j], B3[j], B4[j], B5[j], B6[j], B7[j],
             &C0[i][j], &C1[i][j], &C2[i][j], &C3[i][j],
             &C4[i][j], &C5[i][j], &C6[i][j], &C7[i][j]);
      }
}

void dd_convolute_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7, double **cA )
{
   for(int i=0; i<8*nrows; i++)
      for(int j=0; j<8*ncols; j++) cA[i][j] = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         for(int k=0; k<8; k++) cA[8*i+k][8*j+k] = A0[i][j];
         for(int k=0; k<7; k++) cA[8*i+k+1][8*j+k] = A1[i][j];
         for(int k=0; k<6; k++) cA[8*i+k+2][8*j+k] = A2[i][j];
         for(int k=0; k<5; k++) cA[8*i+k+3][8*j+k] = A3[i][j];
         for(int k=0; k<4; k++) cA[8*i+k+4][8*j+k] = A4[i][j];
         for(int k=0; k<3; k++) cA[8*i+k+5][8*j+k] = A5[i][j];
         for(int k=0; k<2; k++) cA[8*i+k+6][8*j+k] = A6[i][j];
         cA[8*i+7][8*j] = A7[i][j];
      }
}

void dd_stack_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7, double **sA )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         sA[8*i][j] = A0[i][j];
         sA[8*i+1][j] = A1[i][j];
         sA[8*i+2][j] = A2[i][j];
         sA[8*i+3][j] = A3[i][j];
         sA[8*i+4][j] = A4[i][j];
         sA[8*i+5][j] = A5[i][j];
         sA[8*i+6][j] = A6[i][j];
         sA[8*i+7][j] = A7[i][j];
      }
}

void extract_dd_quarters
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7 )
{
   for(int i=0; i<8*nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         int k = i % 8;
         int row = i/8;

         if(k == 0) D0[row][j] = qC[i][j];
         if(k == 1) D1[row][j] = qC[i][j];
         if(k == 2) D2[row][j] = qC[i][j];
         if(k == 3) D3[row][j] = qC[i][j];
         if(k == 4) D4[row][j] = qC[i][j];
         if(k == 5) D5[row][j] = qC[i][j];
         if(k == 6) D6[row][j] = qC[i][j];
         if(k == 7) D7[row][j] = qC[i][j];
      }
}
