// The file unimodular_matrices.cpp defines functions specified
// in unimodular_matrices.h.

#include <iostream>
#include <cstdlib>
#include "random_monomials.h"
#include "unimodular_matrices.h"

using namespace std;

void random_row_and_column ( int dim, int *row, int *col )
{
   *row = rand() % dim;
   *col = *row;

   while(*row == *col) *col = rand() % dim;
}

int nonzero_integer ( int size, bool poscff )
{
   int result = 0;

   while(result == 0)
   {
      if(poscff)
         result = rand() % size;
      else
      {
         result = rand() % 2*size+1;
         result = result - size;
      }
   }
   return result;
}

void unimodular_multiplication
 ( int dim, int row, int col, int nbr, int **mat )
{
   for(int j=0; j<dim; j++)
      mat[row][j] = mat[row][j] + nbr*mat[col][j];
}

int maximum ( int dim, int **mat )
{
   int result = mat[0][0];

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) if(mat[i][j] > result) result = mat[i][j];

   return result;
}

void make_unimodular_matrix
 ( int dim, int size, bool poscff, int itr, int **mat, bool verbose )
{
   int row,col,nbr;

   for(int i=0; i<itr; i++)
   {
      random_row_and_column(dim,&row,&col);
      nbr = nonzero_integer(size,poscff);
      if(verbose)
      {
         cout << "   row index : " << row << endl;
         cout << "column index : " << col << endl;
         cout << "random value : " << nbr << endl;
      }
      unimodular_multiplication(dim,row,col,nbr,mat);
      if(verbose) write_exponent_matrix(dim,mat);
   }
}

void write_exponent_matrix ( int dim, int **mat )
{
   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++) cout << " " << mat[i][j];
      cout << endl;
   }
}

void read_exponent_matrix ( int dim, int **mat )
{
   cout << "reading a " << dim << "-by-" << dim << " matrix ..." << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) cin >> mat[i][j];
}

void lower_triangular_unit ( int dim, int **mat )
{
   for(int i=0 ; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         if(i < j)
            mat[i][j] = 0;  // zeros above the diagonal
         else
            mat[i][j] = 1;  // ones on and below the diagonal
      }
}

void two_variable_monomials ( int dim, int **mat )
{
   for(int i=0 ; i<dim; i++)
      for(int j=0; j<dim; j++) mat[i][j] = 0;

   int offset = dim % 2;  // for odd dimension insert one

   if(offset == 1) mat[0][0] = 1;

   for(int i=0; i<dim/2; i++)
   {
      int idx = offset + 2*i;
      mat[idx][idx] = 1;   mat[idx][idx+1] = 2;
      mat[idx+1][idx] = 2; mat[idx+1][idx+1] = 5;
   }
}

void make_monomial_system
 ( int dim, int size, int posvals, int nbritr,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int **rowsA, int vrblvl )
{
   for(int i=0; i<dim; i++)     // define the identity matrix
   {
      rowsA[i] = new int[dim];
      for(int j=0; j<dim; j++) rowsA[i][j] = 0;
      rowsA[i][i] = 1;
   }
   bool verbose = (vrblvl > 0);
   bool poscff = (posvals > 0);

   if(nbritr == -2)
      two_variable_monomials(dim,rowsA);
   else if(nbritr == -1)
      lower_triangular_unit(dim,rowsA);
   else if(nbritr > 0)
      make_unimodular_matrix(dim,size,poscff,nbritr,rowsA,0); // verbose);
   else
      read_exponent_matrix(dim,rowsA);

   if(verbose)
   {
      cout << "The unimodular matrix :" << endl;
      write_exponent_matrix(dim,rowsA);
      cout << "The largest exponent : " << maximum(dim,rowsA) << endl;
   }
   for(int i=0; i<dim; i++)
   {
      nvr[i] = 0;
      for(int j=0; j<dim; j++) if(rowsA[i][j] != 0) nvr[i]++;
   }
   for(int i=0; i<dim; i++)
   {
      idx[i] = new int[nvr[i]];
      exp[i] = new int[nvr[i]];
      int cnt = 0;
      for(int j=0; j<dim; j++)
         if(rowsA[i][j] != 0)
         {
            idx[i][cnt] = j;
            exp[i][cnt] = rowsA[i][j];
            cnt = cnt + 1;
         }
   }
   for(int i=0; i<dim; i++)
   {
      expfac[i] = new int[nvr[i]];
      common_factors(nvr[i],exp[i],&nbrfac[i],expfac[i]);
   }
   if(verbose)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "monomial " << i << ", nvr : " << nvr[i] << ", idx :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
         cout << ", exp :";
         for(int j=0; j<nvr[i]; j++) cout << " " << exp[i][j];
         if(nbrfac[i] > 0)
         {
            cout << ", expfac :";
            for(int j=0; j<nvr[i]; j++) cout << " " << expfac[i][j];
         } 
         cout << endl;
      }
   }
}

void prompt_dimensions
 ( int *dim, int *deg, int *size, int *posvals, int *vrblvl, int *nbritr,
   int *nbsteps )
{
   cout << "-> give the dimension : ";
   cin >> *dim;
   cout << "-> give the degree of the series : ";
   cin >> *deg;
   cout << "-> verbose level (0 for silent) : ";
   cin >> *vrblvl;
   cout << "MENU for the test problem :" << endl;
   cout << " -5 : all polynomials of the cyclic n-roots system" << endl;
   cout << " -4 : 2-column lower/upper triangular matrix" << endl;
   cout << " -3 : columns of the cyclic n-roots system" << endl;
   cout << " -2 : decoupled two variable monomials" << endl;
   cout << " -1 : lower triangular matrix of ones" << endl;
   cout << "  0 : user input" << endl;
   cout << "  n : number of unimodular multiplications" << endl;
   cout << "-> give the number for the test problem : ";
   cin >> *nbritr;
   if(*nbritr > 0)
   {
      cout << "-> give the size of the numbers : ";
      cin >> *size;
      // cout << "-> positive exponents ? (1 for yes) : ";
      // cin >> *posvals;
      *posvals = 1; // Laurent is not yet supported
   }
   else
   {
      *size = 2;
      *posvals = 1;
   }
   cout << "-> give the number of Newton steps : ";
   cin >> *nbsteps;
}

void copy_integer_matrix ( int dim, int **matfrom, int **matinto )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) matinto[i][j] = matfrom[i][j];
}

void matrix_matrix_multiply ( int dim, int **A, int **B, int **C )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         C[i][j] = 0;
         for(int k=0; k<dim; k++) C[i][j] = C[i][j] + A[i][k]*B[k][j];
      }
}

int posgcd ( int a, int b, int *ka, int *lb, int vrblvl )
{

   if(vrblvl>0) cout << "posgcd(" << a << "," << b << ")";

   int r = a % b;

   if(r == 0)
   {
      *ka = 0; *lb = 1;

      if(vrblvl > 0)
         cout << " = " << b << ", k = " << *ka << ", l = " << *lb << endl;

      return b;
   }
   else
   {
      int aa,bb,k1,l1,k0,l0,q,tmp;

      k0 = 0; l0 = 1;
      bb = r; aa = b;
      k1 = 1; l1 = (-a)/b;

      do
      {
         r = aa % bb;
         if(r == 0) break;
         q = aa / bb;
         tmp = k1;
         k1 = k0 - q*k1;
         k0 = tmp;
         tmp = l1;
         l1 = l0 - q*l1;
         l0 = tmp;
         aa = bb; bb = r;
      }
      while(true);
      *ka = k1; *lb = l1;

      if(vrblvl > 0)
         cout << " = " << b << ", k = " << *ka << ", l = " << *lb << endl;

      return bb;
   }
}

int gcd ( int a, int b, int *ka, int *lb, int vrblvl )
{
   if(vrblvl>0) cout << "gcd(" << a << "," << b << ")";

   if(a == 0)
   {
      if(b < 0)
      {
         *ka = 0; *lb = -1;

         if(vrblvl > 0)
            cout << " = " << -b << ", k = " << *ka << ", l = " << *lb << endl;

         return (-b);
      }
      else
      {
         *ka = 0; *lb = 1;

         if(vrblvl > 0)
            cout << " = " << b << ", k = " << *ka << ", l = " << *lb << endl;
         return b;
      }
   }
   else
   {
      if(b == 0)
      {
         if(a < 0)
         {
            *ka = -1; *lb = 0;

            if(vrblvl > 0)
               cout << " = " << -a << ", k = " << *ka << ", l = " << *lb << endl;

            return (-a);
         }
         else
         {
            *ka = 1; *lb = 0;

            if(vrblvl > 0)
               cout << " = " << a << ", k = " << *ka << ", l = " << *lb << endl;

            return a;
         }
      }
      else
      {
         if(a < b)
         {
            if(b < 0) 
            {
               int d = posgcd(-a,-b,ka,lb,vrblvl);
               *ka = -(*ka);
               *lb = -(*lb);
               return d;
            }
            if(a < 0)
            {
               int d;
               if(b > -a)
                  d = posgcd(b,-a,lb,ka,vrblvl);
               else
                  d = posgcd(-a,b,ka,lb,vrblvl);
               *ka = -(*ka);
               return d;
            }
            return posgcd(b,a,lb,ka,vrblvl);
         }
         else // a >= b
         {
            if(a < 0)
            {
               int d = posgcd(-a,-b,ka,lb,vrblvl);
               *ka = -(*ka);
               *lb = -(*lb);
               return d;
            }
            if(b < 0)
            {
               int d;
               if(a > -b)
                  d = posgcd(a,-b,ka,lb,vrblvl);
               else
                  d = posgcd(-b,a,lb,ka,vrblvl);
               *lb = -(*lb);
            }
            return posgcd(a,b,ka,lb,vrblvl);
         }
      }
   }
}

int lower_triangulate ( int dim, int **mat, int **uni, int vrblvl )
{
   if(vrblvl > 0)
   {
      cout << "the matrix on input :" << endl;
      write_exponent_matrix(dim, mat);
   }
   for(int i=0; i<dim; i++)                // initialize uni
   {
      for(int j=0; j<dim; j++) uni[i][j] = 0;
      uni[i][i] = 1;
   }
   int pivot;
   int col = 0;
   for(int i=0; i<dim; i++)
   {
      pivot = col-1;
      for(int j=col; j<dim; j++)           // look for pivot
      {
         if(vrblvl > 0)
            cout << "mat[" << i << "][" << j << "] : " << mat[i][j] << endl;
         if(mat[i][j] != 0) pivot = j;
         if(pivot >= col) break;
      }
      if(vrblvl > 0) cout << "pivot : " << pivot << endl;
      if(pivot < col)
         return -1;                        // zero column
      else
      {
         if(pivot != col)                  // swap columns
         {
            int temp;
            for(int k=0; k<dim; k++)
            {
               temp = mat[k][col];
               mat[k][col] = mat[k][pivot];
               mat[k][pivot] = temp;
               temp = uni[k][col];
               uni[k][col] = uni[k][pivot];
               uni[k][pivot] = temp;
            }
         }
         if(vrblvl > 0)
            cout << "making zeros, pivot = "
                 << pivot << ", col = " << col << endl;
         for(int j=col+1; j<dim; j++)      // make zeros
         {
            if(mat[i][j] != 0)
            {
               int ka,lb,d,aa,bb;  
               d = gcd(mat[i][col],mat[i][j],&ka,&lb,vrblvl);
               aa = mat[i][col]/d; 
               bb = mat[i][j]/d; 
               if(aa == bb)
               {
                  if(ka == 0)
                  {
                     ka = lb; lb = 0;
                  }
               }
               if(aa == -bb)
               {
                  if(ka == 0)
                  {
                     ka = -lb; lb = 0;
                  }
               }
               for(int k=i; k<dim; k++)
               {
                  int matkcol = mat[k][col];
                  int matkj = mat[k][j];
                  mat[k][col] = matkcol*ka    + matkj*lb;
                  mat[k][j]   = matkcol*(-bb) + matkj*aa;
               }
               for(int k=0; k<dim; k++)
               {
                  int unikcol = uni[k][col];
                  int unikj = uni[k][j];
                  uni[k][col] = unikcol*ka    + unikj*lb;
                  uni[k][j]   = unikcol*(-bb) + unikj*aa;
               }
            }
         }
      }
      col = col + 1;
   }
   return 0;
}

int exponent_forward_substitution ( int dim, int **mat, int *expsol )
{
   for(int i=0; i<dim; i++) expsol[i] = 1;

   if(mat[0][0] == 0) return -1;

   for(int i=1; i<dim; i++) // adjust exponents for component i
   {
      if(mat[i][i] == 0) return -1;

      for(int j=0; j<i; j++)
         expsol[i] = expsol[i] - mat[i][j]*expsol[j];
   }
   return 0;
}

int exponent_unimodular_transformation ( int dim, int **uni, int *expsol )
{
   int *result = new int[dim];

   for(int i=0; i<dim; i++) result[i] = 0;

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         result[i] = result[i] + uni[i][j]*expsol[j];

   for(int i=0; i<dim; i++) expsol[i] = result[i];

   free(result);

   return 0;
}

int exponents_check
 ( int dim, int **rowsA, int *expsol, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The matrix :" << endl;
      write_exponent_matrix(dim, rowsA);
   }
   int **copyA = new int*[dim];  // copy of A
   int **unimd = new int*[dim];  // unimodular transformation

   for(int i=0; i<dim; i++)      // initialize the data
   {
      unimd[i] = new int[dim];
      copyA[i] = new int[dim];
   }
   copy_integer_matrix(dim, rowsA, copyA);

   int sing = lower_triangulate(dim, copyA, unimd, 0);
   if(sing != 0) return sing;

   exponent_forward_substitution(dim, copyA, expsol);

   if(vrblvl > 0)
   {
      cout << "exponents after forward substitution :" << endl;
      for(int i=0; i<dim; i++) cout << " " << expsol[i];
      cout << endl;
   }
   exponent_unimodular_transformation(dim, unimd, expsol);

   if(vrblvl > 0)
   {
      cout << "exponents after unimodular transformation :" << endl;
      for(int i=0; i<dim; i++) cout << " " << expsol[i];
      cout << endl;
   }
   for(int i=0; i<dim; i++)
   {
      free(unimd[i]);
      free(copyA[i]);
   }
   free(unimd); free(copyA);

   return sing;
}

void row_sums ( int dim, int **rowsA, int *sums )
{
   for(int i=0; i<dim; i++)
   {
      sums[i] = 0;

      for(int j=0; j<dim; j++) sums[i] = sums[i] + rowsA[i][j];
   }
}
