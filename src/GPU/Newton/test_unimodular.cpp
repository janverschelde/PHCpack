/* Tests making an integer matrix lower triangular. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing unimodular triangulation ..." << endl;

   cout << "give the dimension of the matrix : ";
   int dim; cin >> dim;

   int **rowsA = new int*[dim];  // rows of a matrix
   int **unimd = new int*[dim];  // unimodular matrix
   int **copyA = new int*[dim];  // copy of A
   int **prodm = new int*[dim];  // product matrix

   for(int i=0; i<dim; i++)      // initialize the data
   {
      rowsA[i] = new int[dim];
      unimd[i] = new int[dim];
      copyA[i] = new int[dim];
      prodm[i] = new int[dim];
   }
   cout << "give bound on the number size (0 for user input) : ";
   int size; cin >> size;

   if(size <= 0)
      read_exponent_matrix(dim, rowsA);
   else
   {
      cout << "generating an integer matrix of dimension " << dim
           << " ..." << endl;
      for(int i=0; i<dim; i++)      // initialize the data
         for(int j=0; j<dim; j++) rowsA[i][j] = rand() % size;
   } 
   cout << "The matrix :" << endl;
   write_exponent_matrix(dim, rowsA);
   copy_integer_matrix(dim, rowsA, copyA);

   int sing = lower_triangulate(dim, rowsA, unimd, 1);

   if(sing < 0)
      cout << "The matrix is singular!" << endl;
   else
   {
      cout << "The lower triangular matrix :" << endl;
      write_exponent_matrix(dim, rowsA);
      cout << "The unimodular transformation :" << endl;
      write_exponent_matrix(dim, unimd);
      matrix_matrix_multiply(dim, copyA, unimd, prodm);
      cout << "The product of the original with unimodular matrix :" << endl;
      write_exponent_matrix(dim, prodm);

      int *expsol = new int[dim];
      exponent_forward_substitution(dim, rowsA, expsol);
      cout << "exponents after forward substitution :" << endl;
      for(int i=0; i<dim; i++) cout << " " << expsol[i];
      cout << endl;
      exponent_unimodular_transformation(dim, unimd, expsol);
      cout << "exponents after unimodular transformation :" << endl;
      for(int i=0; i<dim; i++) cout << " " << expsol[i];
      cout << endl;
   }

   return 0;
}
