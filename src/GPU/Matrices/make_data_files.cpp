/* Make data files in double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "dbl_test_utilities.h"
#include "dbl_data_files.h"

using namespace std;

int make_real_matrix ( int dim );
/*
 * Prompts for a file name and writes the generated
 * real matrix of dimension dim to the file. */

int make_complex_matrix ( int dim );
/*
 * Prompts for a file name and writes the generated
 * complex matrix of dimension dim to the file. */

int main ( void )
{
   srand(time(NULL));

   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << endl;
   cout << "MENU for the data type :" << endl;
   cout << "  1. real numbers" << endl;
   cout << "  2. complex numbers" << endl;
   cout << "Type 1 or 2 to select the type : ";
   int datatype; cin >> datatype;

   cout << endl;

   if(datatype == 1)
      return make_real_matrix(dim);
   else
      return make_complex_matrix(dim);

   return 0;
}

int make_real_matrix ( int dim )
{
   double **A = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      A[i] = new double[dim];

      for(int j=0; j<dim; j++) A[i][j] = 0.0;
   }
   cout << "Give the name of a file : ";
   string filename; cin >> filename;

   cout << "-> making a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   dbl_random_upper_factor(dim,A);

   cout << "-> writing to file " << filename << " ..." << endl;

   dbl_write_matrix(filename,dim,A);

   return 0;
}

int make_complex_matrix ( int dim )
{
   double **Are = new double*[dim];
   double **Aim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Are[i] = new double[dim];
      Aim[i] = new double[dim];

      for(int j=0; j<dim; j++)
      {
         Are[i][j] = 0.0;
         Aim[i][j] = 0.0;
      }
   }
   cout << "Give the name of a file : ";
   string filename; cin >> filename;

   cout << "-> making a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   cmplx_random_upper_factor(dim,Are,Aim);

   cout << "-> writing to file " << filename << " ..." << endl;

   cmplx_write_matrix(filename,dim,Are,Aim);

   return 0;
}
