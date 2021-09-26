/* Make data files in double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "dbl_test_utilities.h"
#include "dbl_data_files.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "Give the dimension : ";
   int dim; cin >> dim;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      A[i] = new double[dim];
      for(int j=0; j<dim; j++) A[i][j] = 0.0;
         
   }
   cout << "Give the name of a file : ";
   string filename; cin >> filename;

   cout << "Making a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   dbl_random_upper_factor(dim,A);

   cout << "Writing to file " << filename << " ..." << endl;

   dbl_write_matrix(filename,dim,A);

   return 0;
}
