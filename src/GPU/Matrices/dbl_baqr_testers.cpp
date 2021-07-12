// The file dbl_baqr_testers.cpp defines the function with prototypes in
// the file dbl_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random_matrices.h"
#include "dbl_test_utilities.h"

using namespace std;

void test_real_blocked_qr ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   // random_dbl_upper_matrix(dim,dim,A);
   dbl_random_upper_factor(dim,A);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
}
