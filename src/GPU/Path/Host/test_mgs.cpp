// Test on the application of the Modified Gram-Schmidt method to solve
// a random linear system in double, double double, and quad double precision
// in the least squares sense.

#include <iostream>
#include <string>
#include <cmath>
#include "complexH.h"
#include "mgs_host.h"

using namespace std;

int double_test ( int rows, int cols );
// test in double precision for a rows-by-cols matrix

int double_double_test ( int rows, int cols );
// test in double double precision for a rows-by-cols matrix

int quad_double_test ( int rows, int cols );
// test in quad double precision for a rows-by-cols matrix

int main ( void )
{
   cout << "Testing the Modified Gram-Schmidt method ..." << endl;

   int rows, cols;

   cout << endl;
   cout << "Give the number of rows : "; cin >> rows;
   cout << "Give the number of columns : "; cin >> cols;
   cout << endl;

   char choice;

   cout << "Choose the precision :" << endl;
   cout << "  0. double precision" << endl;
   cout << "  1. double double precision" << endl;
   cout << "  2. quad double precision" << endl;
   cout << "Type 0, 1, or 2 : "; cin >> choice;
   cout << endl;

   if(choice == '0')
      double_test(rows,cols);
   else if(choice == '1')
      double_double_test(rows,cols);
   else if(choice == '2')
      quad_double_test(rows,cols);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return 0;
}

complexH<double> d_random_number ( void )
// returns a random complex number in double precision
{
   double rnd = (rand()/(double (RAND_MAX)))*2*M_PI;

   return complexH<double>(cos(rnd),sin(rnd));
}

complexH<dd_real> dd_random_number ( void )
// returns a random complex number in double double precision
{
   double rnd = (rand()/(double (RAND_MAX)))*2*M_PI;

   return complexH<dd_real>(cos(rnd),sin(rnd));
}

complexH<qd_real> qd_random_number ( void )
// returns a random complex number in quad double precision
{
   double rnd = (rand()/(double (RAND_MAX)))*2*M_PI;

   return complexH<qd_real>(cos(rnd),sin(rnd));
}

int double_test ( int rows, int cols )
{
   complexH<double> **mat; // coefficient matrix
   complexH<double> *rhs;  // right hand side vector
   complexH<double> *sol;  // solution, mat*sol = rhs
   // generate a random problem in traditional format
   mat = new complexH<double>*[rows];
   for(int i=0; i<rows; i++)
   {
      mat[i] = new complexH<double>[cols];
      for(int j=0; j<cols; j++) mat[i][j] = d_random_number();
   }
   sol = new complexH<double>[cols];
   for(int i=0; i<cols; i++) sol[i] = d_random_number();
   rhs = new complexH<double>[rows];
   for(int i=0; i<rows; i++)
   {
      rhs[i].init(0.0,0.0);
      for(int j=0; j<cols; j++) rhs[i] = rhs[i] + mat[i][j]*sol[j];
   }
   // convert the problem into column wise representation
   complexH<double>** vmt = new complexH<double>*[cols+1];
   for(int j=0; j<cols; j++)
   {
      vmt[j] = new complexH<double>[rows];
      for(int i=0; i<rows; i++) vmt[j][i] = mat[i][j];
   }
   vmt[cols] = new complexH<double>[rows];
   for(int i=0; i<rows; i++) vmt[cols][i] = rhs[i];
   complexH<double>** upp = new complexH<double>*[cols+1];
   for(int j=0; j<cols+1; j++) upp[j] = new complexH<double>[rows];
   complexH<double>* xls = new complexH<double>[cols];
   // apply least squares after Gram-Schmidt orthonormalization
   CPU_mgs2qrls<complexH<double>,double>(vmt,upp,xls,rows,cols+1);
   // check the errors
   for(int i=0; i<cols; i++)
   {
      cout << "sol[" << i << "] : " << sol[i];
      cout << "xls[" << i << "] : " << xls[i];
   }
   return 0;
}

int double_double_test ( int rows, int cols )
{
   complexH<dd_real> **mat; // coefficient matrix
   complexH<dd_real> *rhs; // right hand side vector
   complexH<dd_real> *sol; // solution, mat*sol = rhs
   // generate a random problem in traditional format
   mat = new complexH<dd_real>*[rows];
   for(int i=0; i<rows; i++)
   {
      mat[i] = new complexH<dd_real>[cols];
      for(int j=0; j<cols; j++) mat[i][j] = dd_random_number();
   }
   sol = new complexH<dd_real>[cols];
   for(int i=0; i<cols; i++) sol[i] = dd_random_number();
   rhs = new complexH<dd_real>[rows];
   for(int i=0; i<rows; i++)
   {
      rhs[i].init(0.0,0.0);
      for(int j=0; j<cols; j++) rhs[i] = rhs[i] + mat[i][j]*sol[j];
   }
   // convert the problem into column wise representation
   complexH<dd_real>** vmt = new complexH<dd_real>*[cols+1];
   for(int j=0; j<cols; j++)
   {
      vmt[j] = new complexH<dd_real>[rows];
      for(int i=0; i<rows; i++) vmt[j][i] = mat[i][j];
   }
   vmt[cols] = new complexH<dd_real>[rows];
   for(int i=0; i<rows; i++) vmt[cols][i] = rhs[i];
   complexH<dd_real>** upp = new complexH<dd_real>*[cols+1];
   for(int j=0; j<cols+1; j++) upp[j] = new complexH<dd_real>[rows];
   complexH<dd_real>* xls = new complexH<dd_real>[cols];
   // apply least squares after Gram-Schmidt orthonormalization
   CPU_mgs2qrls<complexH<dd_real>,dd_real>(vmt,upp,xls,rows,cols+1);
   // check the errors
   for(int i=0; i<cols; i++)
   {
      cout << "sol[" << i << "] : " << sol[i];
      cout << "xls[" << i << "] : " << xls[i];
   }
   return 0;
}

int quad_double_test ( int rows, int cols )
{
   complexH<qd_real> **mat; // coefficient matrix
   complexH<qd_real> *rhs;  // right hand side vector
   complexH<qd_real> *sol;  // solution, mat*sol = rhs
   // generate a random problem in traditional format
   mat = new complexH<qd_real>*[rows];
   for(int i=0; i<rows; i++)
   {
      mat[i] = new complexH<qd_real>[cols];
      for(int j=0; j<cols; j++) mat[i][j] = qd_random_number();
   }
   sol = new complexH<qd_real>[cols];
   for(int i=0; i<cols; i++) sol[i] = qd_random_number();
   rhs = new complexH<qd_real>[rows];
   for(int i=0; i<rows; i++)
   {
      rhs[i].init(0.0,0.0);
      for(int j=0; j<cols; j++) rhs[i] = rhs[i] + mat[i][j]*sol[j];
   }
   // convert the problem into column wise representation
   complexH<qd_real>** vmt = new complexH<qd_real>*[cols+1];
   for(int j=0; j<cols; j++)
   {
      vmt[j] = new complexH<qd_real>[rows];
      for(int i=0; i<rows; i++) vmt[j][i] = mat[i][j];
   }
   vmt[cols] = new complexH<qd_real>[rows];
   for(int i=0; i<rows; i++) vmt[cols][i] = rhs[i];
   complexH<qd_real>** upp = new complexH<qd_real>*[cols+1];
   for(int j=0; j<cols+1; j++) upp[j] = new complexH<qd_real>[rows];
   complexH<qd_real>* xls = new complexH<qd_real>[cols];
   // apply least squares after Gram-Schmidt orthonormalization
   CPU_mgs2qrls<complexH<qd_real>,qd_real>(vmt,upp,xls,rows,cols+1);
   // check the errors
   for(int i=0; i<cols; i++)
   {
      cout << "sol[" << i << "] : " << sol[i];
      cout << "xls[" << i << "] : " << xls[i];
   }
   return 0;
}
