// Test on the application of the Modified Gram-Schmidt method to solve
// a random linear system in double, double double, and quad double precision
// in the least squares sense.

#include <iostream>
#include <string>
#include <cmath>
#include "complexH.h"
#include "mgs_host.h"

using namespace std;

template <class ComplexType, class RealType>
ComplexType random_complex ();
// Returns a random complex number.

template <class ComplexType, class RealType>
int test ( int rows, int cols );
// test in double precision for a random rows-by-cols matrix

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
      //double_test(rows,cols);
      test<complexH<double>,double>(rows,cols);
   else if(choice == '1')
      //double_double_test(rows,cols);
      test<complexH<dd_real>,dd_real>(rows,cols);
   else if(choice == '2')
      //quad_double_test(rows,cols);
      test<complexH<qd_real>,qd_real>(rows,cols);
   else
      cout << "Invalid choice " << choice << " for the precision." << endl; 

   return 0;
}

template <class ComplexType, class RealType>
ComplexType random_complex ()
{
   double angle = 2*M_PI*((double) rand())/RAND_MAX;
   double re_part = cos(angle);
   double im_part = sin(angle);

   RealType real_part = RealType(re_part);
   RealType imag_part = RealType(im_part);

   ComplexType result = ComplexType(real_part, imag_part);

   return result;
}

template <class ComplexType, class RealType>
int test ( int rows, int cols )
{
   ComplexType** mat; // coefficient matrix
   ComplexType* rhs;  // right hand side vector
   ComplexType* sol;  // solution, mat*sol = rhs
   // generate a random problem in traditional format
   mat = new ComplexType*[rows];
   for(int i=0; i<rows; i++)
   {
      mat[i] = new ComplexType[cols];
      for(int j=0; j<cols; j++)
         mat[i][j] = random_complex<ComplexType,RealType>();
   }
   sol = new ComplexType[cols];
   for(int i=0; i<cols; i++)
      sol[i] = random_complex<ComplexType,RealType>();
   rhs = new ComplexType[rows];
   for(int i=0; i<rows; i++)
   {
      rhs[i].init(0.0,0.0);
      for(int j=0; j<cols; j++) rhs[i] = rhs[i] + mat[i][j]*sol[j];
   }
   // convert the problem into column wise representation
   ComplexType** vmt = new ComplexType*[cols+1];
   for(int j=0; j<cols; j++)
   {
      vmt[j] = new ComplexType[rows];
      for(int i=0; i<rows; i++) vmt[j][i] = mat[i][j];
   }
   vmt[cols] = new ComplexType[rows];
   for(int i=0; i<rows; i++) vmt[cols][i] = rhs[i];
   ComplexType** upp = new ComplexType*[cols+1];
   for(int j=0; j<cols+1; j++) upp[j] = new ComplexType[rows];
   ComplexType* xls = new ComplexType[cols];
   // apply least squares after Gram-Schmidt orthonormalization
   CPU_mgs2qrls<ComplexType,RealType>(vmt,upp,xls,rows,cols+1);
   // check the errors
   for(int i=0; i<cols; i++)
   {
      cout << "generated sol[" << i << "] : " << sol[i];
      cout << " computed xls[" << i << "] : " << xls[i];
   }
   return 0;
}
