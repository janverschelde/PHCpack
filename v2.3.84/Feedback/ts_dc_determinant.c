/* calculates the determinant of a double complex matrix 
 * based on LU-decomposition */

#include <stdio.h>
#include "dc_matrix.h"
#include "dcmplx.h"

void ts_deter(int n);
/* test the computation of the determinant on a user-given n-by-n matrix */


int main(void)
{
  int n;

  printf("Please give the dimension of the matrix:\n");
  scanf("%d", &n);
  ts_deter(n);

  return 0;
}

void ts_deter(int n)
{
   dcmplx ans, a[n][n];
  
   printf("Reading a %d-by-%d matrix...\n", n, n);
   read_dcmatrix(n,n,a);
   printf("Writing a %d-by-%d matrix : \n", n, n);
   print_dcmatrix(n,n,a);

   ans = determinant( n, a );

   printf("The determinant of the given matrix is:\n"); 
   writeln_dcmplx(ans);
}
