#include <stdio.h>
#include<stdlib.h>
#include<time.h>

#include "dcmplx.h"
#include "poly_dcmplx.h"
#include "poly_matrix.h"

void Call_Inverse( int n);

int main(void)
{
   int n;

   printf("Give the dimension of the matrix: "); scanf("%d", &n);
   printf("%d\n", n);
   srand(time(NULL));
   Call_Inverse(n);
     
   return 0;
}

void Call_Inverse( int n)
{
  POLY ds;
  POLY M1[n][n], t_M1[n][n], product[n][n];
  int k; 

  printf("Please choose the test matrix:  "); 
  printf("1 random matrix.\n");
  printf("2 input your own matrix.\n");
  scanf("%d", &k ); printf("%d\n", k);

  if(k==1) random_matrix ( n, n, M1);
  if(k==2) read_matrix( n, n, M1 );
  printf("the original matrix generated :\n");
  print( n, n, M1);

  copy(n, n, M1, t_M1);
  ds = Inverse_Poly ( n, M1 );
  printf(" The inverse matrix of the matrix is :\n" );
  print(n, n, M1 );
  printf(" The polynomial ds is :\n" );
  Print_Poly( ds.d, ds.p );

  printf("The product (without ds) of the original matrix");
  printf(" and the inverse matrix is:\n");
 
  Multiply(n, n, n, M1, t_M1, product);
  print(n, n, product);
  
  free_matrix(n, n, M1);
  free_matrix(n, n, t_M1);
  free_matrix(n, n, product);
 
}

