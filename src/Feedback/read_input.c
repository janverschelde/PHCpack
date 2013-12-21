#include <stdio.h>
#include "c2ada_dc_matrix.h"
#include "dc_matrix.h"
#include "dc_inverse.h"

void read_input(int n, int m, int p, int q, int nn);

int main ()
{
 int n, m, p, q, nn;
 /* printf("please give the number of the internal states for the given plant (A, B, C)\n"); */
 scanf("%d", &n);
 /* printf("please give the system's input dimension.\n"); */
 scanf("%d", &m);
 /* printf("please give the system's output dimension.\n"); */
 scanf("%d", &p);
 /* printf("please give the number of the internal states for the dynamic compensators.\n"); */
 scanf("%d", &q);
 nn=m*p+q*(m+p);
 read_input(n, m, p, q, nn);

}

void read_input(int n, int m, int p, int q, int nn)
{
  dcmplx A[n][n], B[n][m], C[p][n];
  dcmplx s[nn], Is[n][n], Is_A[n][n], tmp[p][n], M[p][m];  
  double a[nn*2], b[nn*p*m*2];
  int i, j, k, start;

  read_dcmatrix1(n, n, A);
  read_dcmatrix1(n, m, B);
  read_dcmatrix1(p, n, C);
  
  for(i=0; i<nn; i++)
    read_dcmplx(&s[i]);

  j = 0;
  for(i=0; i<nn; i++)
  {
    a[j++] = s[i].re;
    a[j++] = s[i].im;
  }

  start = 0;
  for(k=0; k<nn; k++)
  { 
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
      {
        if(i==j) Is[i][j] = s[k];
        else Is[i][j] = zero;
      }
   sub_dcmatrix(n, n, Is, A, Is_A);
   dcinverse(n, Is_A);
   multiply_dcmatrix(p, n, n, C, Is_A, tmp);
   multiply_dcmatrix(p, n, m, tmp, B, M);
   c2ada_dc_matrix( p, m, M, nn*p*m*2, b, start);
   start = start + p*m*2;  
 
   printf("The M%d matrix is:\n", k);
   print_dcmatrix(p, m, M);
  }

  printf("the double array for s is:\n");
  for(i=0; i<nn*2; i++)
  {
     printf("%.15le ", a[i]);
     if((i+1)%4==0)
        printf("\n");
  }

  printf("the double array for M is:\n");
   for(i=0; i<nn*p*m*2; i++)
   {
      printf("%.15le ", b[i]);
      if((i+1)%4==0)
        printf("\n");
   }


}  

