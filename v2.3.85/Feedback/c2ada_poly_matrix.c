/* file "c2ada__poly_matrix.c" contains definitions of prototypes in c2ada_poly_matrix.h" */
#include "c2ada_poly_matrix.h"

#include<stdlib.h>

double* c2ada_poly_matrix( int n, int m, POLY c[n][m], int l, int a[l], int *sum )
{
  int  i, j, k, s;
  double *b;
 
  k = 0;
  *sum = 0;
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
    {
       a[k]=c[i][j].d;
       *sum=*sum+c[i][j].d+1; 
       k++;
    }
  *sum = *sum *2;
  b = (double*) calloc(*sum, sizeof(double));
  k = 0;
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
      for(s=0; s<=c[i][j].d; s++)
      {  
        b[k++] = c[i][j].p[s].re;
        b[k++] = c[i][j].p[s].im;
      } 

  return b;
}
