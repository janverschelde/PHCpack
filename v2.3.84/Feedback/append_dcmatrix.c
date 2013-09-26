/* file append_dcmatrix.c provides an implementation of routines in append_dcmatrix.h */
#include "dcmplx.h"


void h_append ( int n, int m1, int m2, dcmplx a[n][m1], dcmplx b[n][m2], dcmplx c[n][m1+m2])
{
  int i, j;

  for(i=0; i<n; i++)
  {
    for(j=0; j<m1; j++)
      c[i][j]=a[i][j];  
    for(j=m1; j<m1+m2; j++)
      c[i][j]=b[i][j-m1];
  }  

}

void v_append ( int n1, int n2, int m, dcmplx a[n1][m], dcmplx b[n2][m], dcmplx c[n1+n2][m])
{
  int i, j;
 
  for(i=0; i<n1; i++)
    for(j=0; j<m; j++)
      c[i][j]=a[i][j];
  for(i=n1; i<n1+n2; i++)
    for(j=0; j<m; j++)
      c[i][j]=b[i-n1][j];
}  
