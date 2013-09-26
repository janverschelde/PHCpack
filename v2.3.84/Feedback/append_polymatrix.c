/* file append_polymatrix.c provides an implementation of routines in append_polymatrix.h */
#include "append_polymatrix.h"

void poly_h_append ( int n, int m1, int m2, POLY a[n][m1], POLY b[n][m2], POLY c[n][m1+m2])
{
  int i, j;

  for(i=0; i<n; i++)
  {
    for(j=0; j<m1; j++)
      c[i][j]=assign_poly(a[i][j]);  
    for(j=m1; j<m1+m2; j++)
      c[i][j]=assign_poly(b[i][j-m1]);
  }  

}

void poly_v_append ( int n1, int n2, int m, POLY a[n1][m], POLY b[n2][m], POLY c[n1+n2][m])
{
  int i, j;
 
  for(i=0; i<n1; i++)
    for(j=0; j<m; j++)
      c[i][j]=assign_poly(a[i][j]);
  for(i=n1; i<n1+n2; i++)
    for(j=0; j<m; j++)
      c[i][j]=assign_poly(b[i-n1][j]);
}  
