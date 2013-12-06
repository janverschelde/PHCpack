/* file "c2ada_dc_matrix.c" contains definitions of prototypes in c2ada_dc_matrix.h" */
#include "c2ada_dc_matrix.h"

void c2ada_dc_matrix( int n, int m, dcmplx c[n][m], int db, double b[db], int start)
{
  int  i, j, k;

  k = start;
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
    {
      b[k++] = c[i][j].re;
      b[k++] = c[i][j].im; 
    }
}
 

