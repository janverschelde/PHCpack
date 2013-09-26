/* Implementation of routines in dc_lnverse.h */

#include "dcmplx.h"
#include "dc_matrix.h"
#include "dc_inverse.h"

void dcinverse ( int n, dcmplx a[n][n] )
{
   dcmplx I[n], lu[n][n];
   int ipvt[n], info, i, j;

   copy_dcmatrix(n, n, a, lu);
   info = lufac(n, lu , ipvt);

   if (info == -1)
   {
      for( i=0; i<n; i++ )
      {
         for( j=0; j<n; j++ )
            I[j] = zero;
         I[i] = one;
         lusolve(n, lu, ipvt, I);  
         for( j=0; j<n; j++ )
            a[j][i] = I[j];    
      }
   }

   else
      printf("info = %d, the matrix is singular.\n", info);
}  
  


