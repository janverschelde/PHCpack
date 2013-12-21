#include <stdio.h>
#include<stdlib.h>
#include<time.h>

#include "dcmplx.h"
#include "poly_dcmplx.h"
#include "realization.h"
#include "poly_matrix.h"
#include "dc_matrix.h"
#include "dc_inverse.h"

void Call_Realization ( int m, int p, int q );

int main(void)
{
   int rows, cols, q;
  
   printf("Give the number of input m :\n"); scanf("%d", &rows);
   printf("m=%d\n", rows ); 
   printf("Give the number of output p : \n"); scanf("%d", &cols);
   printf("p=%d\n", cols );
   printf("Give the number of internal states :\n"); scanf("%d", &q);
   printf("q=%d\n", q );
   srand(time(NULL));
   Call_Realization( rows, cols, q );
   
   return 0;
}

void Call_Realization( int m , int p, int q )
{
   POLY M1[p][p], M2[m][p];
   dcmplx Ac[q][q], Bc[q][p], Cc[m][q], Dc[m][p];
   int  k;
   
   int i, j;
   dcmplx Iq[q][q], Iq_F[q][q], t_tmp[m][q], t_tmp1[m][p], t_tmp2[m][p];
   dcmplx s =create2(-0.23423423423, 0);

   printf("Please choose the test matrix:\n");
   printf("1 random matrix.\n");
   printf("2 input your own matrix.\n");
   scanf("%d", &k );
   if(k==1) { 
      random_matrix ( p, p, M1);
      random_matrix ( m, p, M2);
   }

   if(k==2) {
      printf("input M1:\n");
      read_matrix( p, p, M1 );
      printf("input M2:\n");
      read_matrix( m, p, M2 );
   }
   printf("\nthe original matrix M1:\n");
   print( p, p, M1); 
   printf("the original matrix M2:\n");
   print( m, p, M2);

   /******************************/
   realization( p, m, q, M1, M2, Ac, Bc, Cc, Dc); 
   printf("The final realization result:\n"); 
 
   printf("\nThe Ac is:\n");
   print_dcmatrix(q, q, Ac);
  
   printf("\nThe Bc is:\n");
   print_dcmatrix(q, p, Bc);

   printf("\nThe Cc is:\n");
   print_dcmatrix(m, q, Cc);

   printf("\nThe Dc is:\n");
   print_dcmatrix(m, p, Dc);

   /* check the results */ 
   for(i=0; i<q; i++)
     for(j=0; j<q; j++)
      {
        if(i==j) Iq[i][j] = s;
        else Iq[i][j] = zero;
      }
   sub_dcmatrix(q, q, Iq, Ac, Iq_F);
   dcinverse(q, Iq_F);
   multiply_dcmatrix(m, q, q, Cc, Iq_F, t_tmp);

   multiply_dcmatrix(m, q, p, t_tmp, Bc, t_tmp1);
   add_dcmatrix(m, p, t_tmp1, Dc, t_tmp2);
    
   printf("\nThe calculated H matrix at point s is:\n");
   print_dcmatrix(m, p, t_tmp2);  
}  



