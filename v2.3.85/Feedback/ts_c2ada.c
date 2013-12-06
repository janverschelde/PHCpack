#include <stdio.h>
#include "c2ada_poly_matrix.h"
#include "poly_matrix.h"

void call_c2ada(int n, int m);

int main ()
{
 int n, m ;

 /* printf("please give the rows'number of the polynomial matrix.\n"); */
 scanf("%d", &n);
 /* printf("please give the columns'number of the polynomial matrix.\n"); */
 scanf("%d", &m); 
 call_c2ada(n, m);

}

void call_c2ada(int n, int m)
{
   POLY c[n][m];
   int a[n*m], l;
   int sum, i;
   double * b;

   l = n*m;
   read_matrix1(n, m, c);

   b=c2ada_poly_matrix( n, m, c, l, a, &sum);
   printf("%d %d ", n, m );
   for(i=0; i<l; i++)
      printf("%d ", a[i]);
   printf("\n");
   for(i=0; i<sum; i++)
   {
      printf("%.15le ", b[i]);
      if((i+1)%4==0)
        printf("\n");
   }
   free_matrix(n, m, c );
}
