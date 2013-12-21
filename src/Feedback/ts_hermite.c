#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "poly_dcmplx.h"
#include "poly_hermite.h"
#include "poly_matrix.h"
#include "dc_matrix.h"

void Call_Hermite ( int n, int m );

int main(void)
{
   int rows, cols;
   int num, i;
  
   srand(time(NULL)); 

   printf("Give number of rows : "); scanf("%d", &rows); printf("%d\n", rows);
   printf("Give number of columns : "); scanf("%d", &cols); printf("%d\n", cols); 

  /* Testing Hermite algorithm.  */
   Call_Hermite (rows, cols);
   
   /*

   printf("Testing Hermite algorithm.\n");
   printf("How many tests do you want : ");
   scanf("%d", &num);
   printf("Testing %d random cases.\n\n", num);
   for ( i=1; i<=num; i++)
   {
       printf("****************************************\n");
       printf("Now testing the NO. %d case:\n", i); 
	   Call_Hermite (rows, cols);

   }
   */
   return 0;
}

void Call_Hermite ( int n, int m )
{
    POLY a[n][m], a1[n][m];
    POLY p[n][n], t[n][m];
    dcmplx deter, dcmplx_p[n][n], x;
    int k;

    printf("1 random matrix.\n");
    printf("2 input your own matrix.\n");
    printf("Please choose the test matrix:");
    scanf("%d", &k ); printf("%d\n\n", k );

    if(k==1) random_matrix ( n, m, a);
    if(k==2) read_matrix( n, m, a );
    printf("the original matrix generated :\n");
    print(n,m,a); 
    
    zero_matrix ( n, m, a1);
    I_matrix (  n, p);
    copy ( n, m, a, t);

    /* Eliminate_Col(n,m,a,p,0,0); */
    /* Eliminate_Row(n,m,a,q,0,0); */
    Hermite(n, m, a, p);

    printf("The hermite form of matrix a is :\n");
    print(n,m,a);

    /* now begin to test the result */
    Multiply ( n,  n,  m, p, t, a1 ); 
 
    
    printf("The calculated hermite form with p*a is:\n");
    print(n, m, a1);
    printf(" p is:\n");
    print(n, n, p);

    x=create1(1.1);
    evaluate_matrix(n, n, p, dcmplx_p, x);
    deter=determinant(n, dcmplx_p);
    printf("The determinant of the p is: ");
    writeln_dcmplx(deter);
    
}













