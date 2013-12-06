#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "poly_dcmplx.h"
#include "poly_smith.h"
#include "poly_matrix.h"
#include "dc_matrix.h"

void Call_Smith ( int n, int m );

int main(void)
{
   int rows, cols;
   int num, i;
  
   srand(time(NULL)); 

   printf("Give number of rows : "); scanf("%d", &rows);
   printf("Give number of columns : "); scanf("%d", &cols);  

  /* Testing Smith algorithm.  */
   Call_Smith (rows, cols);
   
   /*

   printf("Testing Smith algorithm.\n");
   printf("How many tests do you want : ");
   scanf("%d", &num);
   printf("Testing %d random cases.\n\n", num);
   for ( i=1; i<=num; i++)
   {
       printf("****************************************\n");
       printf("Now testing the NO. %d case:\n", i); 
	   Call_Smith (rows, cols);

   }
   */
   return 0;
}

void Call_Smith ( int n, int m )
{
    POLY a[n][m], a1[n][m], a2[n][m];
    POLY p[n][n];
    POLY q[m][m];
    POLY t[n][m];
    int k;
    dcmplx deter_p, deter_q, pp[n][n], qq[m][m];

    printf("Please choose the test matrix:\n");
    printf("1 random matrix.\n");
    printf("2 input your own matrix.\n");
    scanf("%d", &k );
    if(k==1) random_matrix ( n, m, a);
    if(k==2) read_matrix( n, m, a );
    printf("the original matrix generated :\n");
    print(n,m,a); 
    
    zero_matrix ( n, m, a1);
    zero_matrix ( n, m, a2);
    I_matrix (  n, p);
    I_matrix ( m, q);
    copy ( n, m, a, t);

    Smith(n, m, a, p, q);

    printf("The smith form of matrix a is :\n");
    print(n,m,a);

    /* now begin to test the result */
    Multiply ( n,  n,  m, p, t, a1 ); 
    Multiply ( n,  m,  m, a1, q, a2 );
    
    printf("The calculated smith form with p*a*q is:\n");
    print(n, m, a2);

    printf(" p is:\n");
    print(n, n, p);
    evaluate_matrix( n, n, p, pp, one);
    deter_p=determinant(n, pp);
    printf("the determinant of p is: \n");
    writeln_dcmplx(deter_p);

    printf(" q is:\n");
    print(m, m, q);
    evaluate_matrix( m, m, q, qq, one);
    deter_q=determinant(m, qq);
    printf("the determinant of q is \n");
    writeln_dcmplx(deter_q);        
}













