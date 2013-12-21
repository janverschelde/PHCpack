#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dcmplx.h"
#include "dc_roots.h"
#include "dc_interpolation.h"
#include "poly_gcd.h"

void manual_test ( int n );
void random_test ( int n );
int equal ( dcmplx a, dcmplx b );

int main(void)
{
   int n, m, i;

   printf("How many tests you want: ");
   scanf("%d", &m);
   printf("Give the degree : ");
   scanf("%d", &n);

   manual_test(n); 
   /* for(i=1; i<=m; i++)
   {
     printf("\n******* The No. %d test *********\n", i);
     random_test (n);
     }*/
   return(0);
}

void manual_test ( int n )
{
   int i, j, m[n+1], k;
   dcmplx x[n+1], f[n+1], ff[n+1], tmp;
   double tol=1.0e-8;

   printf("multiple points test:\n");  
   printf("x[i] are multiple points,please give the f value like this:\n");
   printf("x[i], function value at the x[i],\n");
   printf("x[i+1], the first order derivative function value at x[i],\n");
   printf("x[i+2], the second order derivative function value at x[i].\n");
 
   for(i=0; i<=n; i++ ) 
   {
     printf("input x[%d]:\n", i);
     read_dcmplx(&x[i]);
     printf("input the multiplicity of the point\n");
     scanf("%d", &(m[i]));
     printf("input f[%d]:\n", i);
     read_dcmplx(&f[i]);
   }

   group_points(n+1, tol, x, f, m);
   for(i=0; i<=n; i++)
     ff[i]=f[i];
   divided_difference ( n+1,  x, f, m );

   printf("the result polynomial is:\n");
   for(j=0; j<=n; j++)
   {
      printf("f[%d]=", j);
      writeln_dcmplx(f[j]);
   }
   for(i=0; i<=n; i++)
   {
     printf("the point x with multiplicity %d is: ", m[i]);
     writeln_dcmplx(x[i]);
     printf("the given function value at point x is:");
     writeln_dcmplx(ff[i]);
     printf("the result polynomial value at point x is:");
     tmp=horner(n+1, f, x[i]);
     writeln_dcmplx(tmp);
     printf("\n");
     if(m[i]>1)  /*skip the multiple points*/
     {
       k=m[i];
       while(k>1)
       {
         i++;
         k--;
       }  
     }
   }
}

   
void random_test ( int n )
{
   dcmplx p[n+1], x[n+1], f[n+1];
   int i, m[i+1];

   srand(time(NULL));
   
 /* printf(" The original polinomial is: \n");*/

   for(i=0; i<n; i++)
   {      
      p[i] = random_dcmplx1();
    /*  writeln_dcmplx(p[i]); */
   }
   p[n]=create1(1.0);
   /*writeln_dcmplx(p[n]); */
   
   for(i=0; i<=n; i++)
   {
	  x[i] = random_dcmplx1();
          m[i] = 1;
	  f[i] = horner(n+1, p, x[i]);
   }

   divided_difference ( n+1,  x, f, m );

 /*  printf(" The polynomial after interpolation is: \n");
   for(i=0; i<=n; i++)
   {
	   writeln_dcmplx(f[i]);
   }   
*/
   for(i=0; i<=n; i++)
   {
	   if( !equal(p[i], f[i]) )
       {   
		   printf("i=%d\n", i);
		   printf("The original polynomial item is:");
		   writeln_dcmplx(p[i]);
		   printf("The polynomial after interpolation is:");
		   writeln_dcmplx(f[i]);
		   printf(" BUG!!\n ");
		   exit(0);
	   }
   }
   printf(" OK!\n ");
}

int equal ( dcmplx a, dcmplx b )
{
   double tol=1.0e-8; 
   dcmplx c;
   c = sub_dcmplx(a, b);
   if(dabs(c.re)<tol && dabs(c.im)<tol)
	   return 1;
   else
   {
	   printf("The difference is:");
	   writeln_dcmplx(c);
	   return 0;
   }
}




   
