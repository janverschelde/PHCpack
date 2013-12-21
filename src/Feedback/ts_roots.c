#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dc_roots.h"

void interactive_test ( int n );

/* DESCRIPTION :
      Reads the coefficients of a polynomial of degree n 
      and computes the roots. */

void random_test ( int n, int m );

/* DESCRIPTION :
      Generates a random complex polynomial of degree n
      and computes its roots; does this test m times. */

int test_roots ( int n, dcmplx p[n+1] );

/* DESCRIPTION :
      Computes the roots of the polynomial p,
      returns 1 if the accuracy is not reached,
      otherwise 0 is returned. */

int main(void)
{
   int n,m;

   printf("Give the degree : ");
   scanf("%d",&n);

   printf("Give number of tests (0 for interaction) : ");
   scanf("%d",&m);

   if (m == 0)
      interactive_test(n);
   else
      random_test(n,m);

   return 0;
}

void interactive_test ( int n )
{
   dcmplx p[n+1];
   int i,k,bug = 0;

   printf("Reading the coefficients...\n");
   for(i=0; i<=n; i++)
   {
      printf("  give complex coefficient for x^%d : ", i);
      read_dcmplx(&p[i]);
   }
   bug = test_roots(n,p);

   if (bug == 0)
      printf("Computed roots without bug.\n");
   else
      printf("Bug found!?  Accuracy requirement not reached.\n");
} 

int test_roots ( int n, dcmplx p[n+1] )
{
   dcmplx r[n],res;
   double eps = 1.0e-12;
   double tol = 1.0e-8;
   int i,bug = 0;
   int m[n];
   double mul_tol = 1.0e-2;

   multiple_roots(n+1,p,eps,10*n,r,mul_tol,m);

   printf("The roots of p with multiplicities :\n");
   for(i=0; i<n; i++)
   {
      printf("r[%d] = ", i);
      write_dcmplx(r[i]);
      printf("  m = %d\n", m[i]);
   }
   printf("The function evaluated at the roots :\n");
   for(i=0; i<n; i++)
   {
      printf("p(r[%d]) = ", i);
      res = horner(n+1,p,r[i]);
      writeln_dcmplx(res);
      if (dcabs(res) > tol)
      {
         bug = 1;
         printf("BUG FOUND\n");
         break;
      }
      if (bug == 1) break;
   }
   return bug;
}

void random_test ( int n, int m )
{
   dcmplx p[n+1];
   int i,k,bug = 0;

   srand(time(NULL));

   for (k=0; (k<m) && (bug == 0); k++)
   {
      p[n] = create1(1.0);
      for(i=0; i<n; i++)
         p[i] = random_dcmplx1();

      bug = test_roots(n,p);
   }

   if (bug == 0)
      printf("Computed %d random cases without bug.\n", m);
   else
      printf("Bug found at case %d.\n", k);
}
