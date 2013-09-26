#include<stdio.h>
#include<stdlib.h>

extern void _ada_cosupoly ( int n, int *s, int m, double *c );
extern void adainit();
extern void adafinal();

void test_cosup_poly ( void );
/* interactive test on coefficient-support representation of polynomials */

void read_data ( int n, int m );
/* reads coefficients and support for a polynomial in n variables and
   with complex coefficients with m terms */

int main(void)
{
   printf("\nTesting Coefficient Support Polynomial Representations\n");

   test_cosup_poly();
 
   return 0;
}

void test_cosup_poly ( void )
{ 
   int n,m;

   printf("\nCreating a polynomial from coefficients and support.\n\n");

   printf("Give the number of variables : "); scanf("%d", &n);
   printf("Give the number of monomials : "); scanf("%d", &m);

   read_data(n,m);
}

void read_data ( int n, int m )
{
   double *cff;
   int *sup;
   int i,j,indcff,indsup;

   cff = (double*)calloc(2*m,sizeof(double));
   sup = (int*)calloc(n*m,sizeof(int));

   indcff = 0;
   indsup = 0;
   for(i=0; i<m; i++)
   {
      printf("Reading coefficient and exponents for monomial %d ...\n", i+1);
      printf("  coefficient %d : ", i+1);
      scanf("%lf", &cff[indcff++]);
      scanf("%lf", &cff[indcff++]);
      printf("  exponent vector : ");
      for(j=0; j<n; j++)
	 scanf("%d", &sup[indsup++]); 
   }

   printf("The coefficients :\n");
   for(i=0; i<2*m; i++)
      printf("%.15lf\n", cff[i]);
   printf("The support :\n");
   for(i=0; i<n*m; i++)
      printf(" %d", sup[i]);
   printf("\n");

   printf("Calling Ada ...\n");
   adainit();
   _ada_cosupoly(n*m,sup,2*m,cff);
   adafinal();
   printf("... done with the call.\n");
}
