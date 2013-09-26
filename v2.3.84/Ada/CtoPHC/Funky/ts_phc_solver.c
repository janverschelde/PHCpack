#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern void adafinal();
extern void _ada_phc_solver ( int n, int m, int *mc,
                              int ns, int *s, int nc, double *c );

void test_cosup_poly_system ( void );
/* interactive test on coefficient-support representation of poly systems */

void read_data_system ( int n, int m, int monsum, int moncnt[m] );
/* reads coefficients and supports for a system in n variables and m
   equations, where the i-th polynomial has as many monomials as moncnt[i],
   and the total number of monomials equals monsum */

int main(void)
{
   printf("\nCalling the blackbox polynomial system solver of PHCpack\n");

   test_cosup_poly_system();
 
   return 0;
}

void test_cosup_poly_system ( void )
{
   int n,m;

   printf("\nCreating a system from coefficients and supports.\n\n");

   printf("Give the number of variables : "); scanf("%d", &n);
   printf("Give the number of equations : "); scanf("%d", &m);

   {
      int i,moncnt[m],monsum = 0;

      printf("Reading the number of monomials in every polynomial...\n");
      for(i=0; i<m; i++)
      {
	 printf("  Give #monomials in polynomial %d : ", i+1);
         scanf("%d",&moncnt[i]);
      }
      for(i=0; i<m; i++)
         monsum += moncnt[i];

      read_data_system(n,m,monsum,moncnt);
   }
}

void read_data_system ( int n, int m, int monsum, int moncnt[m] )
{
   int dimsup = n*monsum;
   int dimcff = 2*monsum;
   int sup[dimsup];
   double cff[dimcff];
   int i,j,k,indsup,indcff;

   indsup = 0;
   indcff = 0;
   for(i=0; i<m; i++)
   {
      printf("Reading the support and coefficients");
      printf(" of polynomial %d ...\n", i+1);
      for(j=0; j<moncnt[i]; j++)
      {
         printf("  give exponents of monomial %d : ", j+1);
         for(k=0; k<n; k++) scanf("%d", &sup[indsup++]);
         printf("  give two doubles : ");
         scanf("%lf", &cff[indcff++]);
         scanf("%lf", &cff[indcff++]); 
      }
   }

   adainit();
   _ada_phc_solver(n,m,moncnt,dimsup,sup,dimcff,cff);
   adafinal();

}
