#include<stdio.h>
#include<stdlib.h>

/* Important notice :
  Since this routine will be called after already one
  C function has called an Ada routine, the adainit and adafinal
  must not be invoked here, otherwise the program raises an error
  during run time. */

/* extern void adainit(); */
/* extern void adafinal(); */

extern void _ada_getsys3 ( int n, int m, int moncnt[n],
                           int ns, int s[ns], int nc, double c[nc] );

int getsys2 ( int n, int m, int moncnt[n],
              int ns, int s[ns], int nc, double c[nc] )

/* This C function takes as input the support and coefficients of
   a system of n polynomial equations in m variables. */

{
   int i;

   printf("\nThis is your friendly putsys...\n");

   printf("The number of equations : %d.\n", n);
   printf("The number of variables : %d.\n", m);
   printf("#monomials :");
   for(i=0; i<n; i++)
      printf(" %d",moncnt[i]);
   printf("\n");
   printf("Length of the support : %d.\n", ns); 
   printf("The support : ");
   for(i=0; i<ns; i++) printf(" %d", s[i]);  
   printf(".\n");
   printf("Number of coefficients : %d.\n", nc);
   printf("The coefficients :\n");
   for(i=0; i<nc; i++) printf(" %.15lf\n", c[i]);

/*   adainit(); */
   _ada_getsys3(n,m,moncnt,ns,s,nc,c);
/*   adafinal(); */

   printf("... leaving the C function.\n");

   return 0;
}
