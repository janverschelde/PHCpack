#include <stdio.h>

int pieri_sols ( int m, int p, int q, int nbsols,
                 int ns, int s[ns], int nc, double c[nc],
                 int npt, double pt[npt], int npl, double pl[npl],
                 char *filename )

/* This C function takes as input degree q solution maps producing
   p-planes that meet m-planes at specified interpolation points.  */

{
   int i;

   printf("Dimension of input planes : %d.\n", m);
   printf("Dimension of output planes : %d.\n", p);
   printf("Degree of solution maps : %d.\n", q);
   printf("Number of solution maps : %d.\n", nbsols);
   printf("Length of the degrees : %d.\n", ns); 
   printf("The degrees : ");
   for(i=0; i<ns; i++)
      printf(" %d", s[i]);  
   printf(".\n");
   printf("Number of coefficients : %d.\n", nc);
   printf("The coefficients :\n");
   for(i=0; i<nc; i++)
      printf(" %.15le\n", c[i]);

   printf("Length of interpolation points : %d.\n", npt);
   printf("The interpolation points :\n");
   for(i=0; i<npt; i++)
      printf(" %.15le\n", pt[i]);
   printf("Number of coefficient in input planes : %d.\n", npl);
   printf("The coefficients of the input planes :\n");
   for(i=0; i<npl; i++)
      printf(" %.15le\n", pl[i]);

   printf("The name of the file : %s.\n", filename);

   return 0;
}
