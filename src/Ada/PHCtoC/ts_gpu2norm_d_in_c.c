/* Calling the GPU from a main program in C... */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main ( int argc, char* argv[] )
{
   int dim,i;
   double *v;
   double result;

   printf("Give the dimension : "); scanf("%d",&dim);

   v = (double*)calloc(2*dim,sizeof(double));

   printf("The components of a random vector :\n");
   for(i=0; i<dim; i++)
   {
      double r = ((double) rand())/RAND_MAX;
      double angle = 2.0*M_PI*r;
      v[2*i] = cos(angle);
      v[2*i+1] = sin(angle);
      printf("%22.15le  %22.15le\n",v[2*i],v[2*i+1]);
   }

   int rtn = gpu2norm_d(dim,v,&result);

   printf("The result returned by the call : %.15le\n", result); 

   return 0;
}
