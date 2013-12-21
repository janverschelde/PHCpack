#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern double* _ada_phc_sol_rw ( int rw, int size_p, double *p );
extern void adafinal();

#define buffer_size 1000

int main(void)
{
   printf("\nTesting Reading of Solution Lists\n");

   adainit();
   {
      int i;
      double *p;
      double *s = (double*) calloc(buffer_size,sizeof(double));
      s = _ada_phc_sol_rw(0,buffer_size,s);
      fflush(stdout);
      printf("\nHere is the solution list :\n");
      p = _ada_phc_sol_rw(1,buffer_size,s);
   }
   adafinal();

   return 0;
}
