#include <stdio.h>

extern double *_ada_phc_sys_rw ( int rw, int size_p, double *p );

void *call_phc_rw ( int n, int rw, double *x )
{
   int i,buffer_size = n;
   double *p;

   if (rw == 0)
   {
      p = _ada_phc_sys_rw(0,buffer_size,x);
      for (i=0; i<n; i++)
        x[i] = p[i];

      fflush(stdout);
      printf("\nHere is the system :\n");
      _ada_phc_sys_rw(1,buffer_size,x);
   }
   else 
      _ada_phc_sys_rw(1,buffer_size,x);
}
