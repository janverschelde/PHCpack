#include <stdio.h>

extern void _ada_track_paths
              ( int n, int m,
                int *p_ms, int p_ns, int *p_s, int p_nc, double *p_c,
                int *q_ms, int q_ns, int *q_s, int q_nc, double *q_c,
                int nbmul, int *mul, int nbcfs, int *cfs );
extern void adainit();
extern void adafinal();

int path_tracker ( int n, int m,
                   int *p_ms, int p_ns, int *p_s, int p_nc, double *p_c,
                   int *q_ms, int q_ns, int *q_s, int q_nc, double *q_c,
                   int nbmul, int *mul, int nbcfs, int *cfs )

/* This C function takes as input a target system, and a start system
   with start solutions.  It calls the Ada path tracking routine.  */

{
   printf("Passing the data to the Ada routines ...\n");
   /* adainit(); */
   _ada_track_paths(n,m,p_ms,p_ns,p_s,p_nc,p_c,
                        q_ms,q_ns,q_s,q_nc,q_c,nbmul,mul,nbcfs,cfs);
   /* adafinal(); */

   printf("... done with call to Ada routine.\n");

   return 0;
}
