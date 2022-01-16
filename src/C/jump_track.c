/* The file jump_track.c contains the definitions of the functions
 * in jump_track.h */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "jump_track.h"

#define v 0 /* verbose flag */

int read_target_system_without_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(150,a,b,c,0);
   return fail;
}

int read_dobldobl_target_system_without_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(872,a,b,c,0);
   return fail;
}

int read_quaddobl_target_system_without_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(873,a,b,c,0);
   return fail;
}

int read_named_target_without_solutions ( int n, char *s )
{
   int fail,i,b[n];
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(161,&n,b,c,0);

   return fail;
}

int read_start_system_without_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(151,a,b,c,0);
   return fail;
}

int read_named_start_without_solutions ( int n, char *s )
{
   int fail,i,b[n];
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(162,&n,b,c,0);

   return fail;
}

int read_named_linear_product_start_system ( int n, char *s )
{
   int fail,i,b[n];
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(163,&n,b,c,0);

   return fail;
}

int silent_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(155,a,b,c,0);
   *m = b[1];
   *nbstep = a[0]; *nbfail = a[1];
   *nbiter = a[2]; *nbsyst = a[3];
   return fail;
}

int silent_dobldobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(175,a,b,c,0);
   *m = b[1];
   *nbstep = a[0]; *nbfail = a[1];
   *nbiter = a[2]; *nbsyst = a[3];
   return fail;
}

int silent_quaddobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(185,a,b,c,0);
   *m = b[1];
   *nbstep = a[0]; *nbfail = a[1];
   *nbiter = a[2]; *nbsyst = a[3];
   return fail;
}

int reporting_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(156,a,b,c,0);
   *m = b[1];
   *nbstep = a[0]; *nbfail = a[1];
   *nbiter = a[2]; *nbsyst = a[3];
   return fail;
}

int reporting_dobldobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(176,a,b,c,0);
   *m = b[1];
   *nbstep = a[0]; *nbfail = a[1];
   *nbiter = a[2]; *nbsyst = a[3];
   return fail;
}

int reporting_quaddobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(186,a,b,c,0);
   *m = b[1];
   *nbstep = a[0]; *nbfail = a[1];
   *nbiter = a[2]; *nbsyst = a[3];
   return fail;
}

int write_next_solution_with_diagnostics
 ( int *k, int n, int m, double *sol,
   int nbstep, int nbfail, int nbiter, int nbsyst )
{
   int a[5],b[2],fail;

   b[0] = n; b[1] = m;
   a[0] = nbstep; a[1] = nbfail; a[2] = nbiter; a[3] = nbsyst;
   a[4] = *k;
   fail = _ada_use_c2phc4c(157,a,b,sol,0);
   if(fail == 0) (*k)++;
   return fail;
}

int write_next_dobldobl_solution_with_diagnostics
 ( int *k, int n, int m, double *sol,
   int nbstep, int nbfail, int nbiter, int nbsyst )
{
   int a[5],b[2],fail;

   b[0] = n; b[1] = m;
   a[0] = nbstep; a[1] = nbfail; a[2] = nbiter; a[3] = nbsyst;
   a[4] = *k;
   fail = _ada_use_c2phc4c(177,a,b,sol,0);
   if(fail == 0) (*k)++;
   return fail;
}

int write_next_quaddobl_solution_with_diagnostics
 ( int *k, int n, int m, double *sol,
   int nbstep, int nbfail, int nbiter, int nbsyst )
{
   int a[5],b[2],fail;

   b[0] = n; b[1] = m;
   a[0] = nbstep; a[1] = nbfail; a[2] = nbiter; a[3] = nbsyst;
   a[4] = *k;
   fail = _ada_use_c2phc4c(187,a,b,sol,0);
   if(fail == 0) (*k)++;
   return fail;
}

int standard_crude_tracker ( int verbose )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(622,&verbose,b,c,0);

   return fail;
}

int dobldobl_crude_tracker ( int verbose )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(623,&verbose,b,c,0);

   return fail;
}

int quaddobl_crude_tracker ( int verbose )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(624,&verbose,b,c,0);

   return fail;
}
