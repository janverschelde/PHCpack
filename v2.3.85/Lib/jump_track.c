/* The file jump_track.c contains the definitions of the functions
 * in jump_track.h */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#define v 0 /* verbose flag */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

int read_target_system_without_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(150,a,b,c);
   return fail;
}

int read_named_target_without_solutions ( int n, char *s )
{
   int fail,i,b[n];
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(161,&n,b,c);

   return fail;
}

int read_start_system_without_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(151,a,b,c);
   return fail;
}

int read_named_start_without_solutions ( int n, char *s )
{
   int fail,i,b[n];
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(162,&n,b,c);

   return fail;
}

int read_named_linear_product_start_system ( int n, char *s )
{
   int fail,i,b[n];
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(163,&n,b,c);

   return fail;
}

int silent_path_tracker
      ( int n, int *m, double *c,
        int *nbstep, int *nbfail, int *nbiter, int *nbsyst )
{
   int a[4],b[2],fail;
   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc(155,a,b,c);
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
   fail = _ada_use_c2phc(175,a,b,c);
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
   fail = _ada_use_c2phc(185,a,b,c);
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
   fail = _ada_use_c2phc(156,a,b,c);
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
   fail = _ada_use_c2phc(176,a,b,c);
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
   fail = _ada_use_c2phc(186,a,b,c);
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
   fail = _ada_use_c2phc(157,a,b,sol);
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
   fail = _ada_use_c2phc(177,a,b,sol);
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
   fail = _ada_use_c2phc(187,a,b,sol);
   if(fail == 0) (*k)++;
   return fail;
}
