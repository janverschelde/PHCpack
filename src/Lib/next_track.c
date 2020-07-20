/* The file next_track.c defines the functions of next_track.h. */

#include <stdlib.h>
#include <stdio.h>
#include "next_track.h"

int initialize_standard_homotopy
 ( int fixed_gamma, double regamma, double imgamma )
{
   int fail,*b;
   double c[2];

   c[0] = regamma;
   c[1] = imgamma;

   fail = _ada_use_c2phc4c(500,&fixed_gamma,b,c,0);

   return fail;
}

int initialize_dobldobl_homotopy
 ( int fixed_gamma, double regamma, double imgamma )
{
   int fail,*b;
   double c[2];

   c[0] = regamma;
   c[1] = imgamma;

   fail = _ada_use_c2phc4c(501,&fixed_gamma,b,c,0);

   return fail;
}

int initialize_quaddobl_homotopy
 ( int fixed_gamma, double regamma, double imgamma )
{
   int fail,*b;
   double c[2];

   c[0] = regamma;
   c[1] = imgamma;

   fail = _ada_use_c2phc4c(502,&fixed_gamma,b,c,0);

   return fail;
}

int initialize_multprec_homotopy ( int fixed_gamma, int decimals )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(512,&fixed_gamma,&decimals,c,0);

   return fail;
}

int initialize_varbprec_homotopy 
 ( int fixed_gamma, int nc_target, char *target, int nc_start, char *start )
{
   int fail,i;
   const int len = nc_target + nc_start;
   int b[len];
   int a[3];
   double *c;

   a[0] = fixed_gamma;
   a[1] = len;
   a[2] = nc_target; 

   for(i=0; i<nc_target; i++) b[i] = (int) target[i];
   for(i=0; i<nc_start; i++) b[nc_target+i] = (int) start[i];
 
   fail = _ada_use_c2phc4c(516,a,b,c,0);

   return fail;
}

int initialize_standard_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(503,&k,b,c,0);

   return fail;
}

int initialize_dobldobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(504,&k,b,c,0);

   return fail;
}

int initialize_quaddobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(505,&k,b,c,0);

   return fail;
}

int initialize_multprec_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(513,&k,b,c,0);

   return fail;
}

int initialize_varbprec_solution ( int nv, int nc, char *sol )
{
   int fail,i;
   int a[2];
   int b[nc];
   double *c;

   a[0] = nc;
   a[1] = nv;
   for(i=0; i<nc; i++) b[i] = (int) sol[i];

   fail = _ada_use_c2phc4c(517,a,b,c,0);

   return fail;
}

int next_standard_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(506,&k,b,c,0);

   return fail;
}

int next_dobldobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(507,&k,b,c,0);

   return fail;
}

int next_quaddobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(508,&k,b,c,0);

   return fail;
}

int next_multprec_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc4c(514,&k,b,c,0);

   return fail;
}

char *next_varbprec_solution
 ( int want, int maxprc, int maxitr, int verbose, int *nc, int *fail )
{
   int a[4],len;
   double *c;
   char *sol;

   a[0] = want;
   a[1] = maxprc;
   a[2] = maxitr;
   a[3] = verbose;

   *fail = _ada_use_c2phc4c(518,a,&len,c,0);
   {
      int i,b[len];
      *fail = _ada_use_c2phc4c(520,a,b,c,0);
      if(len != a[0])
         *fail = 518;
      else
      {
         sol = (char*) calloc(len,sizeof(char));
         for(i=0; i<a[0]; i++) sol[i] = ' ';
         for(i=0; i<a[0]; i++) sol[i] = (char) b[i];
         for(i=0; i<a[0]; i++) printf("%c",sol[i]);
      }
   }
   return sol;
}

int clear_standard_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(509,a,b,c,0);

   return fail;
}

int clear_dobldobl_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(510,a,b,c,0);

   return fail;
}

int clear_quaddobl_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(511,a,b,c,0);

   return fail;
}

int clear_multprec_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(515,a,b,c,0);

   return fail;
}

int clear_varbprec_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc4c(519,a,b,c,0);

   return fail;
}
