/* This file "syspool.c" contains the definitions of the operations
 * declared in the file "syspool.h". */

#include "syspool.h"

int syspool_standard_initialize ( int n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(300,&n,b,c);
   return fail;
}

int syspool_dobldobl_initialize ( int n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(318,&n,b,c);
   return fail;
}

int syspool_quaddobl_initialize ( int n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(319,&n,b,c);
   return fail;
}

int syspool_standard_size ( int *n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(301,n,b,c);
   return fail;
}

int syspool_dobldobl_size ( int *n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(316,n,b,c);
   return fail;
}

int syspool_quaddobl_size ( int *n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(317,n,b,c);
   return fail;
}

int syspool_standard_read_system ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(302,&k,b,c);
   return fail;
}

int syspool_standard_write_system ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(303,&k,b,c);
   return fail;
}

int syspool_standard_create ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(304,&k,b,c);
   return fail;
}

int syspool_dobldobl_create ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(608,&k,b,c);
   return fail;
}

int syspool_quaddobl_create ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(609,&k,b,c);
   return fail;
}

int syspool_standard_refiner ( int k, int n, int *m, double *c )
{
   int fail,a[2];
   a[0] = k;
   a[1] = n;
   fail = _ada_use_c2phc(305,a,m,c);
   return fail;
}

int syspool_copy_to_standard_container ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(313,&k,b,c);

   return fail;
}

int syspool_copy_to_dobldobl_container ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(314,&k,b,c);

   return fail;
}

int syspool_copy_to_quaddobl_container ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(315,&k,b,c);

   return fail;
}

int syspool_standard_clear ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(697,a,b,c);
}

int syspool_dobldobl_clear ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(698,a,b,c);
}

int syspool_quaddobl_clear ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(699,a,b,c);
}
