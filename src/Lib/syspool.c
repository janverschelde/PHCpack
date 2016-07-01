/* This file "syspool.c" contains the definitions of the operations
 * declared in the file "syspool.h". */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

int syspool_initialize ( int n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(300,&n,b,c);
   return fail;
}

int syspool_size ( int *n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(301,n,b,c);
   return fail;
}

int syspool_read_system ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(302,&k,b,c);
   return fail;
}

int syspool_write_system ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(303,&k,b,c);
   return fail;
}

int syspool_create ( int k )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(304,&k,b,c);
   return fail;
}

int syspool_refine_root ( int k, int n, int *m, double *c )
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
