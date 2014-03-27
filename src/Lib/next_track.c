/* The file next_track.c defines the functions of next_track.h. */

extern void adainit ( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal ( void );

int initialize_standard_homotopy ( int fixed_gamma )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(500,&fixed_gamma,b,c);

   return fail;
}

int initialize_dobldobl_homotopy ( int fixed_gamma )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(501,&fixed_gamma,b,c);

   return fail;
}

int initialize_quaddobl_homotopy ( int fixed_gamma )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(502,&fixed_gamma,b,c);

   return fail;
}

int initialize_multprec_homotopy ( int fixed_gamma, int decimals )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(512,&fixed_gamma,&decimals,c);

   return fail;
}

int initialize_standard_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(503,&k,b,c);

   return fail;
}

int initialize_dobldobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(504,&k,b,c);

   return fail;
}

int initialize_quaddobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(505,&k,b,c);

   return fail;
}

int initialize_multprec_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(513,&k,b,c);

   return fail;
}

int next_standard_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(506,&k,b,c);

   return fail;
}

int next_dobldobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(507,&k,b,c);

   return fail;
}

int next_quaddobl_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(508,&k,b,c);

   return fail;
}

int next_multprec_solution ( int k )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(514,&k,b,c);

   return fail;
}

int clear_standard_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(509,a,b,c);

   return fail;
}

int clear_dobldobl_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(510,a,b,c);

   return fail;
}

int clear_quaddobl_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(511,a,b,c);

   return fail;
}

int clear_multprec_tracker ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(515,a,b,c);

   return fail;
}
