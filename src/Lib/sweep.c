/* The file sweep.c contains the definitions of the functions
 * declared in sweep.h. */

#include "sweep.h"

int sweep_define_parameters_numerically ( int nq, int nv, int np, int *pars )
{
   int fail = 0;
   int nqvp[3];
   double *c;

   nqvp[0] = nq;
   nqvp[1] = nv;
   nqvp[2] = np;
   fail = _ada_use_c2phc(610,nqvp,pars,c);

   return fail;
}

int sweep_get_number_of_equations ( int *n )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(612,n,b,c);

   return fail;
}

int sweep_get_number_of_variables ( int *n )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(613,n,b,c);

   return fail;
}

int sweep_get_number_of_parameters ( int *n )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(614,n,b,c);

   return fail;
}

int sweep_get_indices_numerically ( int *idx )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(615,idx,b,c);

   return fail;
}

int sweep_clear_definitions ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(617,a,b,c);

   return fail;
}
