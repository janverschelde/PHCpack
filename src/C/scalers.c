/* The file scalers.c contains the definitions of the functions
 * declared in scalers.h. */

#include "scalers.h"

int standard_scale_system ( int mode, double *cff )
{
   int fail,*b;

   fail = _ada_use_c2phc4c(590,&mode,b,cff,0);

   return fail;
}

int dobldobl_scale_system ( int mode, double *cff )
{
   int fail,*b;

   fail = _ada_use_c2phc4c(591,&mode,b,cff,0);

   return fail;
}

int quaddobl_scale_system ( int mode, double *cff )
{
   int fail,*b;

   fail = _ada_use_c2phc4c(592,&mode,b,cff,0);

   return fail;
}

int standard_scale_solutions ( int dim, int basis, double *cff )
{
   int fail;

   fail = _ada_use_c2phc4c(594,&dim,&basis,cff,0);

   return fail;
}

int dobldobl_scale_solutions ( int dim, int basis, double *cff )
{
   int fail;

   fail = _ada_use_c2phc4c(595,&dim,&basis,cff,0);

   return fail;
}

int quaddobl_scale_solutions ( int dim, int basis, double *cff )
{
   int fail;

   fail = _ada_use_c2phc4c(596,&dim,&basis,cff,0);

   return fail;
}
