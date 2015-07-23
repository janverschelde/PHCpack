/* The file scalers.c contains the definitions of the functions
 * declared in scalers.h. */

#include "scalers.h"

int standard_scale_system ( int mode, double *cff )
{
   int fail,*b;

   fail = _ada_use_c2phc(590,&mode,b,cff);

   return fail;
}

int dobldobl_scale_system ( int mode, double *cff )
{
   int fail,*b;

   fail = _ada_use_c2phc(591,&mode,b,cff);

   return fail;
}

int quaddobl_scale_system ( int mode, double *cff )
{
   int fail,*b;

   fail = _ada_use_c2phc(592,&mode,b,cff);

   return fail;
}

int standard_scale_solutions ( int dim, int basis, double *cff )
{
   int fail;

   fail = _ada_use_c2phc(594,&dim,&basis,cff);

   return fail;
}

int dobldobl_scale_solutions ( int dim, int basis, double *cff )
{
   int fail;

   fail = _ada_use_c2phc(595,&dim,&basis,cff);

   return fail;
}

int quaddobl_scale_solutions ( int dim, int basis, double *cff )
{
   int fail;

   fail = _ada_use_c2phc(596,&dim,&basis,cff);

   return fail;
}
