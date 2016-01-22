/* The file numbtrop.c contains the definitions of the functions with
 * prototypes documented in numbtrop.h. */

#include "numbtrop.h"

int standard_initialize
 ( int nbt, int dim, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[nbt*(dim+1)];

   nd[0] = nbt;
   nd[1] = dim;
   for(k=0; k<nbt; k++) cff[k] = dir[k];
   for(k=0; k<dim; k++) cff[nbt*dim+k] = err[k];

   fail = _ada_use_c2phc(711,nd,wnd,cff);

   return fail;
}

int standard_size ( int *nbt )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(720,nbt,b,c);

   return fail;
}

int standard_retrieve_tropism
 ( int dim, int idx, int *wnd, double *dir, double *err )
{
   int fail,k;
   int nd[2];
   double cff[dim+1];

   nd[0] = dim;
   nd[1] = idx;

   fail = _ada_use_c2phc(723,nd,wnd,cff);

   for(k = 0; k<dim; k++) dir[k] = cff[k];
   *err = cff[dim];

   return fail;
}

int standard_clear( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(726,a,b,c);

   return fail;
}
